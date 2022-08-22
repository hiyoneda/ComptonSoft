/*************************************************************************
 *                                                                       *
 * Copyright (c) 2011 Hirokazu Odaka                                     *
 *                                                                       *
 * This program is free software: you can redistribute it and/or modify  *
 * it under the terms of the GNU General Public License as published by  *
 * the Free Software Foundation, either version 3 of the License, or     *
 * (at your option) any later version.                                   *
 *                                                                       *
 * This program is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 * GNU General Public License for more details.                          *
 *                                                                       *
 * You should have received a copy of the GNU General Public License     *
 * along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 *                                                                       *
 *************************************************************************/

#include "AllSkyPrimaryGen.hh"
#include "Randomize.hh"
#include "AstroUnits.hh"

#include <healpix_cxx/healpix_map_fitsio.h>
#include <healpix_cxx/fitshandle.h>

namespace
{

/**
 * rotate v as uz will be a new Z-axis.
 * this function is an inverse of G4ThreeVector::rotateUz(G4ThreeVector),
 * which rotates v as the current Z-axis into a uz.
 * The rotation matrix is given as the transpose of G4ThreeVector::rotateUz().
 */
void rotateCoordinateZ(G4ThreeVector& v, G4ThreeVector uz);

} /* anonymous namespace */


using namespace anlnext;
namespace unit = anlgeant4::unit;

namespace comptonsoft
{

AllSkyPrimaryGen::AllSkyPrimaryGen()
  : filename_("allsky.fits"),
    hdu_index_(1),
    origin_longitude_(0.0),
    origin_latitude_(0.0),
    pole_longitude_(0.0),
    pole_latitude_(90.0*unit::degree),
    source_flux_(1.0)
{
  add_alias("AllSkyPrimaryGen");
}

AllSkyPrimaryGen::~AllSkyPrimaryGen() = default;

ANLStatus AllSkyPrimaryGen::mod_define()
{
  ANLStatus status = IsotropicPrimaryGen::mod_define();
  if (status != AS_OK) {
    return status;
  }
 
  define_parameter("filename", &mod_class::filename_);
  define_parameter("hdu_index", &mod_class::hdu_index_);
  define_parameter("origin_longitude_", &mod_class::origin_longitude_);
  define_parameter("origin_latitude_", &mod_class::origin_latitude_);
  define_parameter("pole_longitude_", &mod_class::pole_longitude_);
  define_parameter("pole_latitude_", &mod_class::pole_latitude_);

  return AS_OK;
}

ANLStatus AllSkyPrimaryGen::mod_initialize()
{
  ANLStatus status = IsotropicPrimaryGen::mod_initialize();
  if (status != AS_OK) {
    return status;
  }

  std::unique_ptr<fitshandle> fits(new fitshandle);
  fits->open(filename_);
  fits->goto_hdu(hdu_index_);
  fits->get_key("MAP_MODE", map_mode_);
  fits->get_key("NMAP", num_maps_);

  if (map_mode_ == "multi-band") {
    loadMultiBandImages(fits.get(), num_maps_, status);
  }
  else {
    std::cout << "Undefined map mode: " << map_mode_ << std::endl;
    status = AS_QUIT_ERROR;
  }

  fits->close();

  using std::cos;
  using std::sin;
  G4ThreeVector xaxis();
  G4ThreeVector zaxis();
  // check if x dot z=0
  G4ThreeVector yaxis();
  rotation_matrix_;

  constructMaps();

  return AS_OK;
}

void AllSkyPrimaryGen::makePrimarySetting()
{
  using std::cos;
  using std::sin;
  using std::sqrt;
  using std::atan2;

  std::pair<int, int> currentPixel = samplePixel();
  const int ix = currentPixel.first;
  const int iy = currentPixel.second;
  const double d = Distance();
  G4ThreeVector centerDirection = CenterDirection();
  const G4ThreeVector centerPosition = CenterPosition();
  const double radius = Radius();

  const double cosTheta = cos(-imageDec_[ix][iy] + 0.5*CLHEP::pi);
  const double sinTheta = sqrt(1.0-cosTheta*cosTheta);
  const double phi = imageRA_[ix][iy];
  G4ThreeVector v(d*sinTheta*cos(phi), d*sinTheta*sin(phi), d*cosTheta);
  G4ThreeVector initialPosition = v;
  rotateCoordinateZ(v, -centerDirection);

  G4ThreeVector v2 = v.orthogonal();
  v2.setMag( radius * sqrt(G4UniformRand()) );
  const G4double chi = CLHEP::twopi * G4UniformRand();
  v2.rotate(chi, v);
  G4ThreeVector position = centerPosition + v + v2;

  const double energy = sampleEnergy();
  G4ThreeVector direction = (-v).unit();


  G4ThreeVector uz(0.0, 0.0, 1.0);
  G4ThreeVector uy(0.0, 1.0, 0.0);
  G4ThreeVector newUz = uz;
  rotateCoordinateZ(newUz, -centerDirection);
  G4ThreeVector north(newUz.x(), newUz.y(), 0.0);
  north = north.unit();
  G4ThreeVector dety = north.rotate(detectorRollAngle_, -uz);
  const double detAngle = atan2(dety.y(), dety.x()) - 0.5*CLHEP::pi;
  position.rotate(-detAngle, uz);
  direction.rotate(-detAngle, uz);

  setPrimary(position, energy, direction);

#if 0
  setUnpolarized();
#endif
}

ANLStatus AllSkyPrimaryGen::mod_end_run()
{
  const double radius = Radius();
  const double area = CLHEP::pi*radius*radius;
  const double realTime = TotalEnergy()/(sourceFlux_*area);
  const double pflux = Number()/area/realTime;

  setRealTime(realTime);

  std::cout.setf(std::ios::scientific);
  std::cout << "AllSkyPrimaryGen::mod_end_run \n"
            << "  Number: " << Number() << "\n"
            << "  Flux: " << sourceFlux_/(unit::erg/unit::cm2/unit::s) << " erg/cm2/s\n"
            << "  Total Energy: " << TotalEnergy()/unit::keV << " keV = "
            << TotalEnergy()/unit::erg << " erg\n"
            << "  Area: " << area/unit::cm2 << " cm2\n"
            << "  Real time: " << realTime/unit::s << " s\n"
            << "  Photon flux: " << pflux/(1.0/unit::cm2/unit::s) << " photons/cm2/s\n"
            << std::endl;
  std::cout.unsetf(std::ios::scientific);

  return AS_OK;
}

void AllSkyPrimaryGen::loadMultiBandImages(fitshandle* fits, int num_maps_, ANLStatus& status)
{
  maps_.resize(num_maps_);
  for (int i=0; i<num_maps_; i++) {
    read_Healpix_map_from_fits(*fits, maps_[i], i+1);
  }
}


} /* namespace comptonsoft */

namespace
{

void rotateCoordinateZ(G4ThreeVector& v, G4ThreeVector uz)
{
  uz = uz.unit();
  const double u1 = uz.x();
  const double u2 = uz.y();
  const double u3 = uz.z();
  const double up2 = u1*u1 + u2*u2;

  if (up2>0) {
    const double up = std::sqrt(up2);
    const double vx = v.x();
    const double vy = v.y();
    const double vz = v.z();
    v.setX( u1*u3*vx/up + u2*u3*vy/up - up*vz);
    v.setY(-u2*vx/up    + u1*vy/up           );
    v.setZ( u1*vx       + u2*vy       + u3*vz);
  }
  else if (u3 < 0.) {
    v.setX(-v.x());
    v.setZ(-v.z());
  }
}

} /* anonymous namespace */
