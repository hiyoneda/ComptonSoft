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


namespace cfitsio
{
extern "C" {
#include "fitsio.h"
}
}

using namespace anlnext;
namespace unit = anlgeant4::unit;

namespace comptonsoft
{

AllSkyPrimaryGen::AllSkyPrimaryGen()
  : fitsFilename_("allskymap.fits")
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

  define_parameter("fits_filename", &mod_class::fitsFilename_);
  set_parameter_description("filename (.fits) that contains the allsky differential flux maps at certain energies.");

  define_parameter("hdunum", &mod_class::hdunum_ );
  set_parameter_description("HDU number in healpix file.");

  define_parameter("detector_roll_angle", &mod_class::detectorRollAngle_, unit::degree, "degree");
  set_parameter_description("Roll angle of the detector with respect to north direction.");
  define_parameter("set_polarization", &mod_class::setPolarization_);
  set_parameter_description("If true, polarization is on.");

  undefine_parameter("intensity");

  return AS_OK;
}

ANLStatus AllSkyPrimaryGen::mod_initialize()
{

  ANLStatus status = IsotropicPrimaryGen::mod_initialize();
  if (status != AS_OK) {
    return status;
  }

  energyMin_ = getEnergyMin();
  energyMax_ = getEnergyMax();

  fitshandle* inp = new fitshandle();
  inp->open(fitsFilename_);
  inp->goto_hdu(hdunum_);
    
  int nmap_raw_;
  inp->get_key("NMAP", nmap_raw_);
  nmap_ = 0;
  for(int imap = 0; imap < nmap_raw_; ++imap){
    double this_energy;
    inp->get_key("ENE"+std::to_string(imap+1), this_energy);
    
    if(energyMin_ <= this_energy && this_energy < energyMax_){
      Healpix_Map<double> this_differentialPhotonFluxMap_;
      std::cout << "  Diffrential photon flux map at " << std::to_string(this_energy) << " MeV with an unit of ph/cm2/MeV/s: " << std::endl ;
      std::cout << "    <- solid angle of a single pixel should be already considered in the normalization map " << std::endl ;
      AllSkyPrimaryGenPowerLaw::inputImage( inp, this_differentialPhotonFluxMap_, imap+1, hdunum_, &status );
      differentialPhotonFluxMap_.push_back(this_differentialPhotonFluxMap_);
      energyBand_.push_back(this_energy);
      ++nmap_;
    }
  }

  if(nmap_ < 2){
    status = AS_ERROR;
    return status;
  }
  npix_ = differentialPhotonFluxMap_[0].Npix();
  nside_ = differentialPhotonFluxMap_[0].Nside();
  
  indexMap_.resize(nmap_ - 1);
  integratedPhotonFluxMap_.resize(nmap_ - 1);
  for(int imap = 0; imap < nmap_ - 1; ++imap){
    calcIntegratedPhotonFluxInEnergyBand(differentialPhotonFluxMap_[imap], differentialPhotonFluxMap_[imap+1],
                                         energyBand_[imap], energyBand_[imap+1], imap);
  }

  setCoordinate(&status);
  if (status != AS_OK) {
    return status;
  }
    
  buildEnergyPixelIntegral();

  inp->close();

  return AS_OK;
}


void AllSkyPrimaryGen::makePrimarySetting()
{
  using std::cos;
  using std::sin;
  using std::sqrt;
  using std::atan2;

  const int eindex = sampleEnergyBandIndex();
  const int ix = samplePixel(eindex);

  const double d = Distance();
  G4ThreeVector centerDirection = CenterDirection();
  const G4ThreeVector centerPosition = CenterPosition();
  const double radius = Radius();

  const double cosTheta = cos(-imageDec_[ix] + 0.5*CLHEP::pi); //const double cosTheta = cos(-imageDec_[ix][iy] + 0.5*CLHEP::pi);
  const double sinTheta = sqrt(1.0-cosTheta*cosTheta);
  const double phi = imageRA_[ix]; //const double phi = imageRA_[ix][iy];
  G4ThreeVector v(d*sinTheta*cos(phi), d*sinTheta*sin(phi), d*cosTheta);
  G4ThreeVector initialPosition = v;
  rotateCoordinateZ(v, -centerDirection);

  G4ThreeVector v2 = v.orthogonal();
  v2.setMag( radius * sqrt(G4UniformRand()) );
  const G4double chi = CLHEP::twopi * G4UniformRand();
  v2.rotate(chi, v);
  G4ThreeVector position = centerPosition + v + v2;

  // set spectrum
  const double pindex = indexMap_[eindex][ix]; // photon index at the current pix
  const double energy = sampleFromPowerLaw( pindex, energyBand_[eindex], energyBand_[eindex+1] );

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

  if (!setPolarization_) {
    setPrimary(position, energy, direction);
    return;
  }

#if 0
  std::cout << "position: (" << position.x() << ", " << position.y() << ", " << position.z() << ")" << std::endl;
  std::cout << "direction: (" << direction.x() << ", " << direction.y() << ", " << direction.z() << ")" << std::endl;
#endif
}

ANLStatus AllSkyPrimaryGen::mod_end_run()
{
  double totalPhotonFlux_ = 0;
  for(int imap = 0; imap < nmap_ - 1; ++imap){
    totalPhotonFlux_ += std::accumulate(integratedPhotonFluxMap_[imap].begin(), integratedPhotonFluxMap_[imap].end(), 0.0) * (1.0/unit::cm2/unit::s);
  }
  const double radius = Radius();
  const double area = CLHEP::pi*radius*radius;
  const double realTime = Number()/(totalPhotonFlux_*area);
  const double pflux = Number()/area/realTime;

  setRealTime(realTime);

  std::cout.setf(std::ios::scientific);
  std::cout << "AllSkyPrimaryGen::mod_end_run \n"
            << "  Input file: " << fitsFilename_ << "\n"
            << "  Number (event number; photon number): " << Number() << "\n"
//            << "  Flux (total all-sky integrated) in " << energyMin_ <<" to " << energyMax_ << " MeV: " << totalEnergyFlux_/(unit::MeV/unit::cm2/unit::s) << " MeV/cm2/s = " << totalEnergyFlux_/(unit::erg/unit::cm2/unit::s) << " erg/cm2/s\n"
            << "  Total Energy (=sum of sampled photon energy): " << TotalEnergy()/unit::keV << " keV = "
            << TotalEnergy()/unit::erg << " erg\n"
            << "  Area: " << area/unit::cm2 << " cm2\n"
            << "  Real time : " << realTime/unit::s << " s\n"
            << "  Photon flux in " << energyBand_[0] << " to " << energyBand_[nmap_-1] << " MeV: " << pflux/(1.0/unit::cm2/unit::s) << " photons/cm2/s\n"
            << std::endl;
  std::cout.unsetf(std::ios::scientific);

  return AS_OK;
}

void AllSkyPrimaryGen::setCoordinate(ANLStatus* status)
{
  using std::sin;
  using std::cos;

  *status = AS_OK;

  double z_ ;
  double phi_ ;
  double theta_ ;
  double ra_ ; // or l
  double dec_; // or b

  imageRA_.resize(npix_);
  imageDec_.resize(npix_);

  for(int i = 0; i < npix_; ++i){
  	// ang = ( theta in radian, phi = zenith angle)
  	// z = cos(theta), phi = zenith
  	differentialPhotonFluxMap_[0].pix2zphi( i, z_, phi_ );
  	theta_ = acos( z_ );
  	phi_ = phi_ * 180./CLHEP::pi * unit::degree ;
  	theta_ = theta_ * 180./CLHEP::pi * unit::degree ;

    ra_ = phi_ ;
    dec_ = (theta_-90.*unit::degree )*(-1.0) ;

    // Insert (RA, Dec) or (l,b) to imageRA,Dec
  	imageRA_[i] = ra_;
    imageDec_[i] = dec_;
  	}
}

void AllSkyPrimaryGen::buildEnergyPixelIntegral()
{
  energyIntegral_.resize(nmap_);
  energyIntegral_[0] = 0.0;
  for(int imap = 0; imap < nmap_ - 1; ++imap){
    energyIntegral_[imap+1] = energyIntegral_[imap] + std::accumulate(integratedPhotonFluxMap_[imap].begin(), integratedPhotonFluxMap_[imap].end(), 0.0);
  }
  const double norm = energyIntegral_.back();
  for (auto& v: energyIntegral_) {
    v /= norm;
  }

  pixelIntegral_.resize(nmap_ - 1);
  for(int imap = 0; imap < nmap_ - 1; ++imap){
    std::vector<double> this_pixelIntegral_;
    this_pixelIntegral_.resize(npix_+1);
    this_pixelIntegral_[0] = 0.0;
    for(int ipix = 0; ipix < npix_; ++ipix){
      this_pixelIntegral_[ipix+1] = this_pixelIntegral_[ipix] + integratedPhotonFluxMap_[imap][ipix];
    }
    const double this_norm = this_pixelIntegral_.back();
    for (auto& v: this_pixelIntegral_) {
      v /= this_norm;
    }
    pixelIntegral_[imap] = this_pixelIntegral_;
  }
}

int AllSkyPrimaryGen::sampleEnergyBandIndex()
{
  const double r = G4UniformRand();
  const std::vector<double>::const_iterator it = std::upper_bound(energyIntegral_.begin(), energyIntegral_.end(), r);
  const int r0 = it - energyIntegral_.begin() - 1;
  return r0 ;
}

int AllSkyPrimaryGen::samplePixel(int energyBandIndex)
{
  const double r = G4UniformRand();
  const std::vector<double>::const_iterator it = std::upper_bound(pixelIntegral_[energyBandIndex].begin(), 
                                                                  pixelIntegral_[energyBandIndex].end(), r);
  const int r0 = it - pixelIntegral_[energyBandIndex].begin() - 1;
  return r0 ;
}

void AllSkyPrimaryGen::calcIntegratedPhotonFluxInEnergyBand( Healpix_Map<double> photonfluxmap1, Healpix_Map<double> photonfluxmap2, double e1, double e2, int imap)
{
  std::vector<double> this_indexMap_;
  std::vector<double> this_integratedPhotonFluxMap_;
  this_indexMap_.resize(npix_);
  this_integratedPhotonFluxMap_.resize(npix_);
  for(int ipix = 0; ipix < npix_; ++ipix){
    double n1 = photonfluxmap1[ipix];
    double n2 = photonfluxmap2[ipix];

    double photonIndex_ = -1 * log(n2/n1) / log(e2/e1);
    double intg_1 = n1 * e1 / (-photonIndex_+1);
    double intg_2 = n1 * pow ( 1 / e1, (-1.0)*photonIndex_) * pow (e2, -photonIndex_+1)/(-photonIndex_+1);
    double intg = intg_2 - intg_1;
    
    this_indexMap_[ipix] = photonIndex_;
    this_integratedPhotonFluxMap_[ipix] = intg;
  }
  indexMap_[imap] = this_indexMap_;
  integratedPhotonFluxMap_[imap] = this_integratedPhotonFluxMap_;
}

} /* namespace comptonsoft */
