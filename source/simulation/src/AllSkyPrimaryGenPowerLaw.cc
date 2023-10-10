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

#include "AllSkyPrimaryGenPowerLaw.hh"
#include "Randomize.hh"
#include "AstroUnits.hh"


namespace cfitsio
{
extern "C" {
#include "fitsio.h"
}
}

namespace
{
double calcIntegratedPhotonFlux( double norm_, double photonIndex_, double energy, double E_min, double E_max );
double calcIntegratedEnergyFlux( double norm_, double photonIndex_, double energy, double E_min, double E_max );
} /* anonymous namespace */

using namespace anlnext;
namespace unit = anlgeant4::unit;

namespace comptonsoft
{

AllSkyPrimaryGenPowerLaw::AllSkyPrimaryGenPowerLaw()
  : fitsFilename_("allskymap_powerlaw.fits")
{
  add_alias("AllSkyPrimaryGenPowerLaw");
}

AllSkyPrimaryGenPowerLaw::~AllSkyPrimaryGenPowerLaw() = default;

ANLStatus AllSkyPrimaryGenPowerLaw::mod_define()
{

  ANLStatus status = IsotropicPrimaryGen::mod_define();

  if (status != AS_OK) {
    return status;
  }

  define_parameter("fits_filename", &mod_class::fitsFilename_);
  set_parameter_description("filename (.fits) that contains the allsky maps of the normalization and powerlaw index.");

  define_parameter("hdunum", &mod_class::hdunum_ );
  set_parameter_description("HDU number in healpix file.");

  define_parameter("colnum_norm", &mod_class::colnumNorm_);
  set_parameter_description("Column number of the normalization.");
  define_parameter("colnum_index", &mod_class::colnumIndex_);
  set_parameter_description("Column number of the index.");
  define_parameter("colnum_energyref", &mod_class::colnumEnergyRef_);
  set_parameter_description("Column number of the reference energy.");

  define_parameter("detector_roll_angle", &mod_class::detectorRollAngle_, unit::degree, "degree");
  set_parameter_description("Roll angle of the detector with respect to north direction.");
  define_parameter("set_polarization", &mod_class::setPolarization_);
  set_parameter_description("If true, polarization is on.");

  define_parameter("mono_photon_index", &mod_class::monoPhotonIndex_);
  set_parameter_description("If true, index is constant over the map.");

  undefine_parameter("intensity");

  return AS_OK;
}

ANLStatus AllSkyPrimaryGenPowerLaw::mod_initialize()
{

  ANLStatus status = IsotropicPrimaryGen::mod_initialize();
  if (status != AS_OK) {
    return status;
  }

  energyMin_ = getEnergyMin();
  energyMax_ = getEnergyMax();
  photonIndex_ = getPhotonIndex();

  fitshandle* inp = new fitshandle();
  inp->open(fitsFilename_);
  inp->goto_hdu(hdunum_);

  // Read normalization map
  std::cout << "  Normalization map: " << std::endl ;
  std::cout << "    <- differential flux map at the reference energy with an unit of ph/cm2/MeV/s " << std::endl ;
  std::cout << "    <- solid angle of a single pixel should be already considered in the normalization map " << std::endl ;
  inputImage( inp, imageNorm_, colnumNorm_, hdunum_, &status );
  if (status != AS_OK) {
    return status;
  }
  npix_ = imageNorm_.Npix();
  nside_ = imageNorm_.Nside();

  // Read photon index and flux maps
  if(monoPhotonIndex_){
    imageIndex_.SetNside( nside_, RING );
    for(int ipix = 0; ipix < npix_; ++ipix){
      imageIndex_[ipix] = photonIndex_;
    }
  }else{
    std::cout << "  Index map: " << std::endl ;
    std::cout << "    <- photon index of the powerlaw distribution at each pixel " << std::endl ;
    inputImage( inp, imageIndex_, colnumIndex_, hdunum_ , &status ) ;
  }
  std::cout << "  Reference energy map: " << std::endl ;
  std::cout << "    <- its unit should be MeV " << std::endl ;
  inputImage( inp, imageEnergyRef_, colnumEnergyRef_, hdunum_ , &status ) ;

  if(useFluxForTest_){
    std::cout << "  Flux map (this is not used for the calculation): " << std::endl ;
    inputImage( inp, imageIntegratedFluxForTest_, colnumFluxForTest_, hdunum_ , &status ) ;
  }

  imageIntegratedPhotonFlux_.resize( npix_ );
  imageIntegratedEnergyFlux_.resize( npix_ );
  for(int ipix = 0; ipix < npix_; ++ipix){
    imageIntegratedPhotonFlux_[ipix] = calcIntegratedPhotonFlux( imageNorm_[ipix], imageIndex_[ipix], imageEnergyRef_[ipix] , energyMin_/unit::MeV, energyMax_/unit::MeV ) ;
    imageIntegratedEnergyFlux_[ipix] = calcIntegratedEnergyFlux( imageNorm_[ipix], imageIndex_[ipix], imageEnergyRef_[ipix] , energyMin_/unit::MeV, energyMax_/unit::MeV ) ;
  }
  // Print data
  if(useFluxForTest_){
    std::cout << "  Map data ( norm, photon index, flux (test), calculated photon flux, calculated energy flux, reference energy ): " << std::endl ;
    for(int ipix = 0; ipix < 10; ++ipix){
        std::cout << "    " << ipix << "\t" << imageNorm_[ipix] << ", " << imageIndex_[ipix] << ", " << imageIntegratedFluxForTest_[ipix] << ", " << imageIntegratedPhotonFlux_[ipix] << ", " << imageIntegratedEnergyFlux_[ipix] << ", " << imageEnergyRef_[ipix] << std::endl; 
    }
  }

  setCoordinate(&status);
  if (status != AS_OK) {
    return status;
  }

  buildPixelIntegral();

  inp->close();

  return AS_OK;
}


void AllSkyPrimaryGenPowerLaw::makePrimarySetting()
{
  using std::cos;
  using std::sin;
  using std::sqrt;
  using std::atan2;

  const int ix = samplePixel();

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
  const double pindex = imageIndex_[ix]; // photon index at the current pix
  const double energy = sampleFromPowerLaw( pindex, energyMin_, energyMax_ ) ;

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

ANLStatus AllSkyPrimaryGenPowerLaw::mod_end_run()
{
  const double totalEnergyFlux_ = std::accumulate(imageIntegratedEnergyFlux_.begin(), imageIntegratedEnergyFlux_.end(), 0.0) * (unit::MeV/unit::cm2/unit::s) ;
  const double totalPhotonFlux_ = std::accumulate(imageIntegratedPhotonFlux_.begin(), imageIntegratedPhotonFlux_.end(), 0.0) * (1.0/unit::cm2/unit::s) ;
  const double radius = Radius();
  const double area = CLHEP::pi*radius*radius;
  const double realTime_f = TotalEnergy()/(totalEnergyFlux_*area);
  const double realTime_p = Number()/(totalPhotonFlux_*area);
  const double pflux = Number()/area/realTime_p;

  setRealTime(realTime_p);

  std::cout.setf(std::ios::scientific);
  std::cout << "AllSkyPrimaryGenPowerLaw::mod_end_run \n"
            << "  Input file: " << fitsFilename_ << "\n"
            << "  Number (event number; photon number): " << Number() << "\n"
            << "  Flux (total all-sky integrated) in " << energyMin_ <<" to " << energyMax_ << " MeV: " << totalEnergyFlux_/(unit::MeV/unit::cm2/unit::s) << " MeV/cm2/s = " << totalEnergyFlux_/(unit::erg/unit::cm2/unit::s) << " erg/cm2/s\n"
            << "  Total Energy (=sum of sampled photon energy): " << TotalEnergy()/unit::keV << " keV = "
            << TotalEnergy()/unit::erg << " erg\n"
            << "  Area: " << area/unit::cm2 << " cm2\n"
            << "  Real time (flux): " << realTime_f/unit::s << " s\n"
            << "  Real time (photon): " << realTime_p/unit::s << " s\n"
            << "  Photon flux: " << pflux/(1.0/unit::cm2/unit::s) << " photons/cm2/s\n"
            << std::endl;
  std::cout.unsetf(std::ios::scientific);

  return AS_OK;
}

void AllSkyPrimaryGenPowerLaw::inputImage(std::string filename, Healpix_Map<double>& image , int colnum, int hdunum , ANLStatus* status )
{
  /*
  read_Healpix_map_from_fits	(	const std::string & 	filename,
  Healpix_Map< T > & 	map,
  int 	colnum = 1,
  int 	hdunum = 2
  )
  Opens the FITS file filename, jumps to the HDU hdunum, and reads the column colnum into map.
  */
  *status = AS_OK;
  std::cout<< "    FITS filename = " << filename << std::endl ;
  std::cout<< "    column = " << colnum << std::endl ;
  std::cout<< "    hdunum = " << hdunum << std::endl ;
  read_Healpix_map_from_fits( filename, image, colnum, hdunum ) ;
  std::cout<< "    Nside = " << image.Nside() << std::endl ;
  std::cout<< "    Npix  = " << image.Npix() << std::endl ;
}

void AllSkyPrimaryGenPowerLaw::inputImage(fitshandle* inp, Healpix_Map<double>& image, int colnum, int hdunum , anlnext::ANLStatus* status)
{
  /*
  * void read_Healpix_map_from_fits	(	fitshandle & 	inp,
  Healpix_Map< T > & 	map,
  int 	colnum = 1	 
  )			[inline]
  Reads the map stored in column colnum of the FITS binary table pointed to by inp into map.
  */
  *status = AS_OK;
  inp->goto_hdu(hdunum);
  read_Healpix_map_from_fits( *inp, image, colnum ) ;
  std::cout<< "    colname = " << inp->colname(colnum) << std::endl ;
  std::cout<< "    unit = " << inp->colunit(colnum) << std::endl ;
  std::cout<< "    column = " << colnum << std::endl ;
  std::cout<< "    hdunum = " << hdunum << std::endl ;
  std::cout<< "    Nside = " << image.Nside() << std::endl ;
  std::cout<< "    Npix  = " << image.Npix() << std::endl ;
}

void AllSkyPrimaryGenPowerLaw::setCoordinate(ANLStatus* status)
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
  	imageNorm_.pix2zphi( i, z_, phi_ );
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

void AllSkyPrimaryGenPowerLaw::buildPixelIntegral()
{
  pixelIntegral_.resize(npix_+1);
  pixelIntegral_[0] = 0.0;
  for(int ipix = 0; ipix < npix_; ++ipix){
    pixelIntegral_[ipix+1] = pixelIntegral_[ipix] + imageIntegratedPhotonFlux_[ipix];
  }
  const double norm = pixelIntegral_.back();
  for (auto& v: pixelIntegral_) {
    v /= norm;
  }
}

int AllSkyPrimaryGenPowerLaw::samplePixel()
{
  const double r = G4UniformRand();
  const std::vector<double>::const_iterator it = std::upper_bound(pixelIntegral_.begin(), pixelIntegral_.end(), r);
  const int r0 = it - pixelIntegral_.begin() - 1;
  return r0 ;
}

void AllSkyPrimaryGenPowerLaw::rotateCoordinateZ(G4ThreeVector& v, G4ThreeVector uz)
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

} /* namespace comptonsoft */

namespace
{

double calcIntegratedEnergyFlux( double norm_, double photonIndex_, double energy, double E_min, double E_max ){
  double intg_1 ;
  double intg_2 ;
  double intg_ ;
  if ( (1.999 < photonIndex_) && (photonIndex_ < 2.001) ){
    intg_1 = norm_ * pow ( 1 / energy , (-1.0)*photonIndex_) * log( E_min );
    intg_2 = norm_ * pow ( 1 / energy , (-1.0)*photonIndex_) * log( E_max );
  }
  else{
    intg_1 = norm_ * pow ( 1 / energy , (-1.0)*photonIndex_) * pow (E_min, -photonIndex_+2)/(-photonIndex_+2) ;
    intg_2 = norm_ * pow ( 1 / energy , (-1.0)*photonIndex_) * pow (E_max, -photonIndex_+2)/(-photonIndex_+2) ;
  }
  intg_ = intg_2 - intg_1 ;
  return intg_ ; // MeV/cm2/s/pix
}

double calcIntegratedPhotonFlux( double norm_, double photonIndex_, double energy, double E_min, double E_max ){
  double intg_1 ;
  double intg_2 ;
  double intg_ ;
  if ( (0.999 < photonIndex_) && (photonIndex_ < 1.001) ){
    intg_1 = norm_ * pow ( 1 / energy , (-1.0)*photonIndex_) * log( E_min );
    intg_2 = norm_ * pow ( 1 / energy , (-1.0)*photonIndex_) * log( E_max );
  }
  else{
    intg_1 = norm_ * pow ( 1 / energy , (-1.0)*photonIndex_) * pow (E_min, -photonIndex_+1)/(-photonIndex_+1) ;
    intg_2 = norm_ * pow ( 1 / energy , (-1.0)*photonIndex_) * pow (E_max, -photonIndex_+1)/(-photonIndex_+1) ;
  }
  intg_ = intg_2 - intg_1 ;
  return intg_ ; // ph/cm2/s/pix
}

} /* anonymous namespace */
