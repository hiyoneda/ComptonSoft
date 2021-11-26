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

#ifndef COMPTONSOFT_AllSkyPrimaryGenPowerLaw_H
#define COMPTONSOFT_AllSkyPrimaryGenPowerLaw_H 1

#include <boost/multi_array.hpp>
#include "BasicPrimaryGen.hh"
#include "IsotropicPrimaryGen.hh"
#include "G4ThreeVector.hh"
#include "AstroUnits.hh"

#include <healpix_cxx/healpix_base.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <healpix_cxx/fitshandle.h>

namespace comptonsoft {


/**
 * ComptonSoft PrimaryGen module.
 * Imaging polarimetry realized.
 *
 * @author Naomi Tsuji, 
 * @date 2021-10-22
 *
 */

using image_t = boost::multi_array<double, 2>;

class AllSkyPrimaryGenPowerLaw : public anlgeant4::IsotropicPrimaryGen
{
  DEFINE_ANL_MODULE(AllSkyPrimaryGenPowerLaw, 1.0);
public:
  AllSkyPrimaryGenPowerLaw();
  ~AllSkyPrimaryGenPowerLaw();

  anlnext::ANLStatus mod_define() override;
  anlnext::ANLStatus mod_initialize() override;
  anlnext::ANLStatus mod_end_run() override;

  void makePrimarySetting() override;

protected:
  void inputImage(std::string filename, Healpix_Map<double>& image , int colnum, int hdunum , anlnext::ANLStatus* status);
  void inputImage(fitshandle* inp, Healpix_Map<double>& image, int colnum, int hdunum , anlnext::ANLStatus* status);
  void setCoordinate(anlnext::ANLStatus* status);
  void buildPixelIntegral();
  int samplePixel();

/**
 * rotate v as uz will be a new Z-axis.
 * this function is an inverse of G4ThreeVector::rotateUz(G4ThreeVector),
 * which rotates v as the current Z-axis into a uz.
 * The rotation matrix is given as the transpose of G4ThreeVector::rotateUz().
 */
  void rotateCoordinateZ(G4ThreeVector& v, G4ThreeVector uz);

private:
  std::string fitsFilename_;

  int	hdunum_ = 2; //hdu number to read healpix

  int	colnumNorm_ = 1;
  int	colnumIndex_ = 2;
  int	colnumEnergyRef_ = 3;
  Healpix_Map<double> imageNorm_; // normalization map
  Healpix_Map<double> imageIndex_; // index map
  Healpix_Map<double> imageEnergyRef_; // reference energy map in a unit of MeV.

  std::vector<double> imageIntegratedPhotonFlux_; // flux map, calculated here.
  std::vector<double> imageIntegratedEnergyFlux_; // flux map, calculated here.

  std::vector<double> pixelIntegral_;

  // -- for debug --
  bool useFluxForTest_ = true;
  int colnumFluxForTest_ = 4;
  Healpix_Map<double> imageIntegratedFluxForTest_; // flux map to test the flux integration.
  // ---------------

  int npix_;
  int nside_;
  std::vector<double> imageRA_; // l in Galactic coordinate.
  std::vector<double> imageDec_; // b in Galactic coordinate.

  double detectorRollAngle_ = 0.0;
  bool setPolarization_ = false;
  bool monoPhotonIndex_ = false;

  double energyMin_;
  double energyMax_;
  double photonIndex_ = 2.0;

};

} /* namespace comptonsoft */

#endif /* COMPTONSOFT_AllSkyPrimaryGenPowerLaw_H */
