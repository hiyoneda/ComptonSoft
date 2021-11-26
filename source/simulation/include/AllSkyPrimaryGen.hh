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

#ifndef COMPTONSOFT_AllSkyPrimaryGen_H
#define COMPTONSOFT_AllSkyPrimaryGen_H 1

#include <boost/multi_array.hpp>
#include <string>
#include "BasicPrimaryGen.hh"
#include "IsotropicPrimaryGen.hh"
#include "AllSkyPrimaryGenPowerLaw.hh"
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

class AllSkyPrimaryGen : public comptonsoft::AllSkyPrimaryGenPowerLaw
{
  DEFINE_ANL_MODULE(AllSkyPrimaryGen, 1.0);
public:
  AllSkyPrimaryGen();
  ~AllSkyPrimaryGen();

  anlnext::ANLStatus mod_define() override;
  anlnext::ANLStatus mod_initialize() override;
  anlnext::ANLStatus mod_end_run() override;

  void makePrimarySetting() override;

protected:
  void setCoordinate(anlnext::ANLStatus* status);
  void buildEnergyPixelIntegral();
  int sampleEnergyBandIndex();
  int samplePixel(int energyBandIndex);
  void calcIntegratedPhotonFluxInEnergyBand( Healpix_Map<double> photonfluxmap1, Healpix_Map<double> photonfluxmap2, double e1, double e2, int imap);

private:
  std::string fitsFilename_;

  int nmap_;

  std::vector< Healpix_Map<double> > differentialPhotonFluxMap_; // normalization map
  std::vector< double > energyBand_; // normalization map

  std::vector< std::vector<double> > indexMap_;
  std::vector< std::vector<double> > integratedPhotonFluxMap_;

  std::vector<double> energyIntegral_;
  std::vector< std::vector<double> > pixelIntegral_;

  int hdunum_ = 2; //hdu number to read healpix
  int npix_;
  int nside_;

  std::vector<double> imageRA_; // l in Galactic coordinate.
  std::vector<double> imageDec_; // b in Galactic coordinate.

  double detectorRollAngle_ = 0.0;
  bool setPolarization_ = false;

  double energyMin_;
  double energyMax_;

};

} /* namespace comptonsoft */

#endif /* COMPTONSOFT_AllSkyPrimaryGen_H */
