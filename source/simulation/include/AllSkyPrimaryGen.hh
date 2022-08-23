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

#include <utility>
#include <memory>
#include <boost/multi_array.hpp>
#include "BasicPrimaryGen.hh"
#include "IsotropicPrimaryGen.hh"
#include "G4ThreeVector.hh"

#include <healpix_cxx/healpix_base.h>
#include <healpix_cxx/healpix_map.h>

namespace comptonsoft {


/**
 * ComptonSoft PrimaryGen module for all-sky input via HEALPix.
 *
 * @author Naomi Tsuji, Hiroki Yoneda, Hirokazu Odaka
 * @date 2021-10-22
 * @date 2022-08-11
 *
 */

using image_t = boost::multi_array<double, 2>;

class AllSkyPrimaryGen : public anlgeant4::IsotropicPrimaryGen
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
  loadMultiBandImages();
  constructMaps(anlnext::ANLStatus& status);
  
private:
  /* module parameters */
  std::string filename_;
  int hdu_index_map_;
  int hdu_index_energy_;
  double origin_longitude_;
  double origin_latitude_;
  double pole_longitude_;
  double pole_latitude_;

  /* input from FITS file */
  /**
     map_mode_ is given by FITS header keyword MAPMODE
     1: monochromatic, 2: multiband-images, 3: power law
  */
  std::string map_mode_;

  /**
     num_maps_ is given by FITS header keyword NMAP
  */
  int num_maps_;

  std::vector<Healpix_Map<double>> maps_;
  std::vector<double> energies_;

  /* internal class members */
  
  struct band_intensity {
    double emin;
    double emax;
    double photon_index;
    double intensity;
  };
  std::vector<Healpix_Map<band_intensity>> band_maps_;
  std::vector<double> integrals_;
};

} /* namespace comptonsoft */

#endif /* COMPTONSOFT_AllSkyPrimaryGen_H */
