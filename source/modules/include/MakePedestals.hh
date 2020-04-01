/*************************************************************************
 *                                                                       *
 * Copyright (c) 2019 Hirokazu Odaka                                     *
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

/**
 * MakePedestals.
 *
 * @author Hirokazu Odaka & Tsubasa Tamba
 * @date 2019-05
 * @date 2019-07-18 | merged to comptonsoft
 * @date 2020-04-01 | v1.1
 *
 */

#ifndef COMPTONSOFT_MakePedestals_H
#define COMPTONSOFT_MakePedestals_H 1

#include "VCSModule.hh"

namespace comptonsoft{

class MakePedestals : public VCSModule
{
  DEFINE_ANL_MODULE(MakePedestals, 1.1);
  // ENABLE_PARALLEL_RUN();
public:
  MakePedestals();
  
public:
  anlnext::ANLStatus mod_define() override;
  anlnext::ANLStatus mod_end_run() override;

private:
  std::string filename_;
  int detector_id_ = 0;
};

} /* namespace comptonsoft */

#endif /* COMPTONSOFT_MakePedestals_H */
