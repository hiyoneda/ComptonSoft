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

#include "HistogramARM.hh"

#include "TH1.h"
#include "AstroUnits.hh"
#include "EventReconstruction.hh"
#include "BasicComptonEvent.hh"
#include "DetectorGroupManager.hh"

using namespace anl;

namespace comptonsoft
{

HistogramARM::HistogramARM()
  : detectorGroupManager_(nullptr),
    eventReconstruction_(nullptr),
    numBins_(500), range0_(-25.0), range1_(+25.0)
{
}

ANLStatus HistogramARM::mod_startup()
{
  register_parameter(&numBins_, "number_of_bins");
  register_parameter(&range0_, "range_min", 1.0, "degree");
  register_parameter(&range1_, "range_max", 1.0, "degree");
  return AS_OK;
}

ANLStatus HistogramARM::mod_his()
{
  GetANLModule("DetectorGroupManager", &detectorGroupManager_);
  GetANLModule("EventReconstruction", &eventReconstruction_);
  
  VCSModule::mod_his();
  mkdir();
  
  hist_all_ = new TH1D("arm_all", "ARM (All)",
                       numBins_, range0_, range1_);

  const std::vector<HitPattern>& hitPatterns
    = detectorGroupManager_->getHitPatterns();
  const std::size_t numHitPatterns = hitPatterns.size();
  hist_vec_.resize(numHitPatterns);
  for (std::size_t i=0; i<numHitPatterns; i++) {
    std::string histName = "arm_";
    std::string histTitle = "ARM (";
    histName += hitPatterns[i].ShortName();
    histTitle += hitPatterns[i].Name();
    histTitle += ")";
    hist_vec_[i] = new TH1D(histName.c_str(), histTitle.c_str(),
                            numBins_, range0_, range1_);
  }

  return AS_OK;
}

ANLStatus HistogramARM::mod_ana()
{
  if (!Evs("EventReconstruction:OK")) {
    return AS_OK;
  }
  
  const BasicComptonEvent& event = eventReconstruction_->getComptonEvent();
  const double cosThetaE = event.CosThetaE();
  if (cosThetaE < -1.0 || +1.0 < cosThetaE) {
    return AS_OK;
  }

  const double arm = event.DeltaTheta()/degree;

  hist_all_->Fill(arm);

  for (std::size_t i=0; i<hist_vec_.size(); i++) {
    if (eventReconstruction_->HitPatternFlag(i)) {
      hist_vec_[i]->Fill(arm);
    }
  }

  return AS_OK;
}

} /* namespace comptonsoft */
