#ifndef PAND_SELECT_PARAMS_H
#define PAND_SELECT_PARAMS_H

#include "dunereco/FDSensOpt/FDSensOptData/EnergyRecoOutput.h"

namespace pandselect
{
  class PandSelectParams
  {
  public:

      double selTrackPandizzleScore;
      double selShowerPandrizzleScore;

      dune::EnergyRecoOutput energyRecoNumu;
      dune::EnergyRecoOutput energyRecoNue;
  };
}

#endif
