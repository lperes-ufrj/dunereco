#ifndef HIGHESTENERGYRECOVERTEXSHOWERSELECTOR_H_SEEN
#define HIGHESTENERGYRECOVERTEXSHOWERSELECTOR_H_SEEN
//STL
#include <iostream>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "canvas/Persistency/Common/FindManyP.h"
//LArSoft
#include "larcore/Geometry/Geometry.h"
//#include "larreco/RecoAlg/ShowerEnergyAlg.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
//CUSTOM
#include "RecoShowerSelector.h"

namespace FDSelectionTools{
  class HighestEnergyRecoVertexShowerSelector : public RecoShowerSelector{
    public:
      explicit HighestEnergyRecoVertexShowerSelector(fhicl::ParameterSet const& ps) 
        :
        fShowerModuleLabel(ps.get< std::string> ("ModuleLabels.ShowerModuleLabel")),
        fPFParticleModuleLabel(ps.get< std::string> ("ModuleLabels.PFParticleModuleLabel")) {};
        //fShowerEnergyAlg(ps.get<fhicl::ParameterSet>("ShowerEnergyAlg")) {};

    private:
      art::Ptr<recob::Shower> SelectShower(art::Event const & evt) override;
      std::string fShowerModuleLabel;
      std::string fPFParticleModuleLabel;
      //shower::ShowerEnergyAlg fShowerEnergyAlg;

  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::HighestEnergyRecoVertexShowerSelector)
#endif
