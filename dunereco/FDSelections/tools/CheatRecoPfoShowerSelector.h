#ifndef CHEATRECOPFOSHOWERSELECTOR_H_SEEN
#define CHEATRECOPFOSHOWERSELECTOR_H_SEEN

//STL
#include <iostream>
//ROOT
//ART
#include "art/Utilities/ToolMacros.h" 
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"

//LArSoft
#include "larcore/Geometry/Geometry.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

//DUNE
#include "dunereco/FDSelections/pandrizzle/PandrizzleAlg.h"

//CUSTOM
#include "RecoShowerSelector.h"

namespace FDSelectionTools{
  class CheatRecoPfoShowerSelector : public RecoShowerSelector{
    public:
      explicit CheatRecoPfoShowerSelector(fhicl::ParameterSet const& ps);

    private:
      art::Ptr<recob::Shower> SelectShower(art::Event const & evt) override;

      std::string fNuGenModuleLabel;
      std::string fAllShowerModuleLabel;
      std::string fFittedShowerModuleLabel;
      std::string fPFParticleModuleLabel;
      FDSelection::PandrizzleAlg fPandrizzleAlg;

      // Pandrizzle Stuff
      double fRecoNuVtxX;
      double fRecoNuVtxY;
      double fRecoNuVtxZ;
  };
}

DEFINE_ART_CLASS_TOOL(FDSelectionTools::CheatRecoPfoShowerSelector)
#endif
