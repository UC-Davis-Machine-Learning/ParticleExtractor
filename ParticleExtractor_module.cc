/**
 * @file    ParticleExtractor_module.cc
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief   A module for extracting truth/reco information about G4 particle trajectories
 *          and conducting some standard analysis tasks. 
 * @ingroup ParticleExtractor
**/

// art framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"


// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// necessary ROOT libraries
#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"

// C includes
#include <cmath>

// local includes
#include "Configuration.h"
#include "DetectorGeometry.h"
#include "MCNeutronCaptures.h"

namespace extractor
{
    class ParticleExtractor : public art::EDAnalyzer
    {
    public:
        explicit ParticleExtractor(const Parameters& config);
        ParticleExtractor(const ParticleExtractor&) = delete;
        ParticleExtractor(ParticleExtractor&&) = delete;
        ParticleExtractor& operator=(const ParticleExtractor&) = delete;
        ParticleExtractor& operator=(ParticleExtractor&&) = delete;

        // required EDAnalyzer functions
        void analyze(const art::Event& event) override;
        void beginJob() override;
        void endJob() override;

    private:
        /// Set of parameters
        Parameters fParameters;
        bool    fFillMCNeutronCaptures;
        std::vector<Int_t> fPDGCodes;

        // producer labels
        art::InputTag fLArGeantProducerLabel;

        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        /// TTrees
        TTree *fMetaTree;

        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("ParticleExtractor");
        // MC neutron captures
        MCNeutronCaptures fMCNeutronCaptures;
    };

    // constructor
    ParticleExtractor::ParticleExtractor(const Parameters& config)
    : EDAnalyzer(config)
    , fParameters(config)
    {
        /**
         * Here we check the various parameter settings and ...
         * 
         */
        fFillMCNeutronCaptures = fParameters().FillMCNeutronCaptures();

        // Producer labels
        fLArGeantProducerLabel = fParameters().LArGeantProducerLabel();

        // PDG codes
        fPDGCodes = fParameters().PDGCodes();

    }

    // begin job
    void ParticleExtractor::beginJob()
    {
        fGeometry->FillTTree();
    }

    // analyze
    void ParticleExtractor::analyze(const art::Event& event)
    {
        /**
         * @details For each event, we will look through the various
         * available data products and send event info to the 
         * corresponding submodules that process them, starting with MCParticles
         * 
         */
        art::Handle<std::vector<simb::MCParticle>> particleHandle;
        if (!event.getByLabel(fLArGeantProducerLabel, particleHandle))
        {
            // if there are no particles for the event truth, then
            // we are in big trouble haha.  throw an exception
            throw cet::exception("ParticleExtractor")
                << " No simb::MCParticle objects in this event - "
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        // get the list of MC particles from Geant4
        auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);
        // now pass the list of particles to each of the appropriate submodules
        if (fFillMCNeutronCaptures) {
            fMCNeutronCaptures.process(mcParticles);
        }
    }

    // end job
    void ParticleExtractor::endJob()
    {

    }

}

DEFINE_ART_MODULE(extractor::ParticleExtractor)