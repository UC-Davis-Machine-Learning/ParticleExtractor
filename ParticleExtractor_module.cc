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
#include <algorithm>
#include <chrono>
#include <ctime> 

// local includes
#include "Configuration.h"
#include "DetectorGeometry.h"
#include "MCNeutronCaptures.h"
#include "MCEnergyDeposits.h"
#include "MCVoxels.h"
#include "RecoEnergyDeposits.h"

namespace extractor
{
    /**
     * @brief 
     * 
     */
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
        bool fFillMCNeutronCaptures;
        bool fFillMCEnergyDeposits;
        bool fFillMCVoxels;
        bool fFillRecoEnergyDeposits;
        // MC edep variables
        std::vector<Int_t> fMCEdepPDGCodes;
        std::vector<std::string> fMCEdepPDGLevels;
        // MC voxel variables
        std::vector<Int_t> fMCEdepPDGLabels;
        Double_t fMCVoxelSize;
        std::string fMCVoxelBoundingBox;
        std::string fMCVoxelLabeling;

        // producer labels
        art::InputTag fLArGeantProducerLabel;
        art::InputTag fIonAndScintProducerLabel;
        art::InputTag fSimChannelProducerLabel;
        art::InputTag fSimChannelInstanceProducerLabel;
        art::InputTag fHitProducerLabel;
        art::InputTag fSpacePointProducerLabel;

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
        // MC EnergyDeposits
        MCEnergyDeposits fMCEnergyDeposits;
        // MC Voxels
        MCVoxels fMCVoxels;
        // Reco EnergyDeposits
        RecoEnergyDeposits fRecoEnergyDeposits;
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
        // Producer labels
        fLArGeantProducerLabel = fParameters().LArGeantProducerLabel();
        fIonAndScintProducerLabel = fParameters().IonAndScintProducerLabel();
        fSimChannelProducerLabel = fParameters().SimChannelProducerLabel();
        fSimChannelInstanceProducerLabel = fParameters().SimChannelInstanceProducerLabel();
        fHitProducerLabel = fParameters().HitProducerLabel();
        fSpacePointProducerLabel = fParameters().SpacePointProducerLabel();

        // Which submodules to run
        fFillMCNeutronCaptures = fParameters().FillMCNeutronCaptures();
        fFillMCEnergyDeposits = fParameters().FillMCEnergyDeposits();
        fFillMCVoxels = fParameters().FillMCVoxels();
        fFillRecoEnergyDeposits = fParameters().FillRecoEnergyDeposits();

        // MC edep information
        fMCEdepBoundingBox = fParameters().MCEdepBoundingBox();
        fMCEdepPDGCodes = fParameters().MCEdepPDGCodes();
        fMCEdepPDGLevels = fParameters().MCEdepPDGLevels();
        fMCEdepPDGLabels = fParameters().MCEdepPDGLabels();

        // MC Voxel information
        fMCVoxelSize = fParameters().MCVoxelSize();
        fMCVoxelBoundingBox = fParameters().MCVoxelBoundingBox();
        fMCVoxelLabeling = fParameters().MCVoxelLabeling();

        // check for errors
        if (fMCEdepPDGCodes.size() != fMCEdepPDGLevels.size())
        {
            throw cet::exception("ParticleExtractor")
                << " Configuration parameters 'MCEdepPDGCodes' and 'MCEdepPDGLevels'"
                << " have different numbers of entries, (" << fMCEdepPDGCodes.size() << " and "
                << fMCEdepPDGLevels.size() << ") but must be the same!\n"
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        for (size_t i = 0; i < fMCEdepPDGLevels.size(); i++)
        {
            if (std::find(
                    allowed_mc_edep_levels.begin(), 
                    allowed_mc_edep_levels.end(), 
                    fMCEdepPDGLevels[i]) == allowed_mc_edep_levels.end())
            {
                throw cet::exception("ParticleExtractor")
                << " Parameter i: '" << fMCEdepPDGLevels[i] << "' is not an allowed type for MCEdepPDGLevels!" 
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
            }
        }
        if (fFillMCVoxels and (fMCEdepPDGCodes.size() != fMCEdepPDGLabels.size()))
        {
            throw cet::exception("ParticleExtractor")
                << " Configuration parameters 'MCEdepPDGCodes' and 'MCEdepPDGLabels'"
                << " have different numbers of entries, (" << fMCEdepPDGCodes.size() << " and "
                << fMCEdepPDGLabels.size() << ") but must be the same!\n"
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (fFillMCEnergyDeposits != fFillMCVoxels)
        {
            throw cet::exception("ParticleExtractor")
                << " Both 'FillMCEnergyDeposits' and 'FillMCVoxels' must either both be true or false!"
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (fMCVoxelSize <= 0.0)
        {
            throw cet::exception("ParticleExtractor")
                << " MC Voxel configuration parameter 'MCVoxelSize' cannot be <= 0, "
                << "but was set to: " << fMCVoxelSize << "!\n"
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (std::find(
                allowed_bounding_boxes.begin(), 
                allowed_bounding_boxes.end(), 
                fMCVoxelBoundingBox) == allowed_bounding_boxes.end())
        {
            throw cet::exception("ParticleExtractor")
                << " Parameter 'MCVoxelBoundingBox': '" << fMCVoxelBoundingBox << "' is not an allowed type for MCVoxelBoundingBox!" 
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (std::find(
                allowed_mc_voxel_labeling.begin(), 
                allowed_mc_voxel_labeling.end(), 
                fMCVoxelLabeling) == allowed_mc_voxel_labeling.end())
        {
            throw cet::exception("ParticleExtractor")
                << " Parameter 'MCVoxelLabeling': '" << fMCVoxelLabeling << "' is not an allowed type for MCVoxelLabeling!" 
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        fMCEnergyDeposits.setBoundingBoxType(fMCEdepBoundingBox);
        fMCEnergyDeposits.setPDGCodes(fMCEdepPDGCodes);
        fMCEnergyDeposits.setPDGLevels(fMCEdepPDGLevels);
        fMCVoxels.setPDGCodes(fMCEdepPDGCodes);
        fMCVoxels.setVoxelLabels(fMCEdepPDGLabels);
        fMCVoxels.setVoxelSize(fMCVoxelSize);
        fMCVoxels.setBoundingBox(fMCVoxelBoundingBox);
        fMCVoxels.setVoxelLabeling(fMCVoxelLabeling);

        fMetaTree = fTFileService->make<TTree>("meta", "meta");
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
        auto mcEnergyDeposit = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fIonAndScintProducerLabel);

        // now pass the list of particles to each of the appropriate submodules
        if (fFillMCNeutronCaptures) {
            fMCNeutronCaptures.processEvent(mcParticles);
        }
        if (fFillMCEnergyDeposits) {
            fMCEnergyDeposits.processEvent(mcParticles, mcEnergyDeposit);
        }
        if (fFillMCVoxels) {
            fMCVoxels.processEvent(fMCEnergyDeposits);
        }
        if (fFillRecoEnergyDeposits) 
        {
            auto mcSimChannels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(fSimChannelProducerLabel.label(), fSimChannelInstanceProducerLabel.label())
                );
            auto recoHits =  event.getValidHandle<std::vector<recob::Hit>>(fHitProducerLabel);
            auto recoSpacePoints =  event.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointProducerLabel);
            art::FindManyP<recob::Hit> hitsFromSpacePointsAssn(recoSpacePoints, event, fSpacePointProducerLabel);
            fRecoEnergyDeposits.processEvent(
                mcParticles, 
                mcSimChannels,
                recoSpacePoints,
                hitsFromSpacePointsAssn
            );
        }
    }

    // end job
    void ParticleExtractor::endJob()
    {
        // grab and save system info
        std::string user = str(std::getenv("USER"));
        std::string host = str(std::getenv("HOSTNAME"));
        std::string dir  = str(std::getenv("PWD"));
        // get current time
        auto end = std::chrono::system_clock::now();
        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        auto end_datetime = str(std::ctime(&end_time));

        fMetaTree->Branch("user", &user);
        fMetaTree->Branch("host", &host);
        fMetaTree->Branch("current_dir", &dir);
        fMetaTree->Branch("date", &end_datetime);

        // save configuration parameters
        fMetaTree->Branch("LArGeantProducerLabel", &fLArGeantProducerLabel);
        fMetaTree->Branch("IonAndScintProducerLabel", &fIonAndScintProducerLabel);
        fMetaTree->Branch("SimChannelProducerLabel", &fSimChannelProducerLabel);
        fMetaTree->Branch("SimChannelInstanceProducerLabel", &fSimChannelInstanceProducerLabel);
        fMetaTree->Branch("HitProducerLabel", &fHitProducerLabel);
        fMetaTree->Branch("SpacePointProducerLabel", &fSpacePointProducerLabel);
        fMetaTree->Branch("FillMCNeutronCaptures", &fFillMCNeutronCaptures);
        fMetaTree->Branch("FillMCEnergyDeposits", &fFillMCEnergyDeposits);
        fMetaTree->Branch("FillMCVoxels", &fFillMCVoxels);
        fMetaTree->Branch("FillRecoEnergyDeposits", &fFillRecoEnergyDeposits);
        fMetaTree->Branch("MCEdepBoundingBox", &fMCEdepBoundingBox);
        fMetaTree->Branch("MCEdepPDGCodes", &fMCEdepPDGCodes);
        fMetaTree->Branch("MCEdepPDGLevels", &fMCEdepPDGLevels);
        fMetaTree->Branch("MCEdepPDGLabels", &fMCEdepPDGLabels);
        fMetaTree->Branch("MCVoxelSize", &fMCVoxelSize);
        fMetaTree->Branch("MCVoxelBoundingBox", &fMCVoxelBoundingBox);
        fMetaTree->Branch("MCVoxelLabeling", &fMCVoxelLabeling);
 
        fMetaTree->Fill();
    }

}

DEFINE_ART_MODULE(extractor::ParticleExtractor)