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
#include "RecoVoxels.h"
#include "RawDecoder.h"
#include "RecoTracks.h"
#include "RecoTraining.h"
#include "RecoNeutrons.h"
#include "RecoDBScan3D.h"

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
        bool fFillRawDecoder;
        bool fFillRecoEnergyDeposits;
        bool fFillRecoVoxels;
        bool fFillRecoTracks;
        bool fFillRecoTraining;
        bool fFillRecoNeutrons;
        bool fFillRecoDBScan3D;

        // MC edep variables
        std::string fMCEdepBoundingBox;
        std::vector<Int_t> fMCEdepPDGCodes;
        std::vector<std::string> fMCEdepPDGLevels;
        Double_t fMCEdepEnergyCutoff;

        // MC voxel variables
        std::vector<Int_t> fMCEdepPDGLabels;
        Double_t fMCVoxelSize;
        std::string fMCVoxelBoundingBox;
        std::string fMCVoxelLabeling;
        

        // Reco edep variables
        std::vector<Int_t> fRecoEdepPDGCodes;
        std::string fRecoEdepBoundingBox;
        Double_t fRecoEdepEnergyCutoff;

        // Reco voxel variables
        std::vector<Int_t> fRecoEdepPDGLabels;
        Double_t fRecoVoxelSize;
        std::string fRecoVoxelBoundingBox;
        std::string fRecoVoxelLabeling;

        // producer labels
        art::InputTag fLArGeantProducerLabel;
        art::InputTag fIonAndScintProducerLabel;
        art::InputTag fSimChannelProducerLabel;
        art::InputTag fSimChannelInstanceProducerLabel;
        art::InputTag fHitProducerLabel;
        art::InputTag fSpacePointProducerLabel;
        art::InputTag fPandoraLabel;
        art::InputTag fPandoraTrackLabel;
        art::InputTag fDBScan3DLabel;

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
        // Reco Voxels
        RecoVoxels fRecoVoxels;
        // rawdecoder
        RawDecoder fRawDecoder;
        // recotracks
        RecoTracks fRecoTracks;
        // reco trianing
        RecoTraining fRecoTraining;
        // reco neutrons
        RecoNeutrons fRecoNeutrons;
        // reco DBScan
        RecoDBScan3D fRecoDBScan3D;
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
        fPandoraLabel = fParameters().PandoraLabel();
        fPandoraTrackLabel = fParameters().PandoraTrackLabel();
        fDBScan3DLabel = fParameters().DBScan3DLabel();

        // Which submodules to run
        fFillMCNeutronCaptures = fParameters().FillMCNeutronCaptures();
        fFillMCEnergyDeposits = fParameters().FillMCEnergyDeposits();
        fFillMCVoxels = fParameters().FillMCVoxels();
        fFillRawDecoder = fParameters().FillRawDecoder();
        fFillRecoEnergyDeposits = fParameters().FillRecoEnergyDeposits();
        fFillRecoVoxels = fParameters().FillRecoVoxels();

        //RecoTracks
        fFillRecoTracks = fParameters().FillRecoTracks();

        //RecoTraining
        fFillRecoTraining = fParameters().FillRecoTraining();

        //RecoTraining
        fFillRecoNeutrons = fParameters().FillRecoNeutrons();

        //DBScan3D
        fFillRecoDBScan3D = fParameters().FillRecoDBScan3D();

        // MC edep information
        fMCEdepBoundingBox = fParameters().MCEdepBoundingBox();
        fMCEdepPDGCodes = fParameters().MCEdepPDGCodes();
        fMCEdepPDGLevels = fParameters().MCEdepPDGLevels();
        fMCEdepPDGLabels = fParameters().MCEdepPDGLabels();
        fMCEdepEnergyCutoff = fParameters().MCEdepEnergyCutoff();

        // MC Voxel information
        fMCVoxelSize = fParameters().MCVoxelSize();
        fMCVoxelBoundingBox = fParameters().MCVoxelBoundingBox();
        fMCVoxelLabeling = fParameters().MCVoxelLabeling();

        // Reco edep information
        fRecoEdepBoundingBox = fParameters().RecoEdepBoundingBox();

        // Reco Voxel information
        fRecoVoxelSize = fParameters().RecoVoxelSize();
        fRecoVoxelBoundingBox = fParameters().RecoVoxelBoundingBox();
        fRecoEdepPDGLabels = fParameters().RecoEdepPDGLabels();
        fRecoVoxelLabeling = fParameters().RecoVoxelLabeling();

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
        if (fFillMCVoxels and !fFillMCEnergyDeposits)
        {
            throw cet::exception("ParticleExtractor")
                << " If 'FillMCVoxels' set to true, then 'FillMCEnergyDeposits' must also be true!"
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
        if (fFillRecoVoxels and (fMCEdepPDGCodes.size() != fRecoEdepPDGLabels.size()))
        {
            throw cet::exception("ParticleExtractor")
                << " Configuration parameters 'MCEdepPDGCodes' and 'RecoEdepPDGLabels'"
                << " have different numbers of entries, (" << fMCEdepPDGCodes.size() << " and "
                << fRecoEdepPDGLabels.size() << ") but must be the same!\n"
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (fFillRecoVoxels and !fFillRecoEnergyDeposits)
        {
            throw cet::exception("ParticleExtractor")
                << " If 'FillRecoVoxels' set to true, then 'FillRecoEnergyDeposits' must also be true!"
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (fRecoVoxelSize <= 0.0)
        {
            throw cet::exception("ParticleExtractor")
                << " Reco Voxel configuration parameter 'RecoVoxelSize' cannot be <= 0, "
                << "but was set to: " << fRecoVoxelSize << "!\n"
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        if (std::find(
                allowed_mc_voxel_labeling.begin(), 
                allowed_mc_voxel_labeling.end(), 
                fRecoVoxelLabeling) == allowed_mc_voxel_labeling.end())
        {
            throw cet::exception("ParticleExtractor")
                << " Parameter 'RecoVoxelLabeling': '" << fRecoVoxelLabeling << "' is not an allowed type for MCVoxelLabeling!" 
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }

        fMCEnergyDeposits.setBoundingBoxType(fMCEdepBoundingBox);
        fMCEnergyDeposits.setPDGCodes(fMCEdepPDGCodes);
        fMCEnergyDeposits.setPDGLevels(fMCEdepPDGLevels);
        fMCEnergyDeposits.setEnergyCutoff(fMCEdepEnergyCutoff);

        fMCVoxels.setPDGCodes(fMCEdepPDGCodes);
        fMCVoxels.setVoxelLabels(fMCEdepPDGLabels);
        fMCVoxels.setVoxelSize(fMCVoxelSize);
        fMCVoxels.setBoundingBox(fMCVoxelBoundingBox);
        fMCVoxels.setVoxelLabeling(fMCVoxelLabeling);

        fRecoEnergyDeposits.setBoundingBoxType(fRecoEdepBoundingBox);

        fRecoVoxels.setPDGCodes(fMCEdepPDGCodes);
        fRecoVoxels.setVoxelLabels(fRecoEdepPDGLabels);
        fRecoVoxels.setVoxelSize(fRecoVoxelSize);
        fRecoVoxels.setBoundingBox(fRecoVoxelBoundingBox);
        fRecoVoxels.setVoxelLabeling(fRecoVoxelLabeling);

        fRecoTracks.setBoundingBoxType(fRecoEdepBoundingBox);

        fRecoTraining.setBoundingBoxType(fRecoEdepBoundingBox);

        fRecoNeutrons.setBoundingBoxType(fRecoEdepBoundingBox);

        fRecoDBScan3D.setBoundingBoxType(fRecoEdepBoundingBox);

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
        std::cout << "Collecting MC Particles.." << std::endl;
        auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);

        // now pass the list of particles to each of the appropriate submodules
        if (fFillRawDecoder) {
            std::cout << "Filling Raw Decoder.." << std::endl;
            auto mcSimChannels =
                event.getValidHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(fSimChannelProducerLabel.label(), fSimChannelInstanceProducerLabel.label())
                );
            fRawDecoder.processEvent(mcParticles, mcSimChannels);
        }
        if (fFillMCNeutronCaptures) {
            std::cout << "Filling MC Neutron Captures.." << std::endl;
            auto mcEnergyDeposit = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fIonAndScintProducerLabel);
            fMCNeutronCaptures.processEvent(mcParticles, mcEnergyDeposit);
        }
        if (fFillMCEnergyDeposits) {
            std::cout << "Filling MC Energy Deposits.." << std::endl;
            auto mcEnergyDeposit = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fIonAndScintProducerLabel);
            fMCEnergyDeposits.processEvent(mcParticles, mcEnergyDeposit);
        }
        if (fFillMCVoxels) {
            std::cout << "Filling MC Voxels.." << std::endl;
            fMCVoxels.processEvent(fMCEnergyDeposits);
        }
        if (fFillRecoEnergyDeposits) 
        {
            std::cout << "Filling Reco Energy Deposits.." << std::endl;
            auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)); 
            auto mcSimChannels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(fSimChannelProducerLabel.label(), fSimChannelInstanceProducerLabel.label())
                );
            auto recoHits =  event.getValidHandle<std::vector<recob::Hit>>(fHitProducerLabel);
            auto recoSpacePoints =  event.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointProducerLabel);
            art::FindManyP<recob::Hit> hitsFromSpacePointsAssn(recoSpacePoints, event, fSpacePointProducerLabel);
            fRecoEnergyDeposits.processEvent(
                clockData,
                mcParticles, 
                mcSimChannels,
                recoSpacePoints,
                hitsFromSpacePointsAssn
            );
        }
        if (fFillRecoVoxels) {
            std::cout << "Filling Reco Voxels.." << std::endl;
            fRecoVoxels.processEvent(fRecoEnergyDeposits);
        }

        if (fFillRecoTracks) {
            std::cout << "Filling Reco Tracks.." << std::endl;
            auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)); 
            auto mcSimChannels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(fSimChannelProducerLabel.label(), fSimChannelInstanceProducerLabel.label())
                );
            auto recoSpacePoints = event.getValidHandle<std::vector<recob::SpacePoint>>(fPandoraLabel);
            auto recoTracks = event.getValidHandle< std::vector<recob::Track> >(fPandoraTrackLabel);
            art::FindManyP<recob::Hit> hitsFromSpsPandoraAssn(recoSpacePoints, event, fPandoraLabel); //to associate space point from pandora to hit
            art::FindManyP<recob::Hit> hitsFromTracksAssn(recoTracks, event, fPandoraTrackLabel); // to associate tracks and hits
            fRecoTracks.processEvent(
                clockData,
                mcParticles, 
                mcSimChannels,
                recoSpacePoints,
                recoTracks,
                hitsFromSpsPandoraAssn,
                hitsFromTracksAssn
            );
        }
        if (fFillRecoTraining) {
            std::cout << "Filling Reco Training.." << std::endl;
            auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)); 
            auto mcSimChannels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(fSimChannelProducerLabel.label(), fSimChannelInstanceProducerLabel.label())
                );
            auto recoSpacePoints = event.getValidHandle<std::vector<recob::SpacePoint>>(fPandoraLabel);
            auto recoTracks = event.getValidHandle< std::vector<recob::Track> >(fPandoraTrackLabel);
            art::FindManyP<recob::Hit> hitsFromSpsPandoraAssn(recoSpacePoints, event, fPandoraLabel); //to associate space point from pandora to hit
            art::FindManyP<recob::Hit> hitsFromTracksAssn(recoTracks, event, fPandoraTrackLabel); // to associate tracks and hits
            fRecoTraining.processEvent(
                clockData,
                mcParticles, 
                mcSimChannels,
                recoSpacePoints,
                recoTracks,
                hitsFromSpsPandoraAssn,
                hitsFromTracksAssn
            );
        }
        if (fFillRecoDBScan3D) {
            std::cout << "Filling Reco DBScan3D.." << std::endl;
            auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)); 
            auto mcSimChannels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(fSimChannelProducerLabel.label(), fSimChannelInstanceProducerLabel.label())
                );
            auto recoSpacePoints = event.getValidHandle<std::vector<recob::SpacePoint>>(fSpacePointProducerLabel);
            auto recoSlices = event.getValidHandle< std::vector<recob::Slice> >(fDBScan3DLabel);
            art::FindManyP<recob::Hit> hitsFromSpacePointsAssn(recoSpacePoints, event, fSpacePointProducerLabel);
            art::FindManyP<recob::SpacePoint> spacePointSliceAssn(recoSlices, event, fDBScan3DLabel);
            fRecoDBScan3D.processEvent(
                clockData,
                mcParticles, 
                mcSimChannels,
                recoSpacePoints,
                recoSlices,
                hitsFromSpacePointsAssn,
                spacePointSliceAssn
            );
        }
        if (fFillRecoNeutrons) {
            std::cout << "Filling Reco Neutrons.." << std::endl;
            auto const clockData(art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event)); 
            auto mcSimChannels = 
                event.getValidHandle<std::vector<sim::SimChannel>>(
                    art::InputTag(fSimChannelProducerLabel.label(), fSimChannelInstanceProducerLabel.label())
                );
            auto recoSpacePoints = event.getValidHandle<std::vector<recob::SpacePoint>>(fPandoraLabel);
            art::FindManyP<recob::Hit> hitsFromSpsPandoraAssn(recoSpacePoints, event, fPandoraLabel); //to associate space point from pandora to hit
            fRecoNeutrons.processEvent(
                clockData,
                mcParticles, 
                mcSimChannels,
                recoSpacePoints,
                hitsFromSpsPandoraAssn
            );
        }
    }

    // end job
    void ParticleExtractor::endJob()
    {
        // grab and save system info
        // std::string user = std::to_string(std::getenv("USER"));
        // std::string host = std::to_string(std::getenv("HOSTNAME"));
        // std::string dir  = std::to_string(std::getenv("PWD"));
        // get current time
        // auto end = std::chrono::system_clock::now();
        // std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        // auto end_datetime = std::to_string(std::ctime(&end_time));

        // fMetaTree->Branch("user", &user);
        // fMetaTree->Branch("host", &host);
        // fMetaTree->Branch("current_dir", &dir);
        //fMetaTree->Branch("date", &end_datetime);

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
        fMetaTree->Branch("FillRecoVoxels", &fFillRecoVoxels);

        fMetaTree->Branch("MCEdepBoundingBox", &fMCEdepBoundingBox);
        fMetaTree->Branch("MCEdepPDGCodes", &fMCEdepPDGCodes);
        fMetaTree->Branch("MCEdepPDGLevels", &fMCEdepPDGLevels);
        fMetaTree->Branch("MCEdepPDGLabels", &fMCEdepPDGLabels);

        fMetaTree->Branch("MCVoxelSize", &fMCVoxelSize);
        fMetaTree->Branch("MCVoxelBoundingBox", &fMCVoxelBoundingBox);
        fMetaTree->Branch("MCVoxelLabeling", &fMCVoxelLabeling);

        fMetaTree->Branch("RecoEdepPDGLabels", &fRecoEdepPDGLabels);
        fMetaTree->Branch("RecoVoxelSize", &fRecoVoxelSize);
        fMetaTree->Branch("RecoVoxelLabeling", &fRecoVoxelLabeling);
 
        fMetaTree->Fill();
    }

}

DEFINE_ART_MODULE(extractor::ParticleExtractor)
