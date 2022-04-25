/**
 * @file RecoTraining.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada
 * @author Junying Haung
 * @brief 
 * @version 0.1
 * @date 2022-04-25
 */
#pragma once


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimPhotons.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "larsim/Utils/TruthMatchUtils.h"

#include <vector>
#include <algorithm>
#include <map>

#include "DetectorGeometry.h"
#include "RecoTracks.h"

namespace detinfo {
  class DetectorClocksData;
}

namespace extractor
{
    struct RecoTrainingSet
    {
        std::vector<double> sp_x;
        std::vector<double> sp_y;
        std::vector<double> sp_z;
        std::vector<int>    sp_pdg;
        std::vector<int>    ancestor_pdg;
        std::vector<int>    sp_track_id;
        std::vector<int>    ancestor_track_id;
        std::vector<int>    ancestor_level;
        std::vector<double> summed_adc;
        std::vector<double> mean_adc;
        std::vector<double> peak_adc;
        std::vector<double> sigma_adc;
        std::vector<int>    dbscan_label;
        std::vector<int>    track_label;
    };

    class RecoTraining
    {
    public:
        RecoTraining();
        ~RecoTraining();

        void setBoundingBoxType(std::string volumeType);

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
            const art::ValidHandle<std::vector<recob::Track>>& recoTracks,
            const art::FindManyP<recob::Hit>& hitPandoraSPsAssn, //to associate space points from pandora to hits
            const art::FindManyP<recob::Hit>& hitTrackAssn //to associate Tracks to hits
        );

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fRecoTrainingTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("RecoTraining");

        // pdg codes to construct
        VolumeType fBoundingBoxType;

        // struct for holding event information
        RecoTrainingSet fRecoTrainingSet;
    };
}