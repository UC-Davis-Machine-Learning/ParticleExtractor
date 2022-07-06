/**
 * @file RecoClustering.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada [ymbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
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

#include "DetectorGeometry.h"

namespace detinfo {
  class DetectorClocksData;
}

namespace extractor
{
    /**
     * @brief 
     * 
     */
    struct RecoEdep
    {
        std::vector<std::vector<Int_t>> pdg;
        std::vector<std::vector<Int_t>> track_id;
        std::vector<std::vector<Int_t>> ancestor_id;
        std::vector<std::vector<Int_t>> channel_id;
        std::vector<Double_t> sp_x;
        std::vector<Double_t> sp_y;
        std::vector<Double_t> sp_z;
        std::vector<std::vector<Double_t>> summed_adc;
        std::vector<Double_t> energy;
        std::vector<Int_t> num_electrons;
        std::vector<std::vector<Int_t>> level;
    };

    /**
     * @brief 
     * 
     */
    class RecoClustering
    {
    public:
        RecoClustering();
        ~RecoClustering();

        void setPDGCodes(std::vector<Int_t> PDGCodes) { fPDGCodes = PDGCodes; }
        void setPDGLevels(std::vector<Int_t> PDGLevels) { fPDGLevels = PDGLevels; }
        void setPDGLevels(std::vector<std::string> PDGLevels);

        void setBoundingBoxType(std::string volumeType);

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
            const art::FindManyP<recob::Hit>& hitSpacePointAssn
        );

        RecoEdep getRecoEdep() const { return fRecoEdep; }

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fRecoClusteringTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("RecoClustering");

        // pdg codes to construct
        VolumeType fBoundingBoxType;
        std::vector<Int_t> fPDGCodes;
        std::vector<Int_t> fPDGLevels;

        // struct for holding event information
        RecoEdep fRecoEdep;
    };
}