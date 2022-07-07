/**
 * @file RawTrainingSet.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
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

// LArSoft includes
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"

#include "larsim/Utils/TruthMatchUtils.h"

#include <vector>
#include <algorithm>
#include <map>

#include "DetectorGeometry.h"
#include "TH2.h"

namespace detinfo {
  class DetectorClocksData;
}

namespace extractor
{
    // struct RawTrainingSetStruct
    // {
    //     std::vector<TH2I*> raw_TimeChanU;
    //     std::vector<TH2I*> raw_TimeChanV;
    //     std::vector<TH2I*> raw_TimeChanZ;

    //     std::vector<TH2I*> truth_TimeChanU;
    //     std::vector<TH2I*> truth_TimeChanV;
    //     std::vector<TH2I*> truth_TimeChanZ;
        
    // };

    class RawTrainingSet
    {
    public:
        RawTrainingSet();
        ~RawTrainingSet();

        void setBoundingBoxType(std::string volumeType);

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<raw::RawDigit>>& rawTPC
        );

        void MakeHist();
        Int_t getTrackID(std::vector< sim::TrackIDE > trackInfo, std::map<Int_t, Int_t>& parentDaughterMap);

        

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fRawTrainingSetTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("RawTrainingSet");

        geo::GeometryCore const * fGeom = &*(art::ServiceHandle<geo::Geometry>());

        // pdg codes to construct
        VolumeType fBoundingBoxType;

        // struct for holding event information
        // RawTrainingSetStruct fRawTrainingSetStruct;

        // TPC // Number of channels in each planes
        unsigned int fNUCh;
        unsigned int fNVCh;
        unsigned int fNZCh;

        // find channel boundaries for each view
        unsigned int fUChanMin;
        unsigned int fUChanMax;
        unsigned int fVChanMin;
        unsigned int fVChanMax;
        unsigned int fZChanMin;
        unsigned int fZChanMax;
        unsigned int fNticks;

        unsigned int fNofAPA; //Number of APAs
        unsigned int fChansPerAPA; //Number of channels in each APA


        //unsigned int fMinT, fMaxT, fMaxTimeRange; // unused

        //Stores channel number, TDC and ADC values
        //Each element is for an APA
        std::vector<TH2I*> fRawTimeChanU;
        std::vector<TH2I*> fRawTimeChanV;
        std::vector<TH2I*> fRawTimeChanZ;

        //Stores channel number, TDC and PDG values
        //Each element is for an APA
        std::vector<TH2I*> fTruthTimeChanU;
        std::vector<TH2I*> fTruthTimeChanV;
        std::vector<TH2I*> fTruthTimeChanZ;

        // define nADC counts for uncompressed vs compressed
        unsigned int nADC_uncompPed;
    };
}