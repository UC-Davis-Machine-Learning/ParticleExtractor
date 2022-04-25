/**
 * @file RecoEnergyDeposits.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
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
#include <map>

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
    struct RecoTrackHit
    {
        std::vector<std::vector<Int_t>> label; //-1 - Track hit; > 0 - neutron
        std::vector<Double_t> sp_x;
        std::vector<Double_t> sp_y;
        std::vector<Double_t> sp_z;
        std::vector<std::vector<Double_t>> summed_adc;
    };

    struct hitStruct {
        double PT;
        int cID;

        bool operator==(const hitStruct &o) const {
            return PT == o.PT && cID == o.cID;
        }

        bool operator<(const hitStruct &o)  const {
            return PT < o.PT || (PT == o.PT && cID < o.cID);
        }
    };

    struct gridStruct {
        double gridPT;
        int gridCID;

        bool operator==(const gridStruct &o) const {
            return gridPT == o.gridPT && gridCID == o.gridCID;
        }

        bool operator<(const gridStruct &o)  const {
            return gridPT < o.gridPT || (gridPT == o.gridPT && gridCID < o.gridCID);
        }
    };

    

    /**
     * @brief 
     * 
     */
    class RecoTracks
    {
    public:
        RecoTracks();
        ~RecoTracks();

        void setPDGCodes(std::vector<Int_t> PDGCodes) { fPDGCodes = PDGCodes; }
        void setPDGLevels(std::vector<Int_t> PDGLevels) { fPDGLevels = PDGLevels; }
        void setPDGLevels(std::vector<std::string> PDGLevels);

        void setBoundingBoxType(std::string volumeType);

        void makeGridHitMap(
            std::vector<hitStruct>& List,
            std::map<gridStruct, std::vector<hitStruct>>& Map
        )
        {
            std::cout << "Making a Grid-Hit Map....." << std::endl;
            std::cout << "Is the hit list empty: " << List.empty() << std::endl;
            for(size_t j=0;j< List.size();j++)
            {
                gridStruct grid;
                grid.gridPT = (int) (List[j].cID/50) + 1;
                grid.gridCID = (int) (List[j].PT/250) + 1;
                
                std::map<gridStruct, std::vector<hitStruct>>::iterator gridItr = Map.find(grid);

                if(gridItr != Map.end()){
                    Map[grid].push_back( List[j] );
                } else {
                    Map.insert( make_pair(grid, std::vector<hitStruct>()) );
                    Map[grid].push_back( List[j] );
                }
            }
            if(Map.empty() == 1)
            {
                std::cout << "Failed to make a Grid-Hit Map" << std::endl;
            } else {
                std::cout << "Completed making a Grid-Hit Map" << std::endl;
            }
        }

        bool searchGrid(
            art::Ptr<recob::Hit> hit,
            std::map<gridStruct, std::vector<hitStruct>>& Map
        )
        {
            /*
            ** Returns 0 if it's not a track spt
            ** Returns 1 if it's a track spt
            */
            gridStruct grid;    
            grid.gridPT = (int) (hit->Channel()/50) + 1;
            grid.gridCID = (int) (hit->PeakTime()/250) + 1;

            std::map<gridStruct, std::vector<hitStruct>>::iterator gridItr = Map.find(grid);

            if(gridItr != Map.end())
            {
                for(int i=0; i < (int) gridItr->second.size(); i++)
                {
                    if (gridItr->second[i].cID == (int) hit->Channel() && gridItr->second[i].PT == (int) hit->PeakTime())
                    {
                        return 1;
                    }
                }

                return 0;   
            } else {
                return 0;
            }
        }

        void processEvent(
            detinfo::DetectorClocksData const& clockData,
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
            const art::ValidHandle<std::vector<recob::Track>>& recoTracks,
            const art::FindManyP<recob::Hit>& hitPandoraSPsAssn, //to associate space points from pandora to hits
            const art::FindManyP<recob::Hit>& hitTrackAssn //to associate Tracks to hits
        );

        RecoTrackHit getRecoTrackHit() const { return fRecoTrackHit; }

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fRecoTracksTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("RecoTracks");

        // pdg codes to construct
        VolumeType fBoundingBoxType;
        std::vector<Int_t> fPDGCodes;
        std::vector<Int_t> fPDGLevels;

        // struct for holding event information
        RecoTrackHit fRecoTrackHit;
    };
}