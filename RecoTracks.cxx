/**
 * @file RecoTracks.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "RecoTracks.h"

namespace extractor
{
    RecoTracks::RecoTracks()
    {
        fRecoTracksTree = fTFileService->make<TTree>("reco_track_info", "reco_track_info");
        fRecoTracksTree->Branch("label", &fRecoTrackHit.label);
        fRecoTracksTree->Branch("sp_x", &fRecoTrackHit.sp_x);
        fRecoTracksTree->Branch("sp_y", &fRecoTrackHit.sp_y);
        fRecoTracksTree->Branch("sp_z", &fRecoTrackHit.sp_z);
        fRecoTracksTree->Branch("summed_adc", &fRecoTrackHit.summed_adc);
    }

    RecoTracks::~RecoTracks()
    {}

    void RecoTracks::setBoundingBoxType(std::string volumeType)
    {
        if (volumeType == "TPC" or volumeType == "tpc") { 
            fBoundingBoxType = VolumeType::TPC;
        }
        else if (volumeType == "Cryo" or volumeType == "cryo") {
            fBoundingBoxType = VolumeType::Cryostat;
        }
        else {
            fBoundingBoxType = VolumeType::World;
        }
    }

    void RecoTracks::makeGridHitMap(
        std::vector<hitStruct>& List,
        std::map<int, std::vector<hitStruct>>& Map
    )
    {
        std::cout << "Making a Grid-Hit Map....." << std::endl;
        std::cout << "Is the hit list empty: " << List.empty() << std::endl;
        for(size_t j=0;j< List.size();j++)
        {
            // gridStruct grid;
            // grid.gridPT = (int) (List[j].cID/50) + 1;
            // grid.gridCID = (int) (List[j].PT/250) + 1;

            // Cantor pairing function
            // Pairing function to create a bijective NxN -> N mapping
            int x = (int) (List[j].cID/50) + 1;
            int y = (int) (List[j].PT/250) + 1;
            int gridNum = ( (x + y) * (x + y + 1) )/2 + y;

            std::map<int, std::vector<hitStruct>>::iterator gridItr = Map.find(gridNum);

            if(gridItr != Map.end()){
                Map[grid].push_back( List[j] );
            } else {
                Map.insert( make_pair(gridNum, std::vector<hitStruct>()) );
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

    bool RecoTracks::searchGrid(
        art::Ptr<recob::Hit> hit,
        std::map<gridStruct, std::vector<hitStruct>>& Map
    )
    {
        
        // Returns 0 if it's not a track spt
        // Returns 1 if it's a track spt
        
        // gridStruct grid;    
        // grid.gridPT = (int) (hit->Channel()/50) + 1;
        // grid.gridCID = (int) (hit->PeakTime()/250) + 1;

        // Cantor pairing function
        // Pairing function to create a bijective NxN -> N mapping
        int x = (int) (hit->Channel()/50) + 1;
        int y = (int) (hit->PeakTime()/250) + 1;
        int gridNum = ( (x + y) * (x + y + 1) )/2 + y;

        std::map<int, std::vector<hitStruct>>::iterator gridItr = Map.find(gridNum);

        if(gridItr != Map.end())
        {
            for(int i=0; i < (int) gridItr->second.size(); i++)
            {
                if (gridItr->second[i].cID == (int) hit->Channel() && gridItr->second[i].PT == (double) hit->PeakTime())
                {
                    return 1;
                }
            }
            return 0;   
        } else {
            return 0;
        }
    }

    void RecoTracks::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::ValidHandle<std::vector<recob::Track>>& recoTracks,
        const art::FindManyP<recob::Hit>& hitPandoraSPsAssn, //to associate space points from pandora to hits
        const art::FindManyP<recob::Hit>& hitTrackAssn //to associate Tracks to hits
    )
    {
        RecoTrackHit recoTrackHit;
        if (mcParticles.isValid() and mcChannels.isValid() and recoSpacePoints.isValid() and recoTracks.isValid())
        {
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            std::map<Int_t, Int_t> parentDaughterMap;
            std::map<Int_t, Int_t> particlePDGMap;
            std::map<Int_t, Int_t> ancestorMap;
            for (auto particle : *mcParticles)
            {
                parentDaughterMap[particle.TrackId()] = particle.Mother();
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
            }
            // get a list of hits associated to each space point
            /**
             * We iterate over the list of hits and look for the 
             * SimChannel which matches the hits Channel() value.
             * We then use the Hits PeakTime to find the track IDs
             * that live within that SimChannels PeakTime window, from 
             * which we can associate an MCParticle.
             *      MCParticle.TrackId() -> SimChannel.TrackIDEs(peaktime, peaktime) -> Hit.Channel()  
             * 
             */
            
            //Making a vector of track hits
            std::vector<art::Ptr<recob::Track>> trackList;
            art::fill_ptr_vector(trackList, recoTracks);
            std::vector<hitStruct> trackHitList;
            std::cout << "Making a list of track hits....." << std::endl;
            for (size_t i = 0; i < trackList.size(); i++)
            {
                std::vector<art::Ptr<recob::Hit>> allHits = hitTrackAssn.at(i); //storing hits for ith track
                for (size_t j = 0; j < allHits.size(); j++)
                {
                    hitStruct hit;
                    hit.cID = allHits[j]->Channel();
                    hit.PT = allHits[j]->PeakTime();
                    trackHitList.push_back(hit);
                }  
            }
            std::cout << "Is the track-hit list empty: " << trackHitList.empty() << std::endl;

            //Making a map of grids and hits in them
            std::map<int, std::vector<hitStruct>> GridHitMap;
            makeGridHitMap(trackHitList, GridHitMap);

            std::vector<art::Ptr<recob::SpacePoint>> pointsList;
            art::fill_ptr_vector(pointsList, recoSpacePoints);            
            for (size_t i = 0; i < pointsList.size(); i++)
            {
                std::vector<Int_t> temp_label;
                std::vector<Double_t> temp_summed_adc;
                auto& spsHit = hitPandoraSPsAssn.at(i);
                // auto num_channels = mcChannels->size();
                for (auto hit : spsHit)
                {  
                    
                    Int_t track_id;
                    track_id = TruthMatchUtils::TrueParticleID(
                        clockData, hit, false
                    );
                    
                    Int_t mother = parentDaughterMap[track_id];
                    //Int_t level = 0;
                    while (mother != 0)
                    {
                        //level += 1;
                        track_id = mother;
                        mother = parentDaughterMap[track_id];
                    }

                    // hitStruct HIT;
                    // HIT.cID = hit->Channel();
                    // HIT.PT = hit->PeakTime();

                    bool isTrackHit = searchGrid(hit, GridHitMap);
                    if ( isTrackHit == 1 ) //std::find(trackHitList.begin(), trackHitList.end(), HIT) != trackHitList.end()
                    {
                        temp_label.emplace_back(-2);                                               
                    } else {
                        if (particlePDGMap[track_id] == 2112)
                        {
                            temp_label.emplace_back(track_id);
                        } else {
                            temp_label.emplace_back(-1);
                        }
                    }
                    
                    temp_summed_adc.emplace_back(hit->SummedADC());
                }

                // collect results
                auto xyz = pointsList[i]->XYZ();
                // check if point is in active volume
                // Determine if edep is within the desired volume
                DetectorVolume edep_volume = fGeometry->getVolume(
                    xyz[0], xyz[1], xyz[2]
                );
                if (edep_volume.volume_type != fBoundingBoxType) {
                    continue;
                }  

                recoTrackHit.label.emplace_back(temp_label);
                recoTrackHit.sp_x.emplace_back(xyz[0]);
                recoTrackHit.sp_y.emplace_back(xyz[1]);
                recoTrackHit.sp_z.emplace_back(xyz[2]);
                recoTrackHit.summed_adc.emplace_back(temp_summed_adc);
            }
        }
        fRecoTrackHit = recoTrackHit;
        fRecoTracksTree->Fill();
    }
}