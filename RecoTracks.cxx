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
        std::vector<recob::Hit>& List,
        std::map<gridStruct, std::vector<recob::Hit>>& Map
    )
    {
        for(size_t j=0;j< List.size();j++)
        {
            gridStruct grid;
            grid.gridPT = (int) (List[j].Channel()/50) + 1;
            grid.gridCID = (int) (List[j].PeakTime()/250) + 1;
            
            std::map<gridStruct, std::vector<recob::Hit>>::iterator gridItr = Map.find(grid);

            if(gridItr != Map.end()){
                Map[grid].insert( List[j] );
            } else {
                Map.insert( make_pair(grid, std::vector<recob::Hit>()) );
                Map[grid].insert( List[j] );
            }
        }
    }

    bool Recotracks::searchGrid(
        art::Ptr<recob::Hit> hit,
        std::map<gridStruct, std::vector<recob::Hit>>& Map
    )
    {
        gridStruct grid;
        grid.gridPT = (int) (hit->Channel()/50) + 1;
        grid.gridCID = (int) (hit->PeakTime()/250) + 1;

        std::map<gridStruct, std::vector<recob::Hit>>::iterator gridItr = Map.find(grid);

        if(gridItr != Map.end())
        {
            std::vector<recob::Hit>::iterator ptr;
            for (ptr = gridItr->second.begin(); ptr < gridItr->second.end(); ptr++)
            {
                if (*ptr.Channel() == hit->Channel() && *ptr.PeakTime() == hit->PeakTime())
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
            std::vector<recob::Hit> trackHitList;
            for (size_t i = 0; i < trackList.size(); i++)
            {
                auto& allHits = hitTrackAssn.at(i); //storing hits for ith track
                trackHitList.insert(std::end(trackHitList), std::begin(allHits), std::end(allHits));
            }

            //Making a map of grids and hits in them
            std::map<gridStruct, std::vector<recob::Hit>> GridHitMap;
            makeGridHitMap(trackHitList, GridHitMap);

            std::vector<art::Ptr<recob::SpacePoint>> pointsList;
            art::fill_ptr_vector(pointsList, recoSpacePoints);            
            for (size_t i = 0; i < pointsList.size(); i++)
            {
                std::vector<Int_t> temp_label;
                std::vector<Double_t> temp_summed_adc;
                auto& spsHit = hitPandoraSPsAssn.at(i);
                auto num_channels = mcChannels->size();
                for (auto hit : spsHit)
                {  
                    // find the corresponding sim channels
                    Int_t track_id;
                    // check if hit channel is reached
                    if (hit->Channel() >= num_channels) {
                        break;
                    }
                    auto channel = mcChannels->at(hit->Channel());
                    auto const& trackIDs = channel.TrackIDEs((int)hit->PeakTime(), (int)hit->PeakTime());
                    if (trackIDs.size() == 0) {
                        continue;
                    }
                    // track_id = trackIDs[0].trackID;
                    track_id = TruthMatchUtils::TrueParticleID(
                        clockData, hit, false
                    );
                    // check that track_id is present in parentDaughterMap
                    if (parentDaughterMap.find(track_id) == parentDaughterMap.end())
                    {
                        std::cout << "Track ID: " << track_id << " does not have an associated mother!" << std::endl;
                        continue;
                    }
                    if (particlePDGMap.find(track_id) == particlePDGMap.end())
                    {
                        std::cout << "Track ID: " << track_id << " does not have an associated pdg!" << std::endl;
                        continue;
                    }
                    Int_t mother = parentDaughterMap[track_id];
                    Int_t level = 0;
                    while (mother != 0)
                    {
                        level += 1;
                        track_id = mother;
                        mother = parentDaughterMap[track_id];
                    }

                    bool isTrackHit = searchGrid(hit, GridHitMap);
                    if (isTrackHit == 0)
                    {
                        if (particlePDGMap[track_id] == 2112)
                        {
                            temp_label.emplace_back(track_id);
                        } else {
                            temp_label.emplace_back(-1);
                        }
                        
                    } else {
                        temp_label.emplace_back(-2);
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