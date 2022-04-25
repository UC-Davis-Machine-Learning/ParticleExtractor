/**
 * @file RecoTraining.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "RecoTraining.h"

namespace extractor
{
    RecoTraining::RecoTraining()
    {
        fRecoTrainingTree = fTFileService->make<TTree>("reco_track_info", "reco_track_info");
        fRecoTrainingTree->Branch("sp_x", &fRecoTrainingSet.sp_x);
        fRecoTrainingTree->Branch("sp_y", &fRecoTrainingSet.sp_y);
        fRecoTrainingTree->Branch("sp_z", &fRecoTrainingSet.sp_z);
        fRecoTrainingTree->Branch("sp_pdg", &fRecoTrainingSet.sp_pdg);
        fRecoTrainingTree->Branch("ancestor_pdg", &fRecoTrainingSet.ancestor_pdg);
        fRecoTrainingTree->Branch("sp_track_id", &fRecoTrainingSet.sp_track_id);
        fRecoTrainingTree->Branch("ancestor_track_id", &fRecoTrainingSet.ancestor_track_id);
        fRecoTrainingTree->Branch("ancestor_level", &fRecoTrainingSet.ancestor_level);
        fRecoTrainingTree->Branch("summed_adc", &fRecoTrainingSet.summed_adc);
        fRecoTrainingTree->Branch("mean_adc", &fRecoTrainingSet.mean_adc);
        fRecoTrainingTree->Branch("peak_adc", &fRecoTrainingSet.peak_adc);
        fRecoTrainingTree->Branch("sigma_adc", &fRecoTrainingSet.sigma_adc);
        fRecoTrainingTree->Branch("dbscan_label", &fRecoTrainingSet.dbscan_label);
        fRecoTrainingTree->Branch("track_label", &fRecoTrainingSet.track_label);
    }

    RecoTraining::~RecoTraining()
    {}

    void RecoTraining::setBoundingBoxType(std::string volumeType)
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

    void RecoTraining::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::ValidHandle<std::vector<recob::Track>>& recoTracks,
        const art::FindManyP<recob::Hit>& hitPandoraSPsAssn, //to associate space points from pandora to hits
        const art::FindManyP<recob::Hit>& hitTrackAssn //to associate Tracks to hits
    )
    {
        RecoTrainingSet recoTrainingSet;
        if (mcParticles.isValid() and mcChannels.isValid() and recoSpacePoints.isValid() and recoTracks.isValid())
        {
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            std::map<Int_t, Int_t> parentDaughterMap;
            std::map<Int_t, Int_t> particlePDGMap;
            std::map<Int_t, Int_t> ancestorPDGMap;
            std::map<Int_t, Int_t> ancestorTrackIdMap;
            std::map<Int_t, Int_t> levelMap;
            for (auto particle : *mcParticles)
            {
                parentDaughterMap[particle.TrackId()] = particle.Mother();
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
                Int_t mother = particle.Mother();
                Int_t track_id = particle.TrackId();
                Int_t level = 0;
                while (mother != 0)
                {
                    level += 1;
                    track_id = mother;
                    mother = parentDaughterMap[track_id];
                }
                levelMap[particle.TrackId()] = level;
                ancestorPDGMap[particle.TrackId()] = particlePDGMap[track_id];
                ancestorTrackIdMap[particle.TrackId()] = track_id;
            }
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
                    trackHitList.emplace_back(hit);
                }  
            }
            //Making a map of grids and hits in them
            std::map<gridStruct, std::vector<hitStruct>> GridHitMap;
            makeGridHitMap(trackHitList, GridHitMap);
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
            
            std::vector<art::Ptr<recob::SpacePoint>> pointsList;
            art::fill_ptr_vector(pointsList, recoSpacePoints);            
            for (size_t i = 0; i < pointsList.size(); i++)
            {
                std::vector<Int_t> temp_pdg;
                std::vector<Int_t> temp_track_id;
                std::vector<Int_t> temp_ancestor_id;
                std::vector<Int_t> temp_ancestor_pdg;
                std::vector<Int_t> temp_label;
                std::vector<Double_t> temp_summed_adc;
                std::vector<Double_t> temp_mean_adc;
                std::vector<Double_t> temp_peak_adc;
                std::vector<Double_t> temp_sigma_adc;
                std::vector<Int_t> temp_level;
                auto& spsHit = hitSpacePointAssn.at(i);
                for (auto hit : spsHit)
                {  
                    Int_t track_id = TruthMatchUtils::TrueParticleID(
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
                    // find track hits
                    bool is_track_hit = searchGrid(hit, GridHitMap);
                    if (is_track_hit == 0)
                    {
                        if (ancestorPDGMap[track_id] == 2112)
                        {
                            temp_label.emplace_back(track_id);
                        } else {
                            temp_label.emplace_back(-1);
                        }
                        
                    } else {
                        temp_label.emplace_back(-2);
                    }
                    Int_t pdg = particlePDGMap[track_id];
                    Int_t ancestor = ancestorTrackIdMap[track_id];
                    Int_t ancestorpdg = ancestorPDGMap[track_id];
                    Int_t level = levelMap[track_id];
                    
                    temp_pdg.emplace_back(pdg);
                    temp_track_id.emplace_back(track_id);
                    temp_ancestor_id.emplace_back(ancestor);
                    temp_ancestor_pdg.emplace_back(ancestorpdg);
                    temp_summed_adc.emplace_back(hit->SummedADC());
                    temp_mean_adc.emplace_back(hit->PeakTime());
                    temp_peak_adc.emplace_back(hit->PeakAmplitude());
                    temp_sigma_adc.emplace_back(hit->RMS());
                    temp_level.emplace_back(level);
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
                recoTrainingSet.sp_x.emplace_back(xyz[0]);
                recoTrainingSet.sp_y.emplace_back(xyz[1]);
                recoTrainingSet.sp_z.emplace_back(xyz[2]);

                // now decide how to assign labels
                // take the value with the largest summed adc
                int max_index = std::distance(
                    temp_summed_adc.begin(), 
                    std::max_element(temp_summed_adc.begin(), temp_summed_adc.end())
                );

                recoTrainingSet.sp_pdg.emplace_back(temp_pdg[max_index]);
                recoTrainingSet.sp_track_id.emplace_back(temp_track_id[max_index]);
                recoTrainingSet.ancestor_track_id.emplace_back(temp_ancestor_id[max_index]);
                recoTrainingSet.ancestor_pdg.emplace_back(temp_ancestor_pdg[max_index]);
                recoTrainingSet.summed_adc.emplace_back(temp_summed_adc[max_index]);
                recoTrainingSet.mean_adc.emplace_back(temp_mean_adc[max_index]);
                recoTrainingSet.peak_adc.emplace_back(temp_peak_adc[max_index]);
                recoTrainingSet.sigma_adc.emplace_back(temp_sigma_adc[max_index]);
                recoTrainingSet.ancestor_level.emplace_back(temp_level[max_index]);
                recoTrainingSet.track_label.emplace_back(temp_label[max_index]);
            }
        }
        fRecoTrainingSet = recoTrainingSet;
        fRecoTrainingTree->Fill();
    }
}