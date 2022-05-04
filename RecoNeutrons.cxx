/**
 * @file RecoNeutrons.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "RecoNeutrons.h"

namespace extractor
{
    RecoNeutrons::RecoNeutrons()
    {
        fRecoNeutronsTree = fTFileService->make<TTree>("reco_training", "reco_training");
        fRecoNeutronsTree->Branch("sp_x", &fRecoNeutronsSet.sp_x);
        fRecoNeutronsTree->Branch("sp_y", &fRecoNeutronsSet.sp_y);
        fRecoNeutronsTree->Branch("sp_z", &fRecoNeutronsSet.sp_z);
        fRecoNeutronsTree->Branch("neutron_id", &fRecoNeutronsSet.neutron_id);
        fRecoNeutronsTree->Branch("gamma_id", &fRecoNeutronsSet.gamma_id);
        fRecoNeutronsTree->Branch("summed_adc", &fRecoNeutronsSet.summed_adc);
        fRecoNeutronsTree->Branch("mean_adc", &fRecoNeutronsSet.mean_adc);
        fRecoNeutronsTree->Branch("peak_adc", &fRecoNeutronsSet.peak_adc);
        fRecoNeutronsTree->Branch("sigma_adc", &fRecoNeutronsSet.sigma_adc);
    }

    RecoNeutrons::~RecoNeutrons()
    {}

    void RecoNeutrons::setBoundingBoxType(std::string volumeType)
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

    void RecoNeutrons::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::ValidHandle<std::vector<recob::Track>>& recoTracks,
        const art::FindManyP<recob::Hit>& hitPandoraSPsAssn, //to associate space points from pandora to hits
        const art::FindManyP<recob::Hit>& hitTrackAssn //to associate Tracks to hits
    )
    {
        RecoNeutronsSet RecoNeutronsSet;
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
            std::map<Int_t, bool> ancestorCapture;
            std::map<Int_t, Int_t> levelMap;

            std::vector<int> neutron_captures;
            std::vector<std::vector<int>> gamma_ids;

            for (auto particle : *mcParticles)
            {
                parentDaughterMap[particle.TrackId()] = particle.Mother();
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
            }
            for (auto particle : *mcParticles)
            {
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
                ancestorCapture[particle.TrackId()] = false;
            }
            for (auto particle : *mcParticles)
            {
                // check if the particle is a neutron
                if (particle.PdgCode() == 2112)
                {
                    if (particle.EndProcess() == "nCapture")
                    {
                        neutron_captures.emplace_back(particle.TrackId());
                        ancestorCapture[ancestorTrackIdMap[particle.TrackId()]] = true;
                        gamma_ids.emplace_back(std::vector<int>());
                    }
                }
                // check if the particle is a gamma
                else if (particle.PdgCode() == 22)
                {
                    for(size_t i = 0; i < neutron_captures.size(); i++)
                    {
                        if (neutron_captures[i] == particle.Mother())
                        {
                            gamma_ids[i].emplace_back(particle.TrackId());
                        }
                    }
                }
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
            std::cout << "Iterating over hits..." << std::endl;         
            for (size_t i = 0; i < pointsList.size(); i++)
            {
                auto& spsHit = hitPandoraSPsAssn.at(i);
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
                    if (ancestorCapture[track_id])
                    {
                        Int_t temp_neutron = 0;
                        Int_t temp_gamma = 0;
                        Int_t temp_neutron_index = 0;
                        Int_t temp_track_id = track_id;
                        Int_t mother = parentDaughterMap[temp_track_id];
                        bool neutron_found = false;
                        bool gamma_found = false;
                        while(!neutron_found)
                        {
                            for (size_t i = 0; i < neutron_captures.size(); i++) 
                            {
                                if (neutron_captures[i] == mother) {
                                    temp_neutron = mother;
                                    temp_neutron_index = i;
                                    neutron_found = true;
                                    break;
                                }
                            }
                            temp_track_id = mother;
                            mother = parentDaughterMap[temp_track_id];
                        }
                        if (neutron_found)
                        {
                            Int_t temp_track_id = track_id;
                            Int_t mother = parentDaughterMap[temp_track_id];
                            while(!gamma_found)
                            {
                                for (size_t i = 0; i < gamma_ids[temp_neutron_index].size(); i++) 
                                {
                                    if (gamma_ids[temp_neutron_index][i] == mother) {
                                        temp_gamma = mother;
                                        gamma_found = true;
                                        break;
                                    }
                                }
                                temp_track_id = mother;
                                mother = parentDaughterMap[temp_track_id];
                            }
                        }
                        if (gamma_found)
                        {
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
                            
                            RecoNeutronsSet.sp_x.emplace_back(xyz[0]);
                            RecoNeutronsSet.sp_y.emplace_back(xyz[1]);
                            RecoNeutronsSet.sp_z.emplace_back(xyz[2]);
                            RecoNeutronsSet.neutron_id.emplace_back(temp_neutron);
                            RecoNeutronsSet.gamma_id.emplace_back(temp_gamma);
                            RecoNeutronsSet.summed_adc.emplace_back(hit->SummedADC());
                            RecoNeutronsSet.mean_adc.emplace_back(hit->PeakTime());
                            RecoNeutronsSet.peak_adc.emplace_back(hit->PeakAmplitude());
                            RecoNeutronsSet.sigma_adc.emplace_back(hit->RMS());
                        }
                    }
                }
            }
        }
        fRecoNeutronsSet = RecoNeutronsSet;
        fRecoNeutronsTree->Fill();
    }
}