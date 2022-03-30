/**
 * @file RecoEnergyDeposits.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Junying Huang [jyghuang@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "RecoEnergyDeposits.h"

namespace extractor
{
    RecoEnergyDeposits::RecoEnergyDeposits()
    {
        fRecoEnergyDepositsTree = fTFileService->make<TTree>("reco_energy_deposits", "reco_energy_deposits");
        fRecoEnergyDepositsTree->Branch("pdg", &fRecoEdep.pdg);
        fRecoEnergyDepositsTree->Branch("track_id", &fRecoEdep.track_id);
        fRecoEnergyDepositsTree->Branch("ancestor_id", &fRecoEdep.ancestor_id);
        fRecoEnergyDepositsTree->Branch("channel_id", &fRecoEdep.channel_id);
        fRecoEnergyDepositsTree->Branch("sp_x", &fRecoEdep.sp_x);
        fRecoEnergyDepositsTree->Branch("sp_y", &fRecoEdep.sp_y);
        fRecoEnergyDepositsTree->Branch("sp_z", &fRecoEdep.sp_z);
        fRecoEnergyDepositsTree->Branch("summed_adc", &fRecoEdep.summed_adc);
        fRecoEnergyDepositsTree->Branch("level", &fRecoEdep.level);
    }

    RecoEnergyDeposits::~RecoEnergyDeposits()
    {}

    void RecoEnergyDeposits::setBoundingBoxType(std::string volumeType)
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

    void RecoEnergyDeposits::processEvent(
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::FindManyP<recob::Hit>& hitSpacePointAssn
    )
    {
        RecoEdep recoEdep;
        if (mcParticles.isValid() and mcChannels.isValid() and recoSpacePoints.isValid())
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
            std::vector<art::Ptr<recob::SpacePoint>> pointsList;
            art::fill_ptr_vector(pointsList, recoSpacePoints);            
            for (size_t i = 0; i < pointsList.size(); i++)
            {
                std::vector<Int_t> temp_pdg;
                std::vector<Int_t> temp_track_id;
                std::vector<Int_t> temp_ancestor_id;
                std::vector<Int_t> temp_channel_id;
                std::vector<Double_t> temp_summed_adc;
                std::vector<Int_t> temp_level;
                auto& spsHit = hitSpacePointAssn.at(i);
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
                    track_id = trackIDs[0].trackID;
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
                    temp_channel_id.emplace_back(channel.Channel());
                    temp_track_id.emplace_back(trackIDs[0].trackID);
                    Int_t mother = parentDaughterMap[track_id];
                    Int_t level = 0;
                    while (mother != 0)
                    {
                        level += 1;
                        track_id = mother;
                        mother = parentDaughterMap[track_id];
                    }
                    
                    temp_ancestor_id.emplace_back(track_id);
                    temp_pdg.emplace_back(particlePDGMap[track_id]);
                    temp_summed_adc.emplace_back(hit->SummedADC());
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
                recoEdep.sp_x.emplace_back(xyz[0]);
                recoEdep.sp_y.emplace_back(xyz[1]);
                recoEdep.sp_z.emplace_back(xyz[2]);

                recoEdep.pdg.emplace_back(temp_pdg);
                recoEdep.track_id.emplace_back(temp_track_id);
                recoEdep.ancestor_id.emplace_back(temp_ancestor_id);
                recoEdep.channel_id.emplace_back(temp_channel_id);
                recoEdep.summed_adc.emplace_back(temp_summed_adc);
                recoEdep.level.emplace_back(temp_level);
            }
        }
        fRecoEdep = recoEdep;
        fRecoEnergyDepositsTree->Fill();
    }
}