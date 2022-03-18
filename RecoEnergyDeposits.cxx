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
    }

    RecoEnergyDeposits::~RecoEnergyDeposits()
    {}

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
                std::cout << "point: " << i << std::endl;
                std::vector<Int_t> temp_pdg;
                std::vector<Int_t> temp_track_id;
                std::vector<Int_t> temp_ancestor_id;
                std::vector<Int_t> temp_channel_id;
                std::vector<Double_t> temp_summed_adc;
                Int_t hit_count = 0;
                auto& spsHit = hitSpacePointAssn.at(i);
                Int_t num_channels = mcChannels->size();
                for (auto hit : spsHit)
                {   
                    std::cout << "  hit: " << hit_count << std::endl;
                    hit_count += 1;
                    // If the hit is not in the collection plane,
                    // then just continue.
                    if (hit->WireID().Plane != 2) {
                        std::cout << "wireplane != 2" << std::endl;
                        continue;
                    }
                    // now find the corresponding sim channels
                    Int_t track_id;
                    Int_t channel_count = 0;
                    // check if hit channel is reached
                    if (hit->Channel() >= num_channels) {
                        break;
                    }
                    auto channel = mcChannels->at(hit->Channel());
                    auto const& trackIDs = channel.TrackIDEs((int)hit->PeakTime(), (int)hit->PeakTime());
                    if (trackIDs.size() != 0)
                    {
                        std::cout << "          track size > 0" << std::endl;
                        temp_channel_id.emplace_back(channel.Channel());
                        temp_track_id.emplace_back(trackIDs[0].trackID);
                        track_id = trackIDs[0].trackID;
                        break;
                    }
                    Int_t mother = parentDaughterMap[track_id];
                    std::cout << "      mother: " << mother << std::endl;
                    while (mother != 0)
                    {
                        track_id = mother;
                        mother = parentDaughterMap[track_id];
                    }
                    std::cout << "      mother: " << mother << std::endl;
                    temp_ancestor_id.emplace_back(track_id);
                    temp_pdg.emplace_back(particlePDGMap[track_id]);
                    temp_summed_adc.emplace_back(hit->SummedADC());
                }
                // collect results
                recoEdep.pdg.emplace_back(temp_pdg);
                recoEdep.track_id.emplace_back(temp_track_id);
                recoEdep.ancestor_id.emplace_back(temp_ancestor_id);
                recoEdep.channel_id.emplace_back(temp_channel_id);
                recoEdep.summed_adc.emplace_back(temp_summed_adc);

                auto xyz = pointsList[i]->XYZ();
                recoEdep.sp_x.emplace_back(xyz[0]);
                recoEdep.sp_y.emplace_back(xyz[1]);
                recoEdep.sp_z.emplace_back(xyz[2]);
            }
        }
        fRecoEdep = recoEdep;
        fRecoEnergyDepositsTree->Fill();
    }
}