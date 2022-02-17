/**
 * @file RecoEnergyDeposits.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
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

    }

    RecoEnergyDeposits::~RecoEnergyDeposits()
    {}

    void RecoEnergyDeposits::processEvent(
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<recob::Hit>>& recoHits,
        const art::ValidHandle<std::vector<recob::SpacePoint>>& recoSpacePoints,
        const art::FindManyP<recob::Hit>& hitSpacePointAssn
    )
    {
        RecoEdep recoEdep;
        if (mcParticles.isValid() and mcChannels.isValid() and recoHits.isValid() and recoSpacePoints.isValid())
        {
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            std::map<Int_t, Int_t> parentDaughterMap;
            std::map<Int_t, Int_t> particlePDGMap;
            for (auto particle : *mcParticles)
            {
                parentDaughterMap[particle.TrackId()] = particle.Mother();
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
            }
            // get a list of hits associated to each space point
            

            /**
             * We iterate over the list of hits ...
             * 
             */
            for (auto hit : *recoHits)
            {
                // If the hit is not in the collection plane,
                // then just continue.
                if (hit.WireID().Plane != 2) {
                    continue;
                }
                // now find the corresponding sim channels
                for (auto channel : *mcChannels)
                {
                    // the hit and sim channels must match
                    if (channel.Channel() != hit.Channel()) {
                        continue;
                    }
                    auto trackIDs = channel.TrackIDEs(hit.PeakTime(), hit.PeakTime());
                    if (trackIDs.size() != 0)
                    {

                    }

                }

                for (size_t i = 0; i < recoSpacePoints.size(); i++)
                {
                    for (size_t j = 0; j < hitSpacePointAssn[i].size(); j++)
                    {
                        if (hitSpacePointAssn[i][j]->WireID().Plane != 2) {
                            continue;
                        }
                        for (size_t k = 0; k < recoHits.size(); k++)
                        {
                            if (
                                (hitSpacePointAssn[i][j]->Channel() != recoHits[k].Channel()) or
                                (hitSpacePointAssn[i][j]->PeakTime() != recoHits[k].PeakTime())
                            ) {
                                continue;
                            }
                            auto xyz = spacePoint.XYZ();
                            recoEdep.sp_x.emplace_back(xyz[0]);
                            recoEdep.sp_y.emplace_back(xyz[1]);
                            recoEdep.sp_z.emplace_back(xyz[2]);
                        }
                    }
                }
            }
        }
    }
}