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
        art::ValidHandle<std::vector<simb::MCParticle>> mcParticles,
        art::ValidHandle<std::vector<sim::SimChannel>> mcChannels,
        art::ValidHandle<std::vector<recob::Hit>> recoHits,
        art::ValidHandle<std::vector<recob::SpacePoint>> recoSpacePoints
    )
    {
        if (mcParticles.isValid() and mcChannels.isValid() and recoHits.isValis() and recoSpacePoints.isValid())
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
        }
    }
}