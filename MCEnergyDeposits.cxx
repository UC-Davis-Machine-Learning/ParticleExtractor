/**
 * @file MCEnergyDeposits.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-16
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "MCEnergyDeposits.h"

namespace extractor
{
    MCEnergyDeposits::MCEnergyDeposits()
    {}

    MCEnergyDeposits::~MCEnergyDeposits()
    {}

    void MCEnergyDeposits::processEvent(
        art::ValidHandle<std::vector<simb::MCParticle>> mcParticles,
        art::ValidHandle<std::vector<sim::SimEnergyDeposit>> mcEnergyDeposits
    )
    {
        
    }
}