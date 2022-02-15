/**
 * @file MCNeutronCaptures.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#pragma once

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

namespace extractor
{
    class MCNeutronCaptures
    {
    public:
        MCNeutronCaptures();
        ~MCNeutronCaptures();

        void process(ValidHandle<std::vector<simb::MCParticle>> mcParticles);
    private:
    };
}