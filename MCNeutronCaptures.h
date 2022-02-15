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

#include "DetectorGeometry.h"

namespace extractor
{
    /**
     * @brief 
     * 
     */
    struct MCNeutronStatistics
    {
        std::vector<Double_t> neutron_capture_x;    ///< capture location in x
        std::vector<Double_t> neutron_capture_y;    ///< capture location in y
        std::vector<Double_t> neutron_capture_z;    ///< capture location in z
        std::vector<Int_t> neutron_ids;             ///< track id of each neutron
        std::vector<Int_t> total_number_steps;      ///< total number of steps in the trajectory
        std::vector<Int_t> tpc_number_steps;        ///< number of steps inside TPC
        std::vector<Int_t> lar_number_steps;        ///< number of steps in LAr
    };

    class MCNeutronCaptures
    {
    public:
        MCNeutronCaptures();
        ~MCNeutronCaptures();

        void processEvent(ValidHandle<std::vector<simb::MCParticle>> mcParticles);

    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fMCNeutronCapturesTree;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("MCNeutronCaptures");

        MCNeutronStatistics fMCNeutronStatistics;
        std::vector<Double_t> energy;
        std::vector<Double_t> dEdx;
    };
}