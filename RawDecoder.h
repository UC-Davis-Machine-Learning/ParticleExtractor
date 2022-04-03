/**
 * @file RawDecoder.h
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Junying Huang [jyghuang@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-01
 */
#pragma once
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include <vector>
#include <algorithm>
#include <numeric>
#include <map>

namespace extractor
{
    struct WireTDC
    {
        std::vector<std::vector<Int_t>> scTrackID;
        std::vector<std::vector<Int_t>> scChannelID;
        std::vector<std::vector<Int_t>> scPeakTime;
        std::vector<std::vector<Int_t>> scAncestorPDG;
        std::vector<std::vector<Int_t>> scAncestor;
        std::vector<std::vector<Int_t>> level;
    };

    class RawDecoder
    {
    public:
        RawDecoder();
        ~RawDecoder();

        void processEvent(
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& scs,
        );
    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fRawDecoderTree;
        std::map<int,int> getmother;
        std::map<int,int> getpdg;

        // struct for holding event information
        WireTDC fWireTDC;
    };
}
