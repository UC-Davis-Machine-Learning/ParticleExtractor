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
        std::vector<std::vector<Int_t>> adc;
        std::vector<std::vector<Int_t>> track_id;
        std::vector<std::vector<Int_t>> pdg;
        std::vector<std::vector<Int_t>> ancestor_id;
        std::vector<std::vector<Int_t>> level;
    };

    class RawDecoder
    {
    public:
        RawDecoder();
        ~RawDecoder();
        
        void setPDGCodes(std::vector<Int_t> PDGCodes) { fPDGCodes = PDGCodes; }
        void setPDGLevels(std::vector<Int_t> PDGLevels) { fPDGLevels = PDGLevels; }
        void setPDGLevels(std::vector<std::string> PDGLevels);

        void processEvent(
            const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
            const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
            const art::ValidHandle<std::vector<raw::RawDigit>>& rawDigits,
        );
    private:
        /// ROOT output through art::TFileService
        /** We will save different TTrees to different TFiles specified 
         *  by the directories for each type.
         */ 
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fRawDecoderTree;

        // pdg codes to construct
        std::vector<Int_t> fPDGCodes;
        std::vector<Int_t> fPDGLevels;
        Double_t fEnergyCutoff;

        // struct for holding event information
        WireTDC fWireTDC;
    };
}