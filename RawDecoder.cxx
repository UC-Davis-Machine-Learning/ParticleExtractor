/**
 * @file RawDecoder.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Junying Huang [jyghuang@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-04-01
 */
#include "RawDecoder.h"

namespace extractor
{
    RawDecoder::RawDecoder()
    {
        fWireTDCTree = fTFileService->make<TTree>("raw_decoder", "raw_decoder");
        fWireTDCTree->Branch("pdg", &fWireTDC.pdg);
        fWireTDCTree->Branch("track_id", &fWireTDC.track_id);
        fWireTDCTree->Branch("ancestor_id", &fWireTDC.ancestor_id);
        fWireTDCTree->Branch("level", &fWireTDC.level);
        fWireTDCTree->Branch("adc", &fWireTDC.adc);
    }
    RawDecoder::~RawDecoder()
    {

    }
        
    void RawDecoder::setPDGLevels(std::vector<std::string> PDGLevels)
    {
        std::vector<Int_t> levels;
        for (size_t i = 0; i < PDGLevels.size(); i++)
        {
            if (PDGLevels[i] == "parent") {
                levels.emplace_back(0);
            }
            else if (PDGLevels[i] == "daughters") { 
                levels.emplace_back(1);
            }
            else if (PDGLevels[i] == "electrons") {
                levels.emplace_back(2);
            }
            else if (PDGLevels[i] == "parent_electrons") {
                levels.emplace_back(3);
            }
            else {
                levels.emplace_back(4);
            }
        }
        fPDGLevels = levels;
    }

    void RawDecoder::processEvent(
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<raw::RawDigit>>& rawDigits,
    )
    {
        WireTDC wireTDC;
        if (mcParticles.isValid() and mcChannels.isValid() and rawDigits.isValid())
        {
        }
        fWireTDC = wireTDC;
        fWireTDCTree->Fill();
    }
    
}