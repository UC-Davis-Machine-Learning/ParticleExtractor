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
        fRawDecoderTree = fTFileService->make<TTree>("raw_decoder", "raw_decoder");
        fRawDecoderTree->Branch("scTrackID", &fWireTDC.scTrackID);
        fRawDecoderTree->Branch("scChannelID", &fWireTDC.scChannelID);
        fRawDecoderTree->Branch("ancestor_id", &fWireTDC.scPeakTime);
        fRawDecoderTree->Branch("scAncestorPDG", &fWireTDC.scAncestorPDG);
        fRawDecoderTree->Branch("scAncestor", &fWireTDC.scAncestor);
        fRawDecoderTree->Branch("level", &fWireTDC.level);
    }
    RawDecoder::~RawDecoder()
    {

    }
        

    void RawDecoder::processEvent(
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& scs
    )
    {
        getmother.clear();
        getpdg.clear();
        WireTDC wireTDC;
        if (mcParticles.isValid() and scs.isValid())
        {
            for(auto &trueParticle : *mcParticles) {
                getmother.insert(std::pair<int,int>(trueParticle.TrackId(),trueParticle.Mother()));
                getpdg.insert(std::pair<int,int>(trueParticle.TrackId(),trueParticle.PdgCode()));
            }
            int mother=-1;
            int mothertemp=-1;
            int pdg=-1;
            int scid=-1;
            for(auto &sc : *scs) {
                for(int pt=0;pt<6000;pt++){
                    auto simChannelNumber = sc.Channel();
                    auto const& trackInfo=sc.TrackIDEs(pt, pt);
                    if(trackInfo.size()!=0){
                        fWireTDC.scChannelID.push_back(simChannelNumber);
                        fWireTDC.scTrackID.push_back(trackInfo[0].trackID);
                        fWireTDC.scPeakTime.push_back(pt);
                        scid=trackInfo[0].trackID;
                    }
                    /*std::cout<<"tracksize: "<<trackInfo.size()<<std::endl;
                    for(int j = 0; j < (int) trackInfo.size(); ++j){
                        std::cout<<"trackid: "<<trackInfo[j].trackID<<std::endl;
                    }*/
                    if(scid!=-1){
                        for(auto &trueParticle : *mcParticles) {
                            auto mcid=trueParticle.TrackId();
                            if (mcid != scid) continue;
                            mother = trueParticle.Mother();
                            mothertemp=scid;
                            while (mother != 0)
                            {
                                mothertemp=mother;
                                mother=getmother[mother];
                            }
                            pdg=getpdg[mothertemp];
                            fWireTDC.scAncestor.push_back(mothertemp);
                            fWireTDC.scAncestorPDG.push_back(pdg);
                            break;
                        }
                    }
                    
                }
            }
        fWireTDC = wireTDC;
        fRawDecoderTree->Fill();
    }
    
}
}
