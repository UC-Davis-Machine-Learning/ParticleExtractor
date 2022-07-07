/**
 * @file RawTrainingSet.cxx
 * @author Nicholas Carrara [nmcarrara@ucdavis.edu]
 * @author Yashwanth Bezawada [ysbezawada@ucdavis.edu]
 * @brief 
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "RawTrainingSet.h"

namespace extractor
{
    RawTrainingSet::RawTrainingSet()
    {
        fRawTrainingSetTree = fTFileService->make<TTree>("raw_training_set", "raw_training_set");
        fRawTrainingSetTree->Branch("RawTimeChanU", &fRawTimeChanU);
        fRawTrainingSetTree->Branch("RawTimeChanV", &fRawTimeChanV);
        fRawTrainingSetTree->Branch("RawTimeChanZ", &fRawTimeChanZ);
        fRawTrainingSetTree->Branch("TruthTimeChanU", &fTruthTimeChanU);
        fRawTrainingSetTree->Branch("TruthTimeChanV", &fTruthTimeChanV);
        fRawTrainingSetTree->Branch("TruthTimeChanZ", &fTruthTimeChanZ);
    }

    RawTrainingSet::~RawTrainingSet()
    {}

    void RawTrainingSet::setBoundingBoxType(std::string volumeType)
    {
        if (volumeType == "TPC" or volumeType == "tpc") { 
            fBoundingBoxType = VolumeType::TPC;
        }
        else if (volumeType == "Cryo" or volumeType == "cryo") {
            fBoundingBoxType = VolumeType::Cryostat;
        }
        else {
            fBoundingBoxType = VolumeType::World;
        }
    }

    //Histograms to store raw and truth data 
    void RawTrainingSet::MakeHist()
    {
        // place to define the histograms

        art::ServiceHandle<art::TFileService> tfs;
        //Histogram names and titles                                                                                                                                                         
        std::stringstream  name, title;

        //Channel boundries
        unsigned int UChMin;
        unsigned int UChMax;
        unsigned int VChMin;
        unsigned int VChMax;
        unsigned int ZChMin;
        unsigned int ZChMax;
        TH2I* TempHisto;


        // Accquiring geometry data
        fNofAPA=fGeom->NTPC()*fGeom->Ncryostats()/2; //No. of APAs
        fChansPerAPA = fGeom->Nchannels()/fNofAPA; //No. of channels per APA

        //To get max TDC
        auto const *fDetProp = lar::providerFrom<detinfo::DetectorPropertiesService>();
        fNticks = fDetProp->NumberTimeSamples();

        // taken from dune35t module a way to organise the channel mapping:
        // loop through channels in the first APA to find the channel boundaries for each view
        // will adjust for desired APA after
        fUChanMin = 0;
        fZChanMax = fChansPerAPA - 1;
        for ( unsigned int c = fUChanMin + 1; c < fZChanMax; c++ ){
            if ( fGeom->View(c) == geo::kV && fGeom->View(c-1) == geo::kU ){
                fVChanMin = c;
                fUChanMax = c - 1;
            }
            if ( fGeom->View(c) == geo::kZ && fGeom->View(c-1) == geo::kV ){
                fZChanMin = c;
                fVChanMax = c-1;
            }
        }

        //Number of channels in each view
        fNUCh=fUChanMax-fUChanMin+1; //U
        fNVCh=fVChanMax-fVChanMin+1; //V
        fNZCh=fZChanMax-fZChanMin+1; //Z (collection plane)

        
        unsigned int minT = 0;
        unsigned int maxT = 0;
        minT = 0;
        maxT = fNticks;
        unsigned int binT = (maxT-minT); //Bin width for TDC

        for(unsigned int i=0;i<fNofAPA;i++){
            UChMin=fUChanMin + i*fChansPerAPA;
            UChMax=fUChanMax + i*fChansPerAPA;
            VChMin=fVChanMin + i*fChansPerAPA;
            VChMax=fVChanMax + i*fChansPerAPA;
            ZChMin=fZChanMin + i*fChansPerAPA;
            ZChMax=fZChanMax + i*fChansPerAPA;

            // construct the histograms to store Raw data; TH2 constructors: ("Name", "Title", NxBin, xMin, xMax, NyBin, yMin, yMax)
            name.str("");
            name << "fRawTimeChanU";
            name <<  i;
            title.str("");
            title << "Raw Time vs Channel(Plane U, APA";
            title << i<<")";
            TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), UChMax - UChMin + 1, UChMin, UChMax, binT, minT, maxT);
            fRawTimeChanU.push_back(TempHisto);

            name.str("");
            name << "fRawTimeChanV";
            name << i;
            title.str("");
            title << "Raw Time vs Channel(Plane V, APA";
            title << i<<")";
            TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), VChMax - VChMin + 1, VChMin, VChMax, binT, minT, maxT);
            fRawTimeChanV.push_back(TempHisto);

            name.str("");
            name << "fRawTimeChanZ";
            name << i;
            title.str("");
            title << "Raw Time vs Channel(Plane Z, APA";
            title <<i<<")";
            TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), ZChMax - ZChMin + 1, ZChMin, ZChMax, binT, minT, maxT);
            fRawTimeChanZ.push_back(TempHisto);


            fRawTimeChanU[i]->SetStats(0);
            fRawTimeChanV[i]->SetStats(0);    
            fRawTimeChanZ[i]->SetStats(0);    


            fRawTimeChanU[i]->GetXaxis()->SetTitle("Channel"); fRawTimeChanU[i]->GetYaxis()->SetTitle("TDC");
            fRawTimeChanV[i]->GetXaxis()->SetTitle("Channel"); fRawTimeChanV[i]->GetYaxis()->SetTitle("TDC");
            fRawTimeChanZ[i]->GetXaxis()->SetTitle("Channel"); fRawTimeChanZ[i]->GetYaxis()->SetTitle("TDC");

            //////////////////////////////////////////////////////////////////

            // construct the histograms to store Truth data; TH2 constructors: ("Name", "Title", NxBin, xMin, xMax, NyBin, yMin, yMax)
            name.str("");
            name << "fTruthTimeChanU";
            name <<  i;
            title.str("");
            title << "Truth Time vs Channel(Plane U, APA";
            title << i<<")";
            TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), UChMax - UChMin + 1, UChMin, UChMax, binT, minT, maxT);
            fTruthTimeChanU.push_back(TempHisto);

            name.str("");
            name << "fTruthTimeChanV";
            name << i;
            title.str("");
            title << "Truth Time vs Channel(Plane V, APA";
            title << i<<")";
            TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), VChMax - VChMin + 1, VChMin, VChMax, binT, minT, maxT);
            fTruthTimeChanV.push_back(TempHisto);

            name.str("");
            name << "fTruthTimeChanZ";
            name << i;
            title.str("");
            title << "Truth Time vs Channel(Plane Z, APA";
            title <<i<<")";
            TempHisto = tfs->make<TH2I>(name.str().c_str(),title.str().c_str(), ZChMax - ZChMin + 1, ZChMin, ZChMax, binT, minT, maxT);
            fTruthTimeChanZ.push_back(TempHisto);


            fTruthTimeChanU[i]->SetStats(0);
            fTruthTimeChanV[i]->SetStats(0);    
            fTruthTimeChanZ[i]->SetStats(0);    


            fTruthTimeChanU[i]->GetXaxis()->SetTitle("Channel"); fTruthTimeChanU[i]->GetYaxis()->SetTitle("TDC");
            fTruthTimeChanV[i]->GetXaxis()->SetTitle("Channel"); fTruthTimeChanV[i]->GetYaxis()->SetTitle("TDC");
            fTruthTimeChanZ[i]->GetXaxis()->SetTitle("Channel"); fTruthTimeChanZ[i]->GetYaxis()->SetTitle("TDC");
        }
    }

    Int_t RawTrainingSet::getTrackID(std::vector< sim::TrackIDE > trackInfo, std::map<Int_t, Int_t>& parentDaughterMap){
        Int_t track_id;
        float_t track_energy;
        size_t track_index;
        track_energy = trackInfo[0].energy;
        track_index = 0;
        for (size_t i = 1; i < trackInfo.size(); i++)
        {
            if (trackInfo[i].energy > track_energy)
            {
                track_energy = trackInfo[i].energy;
                track_index = i;
            }
        }

        track_id = trackInfo[track_index].trackID;

        Int_t mother = parentDaughterMap[track_id];
        while (mother != 0)
        {
            track_id = mother;
            mother = parentDaughterMap[track_id];
        }
        return track_id;
    }

    void RawTrainingSet::processEvent(
        detinfo::DetectorClocksData const& clockData,
        const art::ValidHandle<std::vector<simb::MCParticle>>& mcParticles,
        const art::ValidHandle<std::vector<sim::SimChannel>>& mcChannels,
        const art::ValidHandle<std::vector<raw::RawDigit>>& rawTPC
    )
    {
        // RawTrainingSetStruct RawTrainingSet;
        if (mcParticles.isValid() and mcChannels.isValid() and rawTPC.isValid())
        {
            /**
             * We first iterate through all particles and create a map of 
             * parent-daughter pairs for track ids.  This way we can search
             * recursively for ancestors of each particle.
             */
            std::map<Int_t, Int_t> parentDaughterMap;
            std::map<Int_t, Int_t> particlePDGMap;
            std::map<Int_t, Int_t> secondaryPDGMap;
            std::map<Int_t, Int_t> ancestorPDGMap;
            std::map<Int_t, Int_t> secondaryTrackIdMap;
            std::map<Int_t, Int_t> ancestorTrackIdMap;
            std::map<Int_t, Int_t> levelMap;
            for (auto particle : *mcParticles)
            {
                parentDaughterMap[particle.TrackId()] = particle.Mother();
                particlePDGMap[particle.TrackId()] = particle.PdgCode();
            }
            for (auto particle : *mcParticles)
            {
                Int_t mother = particle.Mother();
                Int_t track_id = particle.TrackId();
                Int_t prev_track_id = 0;
                Int_t level = 0;
                while (mother != 0)
                {
                    level += 1;
                    prev_track_id = track_id;
                    track_id = mother;
                    mother = parentDaughterMap[track_id];
                }
                levelMap[particle.TrackId()] = level;
                secondaryPDGMap[particle.TrackId()] = particlePDGMap[prev_track_id];
                ancestorPDGMap[particle.TrackId()] = particlePDGMap[track_id];
                secondaryTrackIdMap[particle.TrackId()] = prev_track_id;
                ancestorTrackIdMap[particle.TrackId()] = track_id;
            }

            MakeHist();

            // Fill pointer vectors - more useful form for the raw data
            // a more usable form
            std::vector< art::Ptr<raw::RawDigit> > RawDigits;
            art::fill_ptr_vector(RawDigits, rawTPC);

            for(auto const & dptr : RawDigits) {
                const raw::RawDigit & digit = *dptr;

                // Get the channel number for this digit
                uint32_t chan = digit.Channel();
                // number of samples in uncompressed ADC
                int nSamples = digit.Samples();
                unsigned int apa = std::floor( chan/fChansPerAPA );	  
                int pedestal = (int)digit.GetPedestal();
                
                std::vector<Int_t> uncompressed(nSamples);
                // with pedestal	  
                raw::Uncompress(digit.ADCs(), uncompressed, pedestal, digit.Compression());
                // subtract pedestals
                std::vector<Int_t> uncompPed(nSamples);
                for (int i=0; i<nSamples; i++) uncompPed.at(i)=uncompressed.at(i)-pedestal;
                
                // number of ADC uncompressed without pedestal
                nADC_uncompPed=uncompPed.size();	  

                //Truth Channel
                auto truth_channel = mcChannels->at(chan);

                //Induction Plane   
                if( fGeom->View(chan) == geo::kU){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(uncompPed.at(l)!=0){
                            //Truth
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }
                            Int_t track_id = getTrackID(trackIDs, parentDaughterMap);
                            fTruthTimeChanU[apa]->Fill(truth_channel, l, particlePDGMap[track_id]);
                            //Raw data
                            fRawTimeChanU[apa]->Fill(chan, l, (Int_t) uncompPed.at(l));
                        }
                    }
                }// end of U View

                //Induction Plane   
                if( fGeom->View(chan) == geo::kV){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(uncompPed.at(l)!=0){
                            //Truth
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }
                            Int_t track_id = getTrackID(trackIDs, parentDaughterMap);
                            fTruthTimeChanV[apa]->Fill(truth_channel, l, particlePDGMap[track_id]);
                            //Raw data
                            fRawTimeChanV[apa]->Fill(chan,l, (Int_t) uncompPed.at(l));
                        }
                    }
                }// end of V View

                //Collection Plane
                if ( fGeom->View(chan) == geo::kZ){
                    for(unsigned int l=0;l<nADC_uncompPed;l++) {
                        if(uncompPed.at(l)!=0){
                            //Truth
                            auto const& trackIDs = truth_channel.TrackIDEs(l, l);
                            if (trackIDs.size() == 0) {
                                continue;
                            }
                            Int_t track_id = getTrackID(trackIDs, parentDaughterMap);
                            fTruthTimeChanZ[apa]->Fill(truth_channel, l, particlePDGMap[track_id]);
                            //Raw data
                            fRawTimeChanZ[apa]->Fill(chan,l, (Int_t) uncompPed.at(l));
                        }
                    }
                }// end of Z View

            } //End of loop over raw digits

        }
        // fRawTrainingSet = RawTrainingSet;
        fRawTrainingSetTree->Fill();
    }
}