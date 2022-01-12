/**
 * @file    ParticleExtractor_module.cc
 * @brief   A module for extracting truth/reco information about G4 particle trajectories
 *          and conducting some standard analysis tasks. 
 *          Generated at Mon Oct 11 11:21:12 2021 using cetskelgen
 * @ingroup ParticleExtractor
 * @author  Nicholas Carrara (nmcarrara@ucdavis.edu),
**/

// art framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Utilities/Exception.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"


// special utility includes
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "art_root_io/TFileService.h"


// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimEnergyDeposit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// necessary ROOT libraries
#include <TTree.h>
#include <TH1.h>
#include "TH1F.h"
#include "TGeoMaterial.h"
#include "TGeoElement.h"

// C includes
#include <cmath>

// local includes
#include "DetectorGeometry.h"

namespace extractor {
    class ParticleExtractor;
}

namespace extractor
{
    struct EventList
    {
        Int_t event_id;
        std::vector<Int_t> ids;     // list holding all unique id values
        std::vector<std::vector<Int_t>> particle_daughters;
        std::vector<Int_t> particle_pdgs;
        std::vector<Int_t> particle_ids;
        std::vector<Int_t> particle_parent_ids;
        std::vector<Double_t> particle_x;
        std::vector<Double_t> particle_y;
        std::vector<Double_t> particle_z;
        std::vector<Double_t> particle_edep_energy; 
        std::vector<Int_t> particle_edep_num_electrons;

        EventList(Int_t event) : event_id(event){}
    };

    struct ParentParticleList
    {
        std::vector<Int_t> tracks;
        std::vector<Int_t> mothers;

        ParticleParentList(){}
    };

    class ParticleExtractor : public art::EDAnalyzer
    {
    public:
        struct Config
        {
            fhicl::Sequence<int> PDGCodes
            {
                fhicl::Name("PDGCodes"),
                fhicl::Comment("particle type pdg codes to look for")
            };
            fhicl::Atom<art::InputTag> LArGeantProducerLabel
            {
                fhicl::Name("LArGeantProducerLabel"),
                fhicl::Comment("tag of the input data product with the largeant side of the simulation")
            };
            fhicl::Atom<art::InputTag> LArGeantEnergyDepositProducerLabel
            {
                fhicl::Name("LArGeantEnergyDepositProducerLabel"),
                fhicl::Comment("tag of the input data product with the largeant side of the simulation")
            };
            fhicl::Atom<art::InputTag> IonAndScintProducerLabel
            {
                fhicl::Name("IonAndScintProducerLabel"),
                fhicl::Comment("tag of the input data product with the ionization and scintillation simulation")
            };
            fhicl::Atom<art::InputTag> OutputFile
            {
                fhicl::Name("OutputFile"),
                fhicl::Comment("name of the file to output the neutron statistics to")
            };
            fhicl::Atom<bool> CollectDaughters
            {
                fhicl::Name("CollectDaughters"),
                fhicl::Comment("whether to keep all daughter particle information")
            };
            fhicl::Atom<bool> CollectAll
            {
                fhicl::Name("CollectAll"),
                fhicl::Comment("whether to collect all particles with listed pdg codes or only primaries")
            };
        };
    public:
        using Parameters = art::EDAnalyzer::Table<Config>;
        explicit ParticleExtractor(Parameters const& config);
        ParticleExtractor(ParticleExtractor const&) = delete;
        ParticleExtractor(ParticleExtractor&&) = delete;
        ParticleExtractor& operator=(ParticleExtractor const&) = delete;
        ParticleExtractor& operator=(ParticleExtractor&&) = delete;

        // required EDAnalyzer functions
        void analyze(art::Event const& event) override;
        void beginJob() override;
        void endJob() override;

        // special functions
        void FillTTree();
        bool checkEventIds(EventList eventList, Int_t trackId);
        int findEdepPosition(EventList eventList, Int_t trackId, Double_t x, Double_t y, Double_t z);

    private:
        std::vector<Int_t> fPdgCodes;
        art::InputTag fLArGeantProducerLabel;
        art::InputTag fLArGeantEnergyDepositProducerLabel;
        art::InputTag fIonAndScintProducerLabel;
        art::InputTag fOutputFileArt;
        bool fCollectDaughters;
        bool fCollectAll;
        // geometry information
        DetectorGeometry* fGeometry = DetectorGeometry::getInstance("ParticleExtractor");
        // ROOT
        art::ServiceHandle<art::TFileService> fTFileService;
        TTree *fMetaTree;
        TTree *fParticleTree;
        TTree *fParentParticleTree;
        // event variables
        int fRun;
        int fSubRun;
        int fEvent;

        // number of events
        Int_t fNumberOfEvents;
        EventList fTempEventList(0);
        ParticleParentList fTempParticleParentList;
    };

    // constructor
    ParticleExtractor::ParticleExtractor(Parameters const& config)
    : EDAnalyzer(config)
    , fPdgCodes(config().PDGCodes())
    , fLArGeantProducerLabel(config().LArGeantProducerLabel())
    , fLArGeantEnergyDepositProducerLabel(config().LArGeantEnergyDepositProducerLabel())
    , fIonAndScintProducerLabel(config().IonAndScintProducerLabel())
    , fOutputFileArt(config().OutputFile())
    , fCollectDaughters(config().CollectDaughters())
    , fCollectAll(config().CollectAll())
    {
        // here we initiate the TFile services for each of the trees
        // we're going to create.  
        fMetaTree = fTFileService->make<TTree>("meta", "meta");
        fParticleTree = fTFileService->make<TTree>("particle", "particle");
        fParticleParentTree = fFileService->make<TTree>("particle_parent", "particle_parent");
        consumes<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);

        fParticleTree->Branch("event_id", &fTempEventList.event_id);
        fParticleTree->Branch("ids", &fTempEventList.ids);
        fParticleTree->Branch("particle_daughters", &fTempEventList.particle_daughters);
        fParticleTree->Branch("particle_pdgs", &fTempEventList.particle_pdgs);
        fParticleTree->Branch("particle_ids", &fTempEventList.particle_ids);
        fParticleTree->Branch("particle_parent_ids", &fTempEventList.particle_parent_ids);
        fParticleTree->Branch("particle_x", &fTempEventList.particle_x);
        fParticleTree->Branch("particle_y", &fTempEventList.particle_y);
        fParticleTree->Branch("particle_z", &fTempEventList.particle_z);
        fParticleTree->Branch("particle_edep_energy", &fTempEventList.particle_edep_energy);
        fParticleTree->Branch("particle_edep_num_electrons", &fTempEventList.particle_edep_num_electrons);

        fParticleParentTree->Branch("track_ids", &fTempParticleParentList.tracks);
        fParticleParentTree->Branch("mothers", &fTempParticleParentList.mothers);
    }

    bool ParticleExtractor::checkEventIds(EventList eventList, Int_t trackId)
    {
        for (size_t k = 0; k < eventList.ids.size(); k++)
        {
            if (eventList.ids[k] == trackId) {
                return true;
            }
        }
        return false;
    }

    int ParticleExtractor::findEdepPosition(EventList eventList, Int_t trackId, Double_t x, Double_t y, Double_t z)
    {
        for (size_t k = 0; k < eventList.particle_ids.size(); k++)
        {
            if (eventList.particle_ids[k] == trackId) {
                if((eventList.particle_x[k] == x) &&
                   (eventList.particle_y[k] == y) &&
                   (eventList.particle_z[k] == z)
                )
                {
                    return k;
                }
            }
        }
        return -1;
    }

    // analyze function
    void ParticleExtractor::analyze(art::Event const& event)
    {
        if (event.isRealData())
        {
            // If we are looking at real data, then we need to stop the analysis
            // and back out.
            throw cet::exception("ParticleExtractor")
                << " Event contains real data - "
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        // define a "handle" to point to a vector of the objects
        art::Handle<std::vector<simb::MCParticle>> particleHandle;
        if (!event.getByLabel(fLArGeantProducerLabel, particleHandle))
        {
            // if there are no particles for the event truth, then
            // we are in big trouble haha.  throw an exception
            throw cet::exception("ParticleExtractor")
                << " No simb::MCParticle objects in this event - "
                << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
        }
        // get the event meta data 
        fRun    = event.run();
        fSubRun = event.subRun();
        fEvent  = event.id().event();

        fNumberOfEvents++;
        // create a new event list
        EventList eventList(fNumberOfEvents-1);
        ParticleParentList particleParentList();
        // get the list of MC particles from Geant4
        auto mcParticles = event.getValidHandle<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);
        // iterate over all MC particles and grab all particles specified 
        // in the pdg code list.  
        if (mcParticles.isValid())
        {
            for (auto particle : *mcParticles)
            {
                // add particle/parent pair to list
                particleParentList.tracks.emplace_back(particle.TrackId());
                particleParentList.mothers.emplace_back(particle.Mother());
                // check if the particle is the right pdg code
                for (size_t j = 0; j < fPdgCodes.size(); j++)
                {
                    if (particle.PdgCode() == fPdgCodes[j])
                    {
                        // check if particle should be recorded
                        if (particle.Mother() == 0 ||
                            (fCollectDaughters && checkEventIds(eventList,particle.TrackId())) ||
                            fCollectAll
                        )
                        {
                            eventList.ids.emplace_back(particle.TrackId());
                            std::vector<Int_t> daughters;
                            if (fCollectDaughters)
                            {
                                for (int k = 0; k < particle.NumberDaughters(); k++) {
                                    eventList.ids.emplace_back(particle.Daughter(k));
                                    daughters.emplace_back(particle.Daughter(k));
                                }
                            }
                            eventList.particle_daughters.emplace_back(daughters);
                            for (size_t k = 0; k < particle.NumberTrajectoryPoints(); k++)
                            {
                                DetectorVolume currentVolume = fGeometry->getVolume(
                                    particle.Vx(), particle.Vy(), particle.Vz()
                                );
                                if (currentVolume.material_name == "LAr")
                                {
                                    eventList.particle_pdgs.emplace_back(particle.PdgCode());
                                    eventList.particle_ids.emplace_back(particle.TrackId());
                                    eventList.particle_parent_ids.emplace_back(particle.Mother());
                                    eventList.particle_x.emplace_back(particle.Vx(k));
                                    eventList.particle_y.emplace_back(particle.Vy(k));
                                    eventList.particle_z.emplace_back(particle.Vz(k));
                                    eventList.particle_edep_energy.emplace_back(-1.);
                                    eventList.particle_edep_num_electrons.emplace_back(-1);
                                }
                            }
                        }
                    }
                }
            }
        }
        auto mcEnergyDeposit = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fIonAndScintProducerLabel);
        if (mcEnergyDeposit.isValid())
        {
            for (auto energyDeposit : *mcEnergyDeposit)
            {
                if (checkEventIds(eventList, energyDeposit.TrackID()))
                {
                    int index = findEdepPosition(
                        eventList, 
                        energyDeposit.TrackID(), 
                        energyDeposit.StartX(),
                        energyDeposit.StartY(),
                        energyDeposit.StartZ()
                    );
                    if (index != -1)
                    {
                        eventList.particle_edep_energy[index] = energyDeposit.Energy();
                        eventList.particle_edep_num_electrons[index] = energyDeposit.NumElectrons();
                    }
                }   
            }
        }
        fTempEventList.event_id = eventList[i].event_id;
        fTempEventList.ids = eventList[i].ids;
        fTempEventList.particle_daughters = eventList[i].particle_daughters;
        fTempEventList.particle_pdgs = eventList[i].particle_pdgs;
        fTempEventList.particle_ids = eventList[i].particle_ids;
        fTempEventList.particle_parent_ids = eventList[i].particle_parent_ids;
        fTempEventList.particle_x = eventList[i].particle_x;
        fTempEventList.particle_y = eventList[i].particle_y;
        fTempEventList.particle_z = eventList[i].particle_z;
        fTempEventList.particle_edep_energy = eventList[i].particle_edep_energy;
        fTempEventList.particle_edep_num_electrons = eventList[i].particle_edep_num_electrons;
        fParticleTree->Fill();

        fTempParticleParentList.tracks = particleParentList[i].tracks;
        fTempParticleParentList.mothers = particleParentList[i].mothers;
        fParticleParentTree->Fill();
    }
    // begin job
    void ParticleExtractor::beginJob()
    {
        fGeometry->FillTTree();
        fNumberOfEvents = 0;
    }
    // end job
    void ParticleExtractor::endJob()
    {
        // global neutron info
        fMetaTree->Branch("number_of_events", &fNumberOfEvents);
        fMetaTree->Fill();
    }
}

DEFINE_ART_MODULE(extractor::ParticleExtractor)