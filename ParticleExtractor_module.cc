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
        // event level information
        Int_t event_id;
        // all the unique particle ids for this event
        std::vector<Int_t> ids;
        // all the pdg codes for each unique particle id
        std::vector<Int_t> particle_pdgs;
        // all of the parent ids for each unique particle id
        std::vector<Int_t> particle_parent_ids;
        // the daughter particle ids for each unique particle id
        std::vector<std::vector<Int_t>> particle_daughters;
        // the list of ancestors for each particle
        std::vector<std::vector<Int_t>> particle_ancestors;

        // individual trajectory information
        std::vector<Int_t> particle_ids;
        std::vector<Double_t> particle_x;
        std::vector<Double_t> particle_y;
        std::vector<Double_t> particle_z;
        std::vector<Double_t> particle_edep_energy; 
        std::vector<Int_t> particle_edep_num_electrons;

        EventList(Int_t event) : event_id(event){}
    };

    struct ParticleParentList
    {
        Int_t event_id;
        std::vector<Int_t> tracks;
        std::vector<Int_t> mothers;
        ParticleParentList(Int_t event) : event_id(event){}
    };

    struct ParticleTrajectory
    {
        Int_t track_id;
        Int_t mother;
        Int_t pdg;
        std::vector<Double_t> t;
        std::vector<Double_t> x;
        std::vector<Double_t> y;
        std::vector<Double_t> z;
        std::vector<Double_t> edep_energy;
        std::vector<Double_t> edep_num_electrons;
        ParticleTrajectory(Int_t track_id, Int_t mother, Int_t pdg)
        : track_id(track_id), mother(mother), pdg(pdg)
        {}
    };

    struct ParticleTree
    {
        Int_t event_id;
        std::vector<ParticleTrajectory> particles;
        std::vector<Int_t> edge_start;
        std::vector<Int_t> edge_end;
        ParticleTree(Int_t event) : event_id(event){}
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
        Int_t findTrackId(EventList eventList, Int_t trackId);
        Int_t findEdepPosition(EventList eventList, Int_t trackId, Double_t x, Double_t y, Double_t z);
        Int_t getParent(ParticleParentList particleParentList, Int_t trackId);
        Int_t findParentLocation(ParticleTree particleTree, Int_t track_id, Double_t t, Double_t x, Double_t y, Double_t z);        
        Int_t findParentTree(std::vector<ParticleTree> particleTree, Int_t track_id);
        Int_t findParentTrajectory(std::vector<ParticleTree> particleTree, Int_t track_id);
        std::vector<Double_t>::iterator findClosestTime(ParticleTrajectory particleTrajectory, Double_t time);
    
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
        TTree *fEventTree;
        TTree *fParticleTree;
        TTree *fParticleParentTree;
        
        // event variables
        int fRun;
        int fSubRun;
        int fEvent;

        // number of events
        Int_t fNumberOfEvents;
        Int_t fNumberOfPrimaries;

        EventList fTempEventList;
        ParticleParentList fTempParticleParentList;
        ParticleTree fTempParticleTree;
        std::vector<ParticleTree> fParticleTreeList;
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
    , fTempEventList(0)
    , fTempParticleParentList(0)
    , fTempParticleTree(0)
    {
        // here we initiate the TFile services for each of the trees
        // we're going to create.  
        fMetaTree = fTFileService->make<TTree>("meta", "meta");
        fEventTree = fTFileService->make<TTree>("event", "event");
        fParticleTree = fTFileService->make<TTree>("particle", "particle");
        fParticleParentTree = fTFileService->make<TTree>("particle_parent", "particle_parent");
        consumes<std::vector<simb::MCParticle>>(fLArGeantProducerLabel);

        fEventTree->Branch("event_id", &fTempEventList.event_id);
        fEventTree->Branch("ids", &fTempEventList.ids);
        fEventTree->Branch("particle_daughters", &fTempEventList.particle_daughters);
        fEventTree->Branch("particle_pdgs", &fTempEventList.particle_pdgs);
        fEventTree->Branch("particle_ids", &fTempEventList.particle_ids);
        fEventTree->Branch("particle_parent_ids", &fTempEventList.particle_parent_ids);
        fEventTree->Branch("particle_x", &fTempEventList.particle_x);
        fEventTree->Branch("particle_y", &fTempEventList.particle_y);
        fEventTree->Branch("particle_z", &fTempEventList.particle_z);
        fEventTree->Branch("particle_edep_energy", &fTempEventList.particle_edep_energy);
        fEventTree->Branch("particle_edep_num_electrons", &fTempEventList.particle_edep_num_electrons);
        fEventTree->Branch("particle_ancestors", &fTempEventList.particle_ancestors);

        // fParticleTree->Branch("event_id", &fTempParticleTree.event_id);
        // fParticleTree->Branch("track_id", &fTempParticleTree.track_id);
        // fParticleTree->Branch("mother", &fTempParticleTree.mother);
        // fParticleTree->Branch("pdg", &fTempParticleTree.pdg);
        // fParticleTree->Branch("t", &fTempParticleTree.t);
        // fParticleTree->Branch("x", &fTempParticleTree.x);
        // fParticleTree->Branch("y", &fTempParticleTree.y);
        // fParticleTree->Branch("z", &fTempParticleTree.z);
        // fParticleTree->Branch("process", &fTempParticleTree.process);
        // fParticleTree->Branch("edep_energy", &fTempParticleTree.edep_energy);
        // fParticleTree->Branch("edep_num_electrons", &fTempParticleTree.edep_num_electrons);
        // fParticleTree->Branch("edge_start", &fTempParticleTree.edge_start);
        // fParticleTree->Branch("edge_end", &fTempParticleTree.edge_end);

        fParticleParentTree->Branch("track_ids", &fTempParticleParentList.tracks);
        fParticleParentTree->Branch("mothers", &fTempParticleParentList.mothers);
    }

    bool ParticleExtractor::checkEventIds(EventList eventList, Int_t trackId)
    {
        for (size_t k = 0; k < eventList.ids.size(); k++) {
            if (eventList.ids[k] == trackId)  {
                return true;
            }
        }
        return false;
    }

    Int_t ParticleExtractor::findTrackId(EventList eventList, Int_t trackId)
    {
        for (size_t k = 0; k < eventList.ids.size(); k++) {
            if (eventList.ids[k] == trackId)  {
                return k;
            }
        }
        return -1;
    }

    Int_t ParticleExtractor::findEdepPosition(EventList eventList, Int_t trackId, Double_t x, Double_t y, Double_t z)
    {
        for (size_t k = 0; k < eventList.particle_ids.size(); k++) {
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

    Int_t ParticleExtractor::getParent(ParticleParentList particleParentList, Int_t trackId)
    {
        for (size_t k = 0; k < particleParentList.tracks.size(); k++) {
            if (particleParentList.tracks[k] == trackId) {
                return particleParentList.mothers[k];
            }
        }
        return -1;
    }

    Int_t ParticleExtractor::findParentLocation(ParticleTree particleTree, 
        Int_t track_id, Double_t t, Double_t x, Double_t y, Double_t z)
    {
        for (size_t k = 0; k < particleTree.particles.size(); k++) {
            if (particleTree.particles[k].track_id == track_id) {
                for (size_t j = 0; j < particleTree.particles[k].x.size(); j++)
                {
                    if((particleTree.particles[k].x[j] == x) &&
                        (particleTree.particles[k].y[j] == y) &&
                        (particleTree.particles[k].z[j] == z)
                    )
                    {
                        return k;
                    }
                    if (particleTree.particles[k].t[j] == t)
                    {
                        return k;
                    }
                }
            }
        }
        return -1;
    }

    Int_t ParticleExtractor::findParentTrajectory(std::vector<ParticleTree> particleTree,
        Int_t track_id)
    {
        for (size_t k = 0; k < particleTree.size(); k++) {
            for (size_t j = 0; j < particleTree[k].particles.size(); j++) {
                if (particleTree[k].particles[j].track_id == track_id) {
                    return j;
                }
            }
        }
        return -1;
    }

    Int_t ParticleExtractor::findParentTree(std::vector<ParticleTree> particleTree,
        Int_t track_id)
    {
        for (size_t k = 0; k < particleTree.size(); k++) {
            for (size_t j = 0; j < particleTree[k].particles.size(); j++) {
                if (particleTree[k].particles[j].track_id == track_id) {
                    return k;
                }
            }
        }
        return -1;
    }

    std::vector<Double_t>::iterator ParticleExtractor::findClosestTime(ParticleTrajectory particleTrajectory, Double_t time)
    {
        for (size_t k = 1; k < particleTrajectory.t.size(); k++)
        {
            if (time < particleTrajectory.t[k])
            {
                return k-1;
            }
        }
        return particleTrajectory.t.end();
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
        fNumberOfPrimaries=-1;
        // create a new event list
        EventList eventList(fNumberOfEvents-1);
        ParticleParentList particleParentList(fNumberOfEvents-1);
        fParticleTreeList.clear();
        // main loop
        /*  this part of the code loops over several data products starting
         *  with simb::MCParticle, from which it extracts trajectory information
         *  for all particles with pdg codes listed in the PDGCodes config variable.
         *  
         */
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
                // generate vertices and edges
                fParticleTreeList.emplace_back(ParticleTree(fEvent));
                fNumberOfPrimaries++;
                // fill primary tree with vertices
                fParticleTreeList[fNumberOfPrimaries].particles.emplace_back(
                    ParticleTrajectory(particle.TrackId(), particle.Mother(), particle.PdgCode())
                );
                if (particle.Mother() == 0)
                {
                    for (size_t k = 0; k < particle.NumberTrajectoryPoints(); k++)
                    {
                        fParticleTreeList[fNumberOfPrimaries].particles[0].t.emplace_back(particle.T(k));
                        fParticleTreeList[fNumberOfPrimaries].particles[0].x.emplace_back(particle.Vx(k));
                        fParticleTreeList[fNumberOfPrimaries].particles[0].y.emplace_back(particle.Vy(k));
                        fParticleTreeList[fNumberOfPrimaries].particles[0].z.emplace_back(particle.Vz(k));
                    }
                }
                else
                {
                    // find parent tree
                    Int_t parent_tree = findParentTree(fParticleTreeList, particle.Mother());
                    // check to see if this interaction is in the parent trajectory
                    Int_t traj_point = findParentLocation(
                        fParticleTreeList[parent_tree],
                        particle.Mother(),
                        particle.T(0),
                        particle.Vx(0),
                        particle.Vy(0),
                        particle.Vz(0)
                    );
                    // if not, then we need to add it
                    if (traj_point == -1) {
                        Int_t traj = findParentTrajectory(fParticleTreeList, particle.Mother());
                        // find closest time to start time of particle
                        // this will return the location of t < T
                        std::vector<Double_t>::iterator insert_time = findClosestTime(fParticleTreeList[parent_tree].particles[traj],
                            particle.T(0)
                        );
                        fParticleTreeList[parent_tree].particles[traj].t.insert(insert_time+1,particle.T(0),1);
                        fParticleTreeList[parent_tree].particles[traj].x.insert(insert_time+1,particle.Vx(0),1);
                        fParticleTreeList[parent_tree].particles[traj].y.insert(insert_time+1,particle.Vy(0),1);
                        fParticleTreeList[parent_tree].particles[traj].z.insert(insert_time+1,particle.Vz(0),1);
                    } 
                    for (size_t k = 0; k < particle.NumberTrajectoryPoints(); k++)
                    {
                        fParticleTreeList[parent_tree].particles[-1].t.emplace_back(particle.T(k));
                        fParticleTreeList[parent_tree].particles[-1].x.emplace_back(particle.Vx(k));
                        fParticleTreeList[parent_tree].particles[-1].y.emplace_back(particle.Vy(k));
                        fParticleTreeList[parent_tree].particles[-1].z.emplace_back(particle.Vz(k));
                    }
                }
                // check if the particle is the right pdg code
                for (size_t j = 0; j < fPdgCodes.size(); j++)
                {
                    if (particle.PdgCode() == fPdgCodes[j])
                    {
                        // check if particle should be recorded
                        // either the particle is a primary (Mother() == 0)
                        // or its a daughter particle of a primary and
                        // fCollectDaughters has been set to True
                        // or fCollectAll has been set to True.
                        if (particle.Mother() == 0 ||
                            (fCollectDaughters && checkEventIds(eventList,particle.TrackId())) ||
                            fCollectAll
                        )
                        {
                            // first we collect the trackid of the particle
                            // we are keeping and the pdg code and parent ids
                            // as well as all of the trackid's of its daughters.
                            eventList.ids.emplace_back(particle.TrackId());
                            eventList.particle_pdgs.emplace_back(particle.PdgCode());
                            eventList.particle_parent_ids.emplace_back(particle.Mother());
                            std::vector<Int_t> daughters;
                            if (fCollectDaughters)
                            {
                                for (int k = 0; k < particle.NumberDaughters(); k++) {
                                    eventList.ids.emplace_back(particle.Daughter(k));
                                    daughters.emplace_back(particle.Daughter(k));
                                }
                            }
                            eventList.particle_daughters.emplace_back(daughters);
                            // now we loop over all trajectory points and add them to the 
                            // the event list
                            for (size_t k = 0; k < particle.NumberTrajectoryPoints(); k++)
                            {
                                DetectorVolume currentVolume = fGeometry->getVolume(
                                    particle.Vx(), particle.Vy(), particle.Vz()
                                );
                                if (currentVolume.material_name == "LAr")
                                {
                                    eventList.particle_ids.emplace_back(particle.TrackId());
                                    eventList.particle_x.emplace_back(particle.Vx(k));
                                    eventList.particle_y.emplace_back(particle.Vy(k));
                                    eventList.particle_z.emplace_back(particle.Vz(k));
                                    eventList.particle_edep_energy.emplace_back(-1.);
                                    eventList.particle_edep_num_electrons.emplace_back(-1);
                                }
                            }
                            // now we loop over the set of particle ids to 
                            // find the ancestors of each particle until we
                            // reach a primary.
                            std::vector<Int_t> ancestors;
                            Int_t mother = particle.Mother();
                            // search until we reach a primary
                            while (mother != 0)
                            {
                                // first record the mother
                                ancestors.emplace_back(mother);
                                // now look for the next mother
                                Int_t particle = getParent(particleParentList, mother);
                                mother = particle;
                            }
                            ancestors.emplace_back(mother);
                            eventList.particle_ancestors.emplace_back(ancestors);
                        }
                    }
                }
            }
        }
        // now we loop over all of the energy deposits and check to see
        // which particles correspond to each deposit.
        auto mcEnergyDeposit = event.getValidHandle<std::vector<sim::SimEnergyDeposit>>(fIonAndScintProducerLabel);
        if (mcEnergyDeposit.isValid())
        {
            for (auto energyDeposit : *mcEnergyDeposit)
            {
                if (checkEventIds(eventList, energyDeposit.TrackID()))
                {
                    // if the particle is an electron, find the closest ancestor which
                    // is of another type.
                    Int_t eventIndex = findTrackId(eventList, energyDeposit.TrackID());
                    if (eventList.particle_pdgs[eventIndex] == 11)
                    {

                    }
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
        fTempEventList.event_id = eventList.event_id;
        fTempEventList.ids = eventList.ids;
        fTempEventList.particle_daughters = eventList.particle_daughters;
        fTempEventList.particle_pdgs = eventList.particle_pdgs;
        fTempEventList.particle_ids = eventList.particle_ids;
        fTempEventList.particle_parent_ids = eventList.particle_parent_ids;
        fTempEventList.particle_x = eventList.particle_x;
        fTempEventList.particle_y = eventList.particle_y;
        fTempEventList.particle_z = eventList.particle_z;
        fTempEventList.particle_edep_energy = eventList.particle_edep_energy;
        fTempEventList.particle_edep_num_electrons = eventList.particle_edep_num_electrons;
        fEventTree->Fill();

        // for (size_t k = 0; k < fParticleTreeList.size(); k++)
        // {
        //     fTempParticleTree.event_id = fParticleTreeList[k].event_id;
        //     fTempParticleTree.track_id = fParticleTreeList[k].track_id;
        //     fTempParticleTree.mother = fParticleTreeList[k].mother;
        //     fTempParticleTree.pdg = fParticleTreeList[k].pdg;
        //     fTempParticleTree.t = fParticleTreeList[k].t;
        //     fTempParticleTree.x = fParticleTreeList[k].x;
        //     fTempParticleTree.y = fParticleTreeList[k].y;
        //     fTempParticleTree.z = fParticleTreeList[k].z;
        //     fTempParticleTree.edep_energy = fParticleTreeList[k].edep_energy;
        //     fTempParticleTree.edep_num_electrons = fParticleTreeList[k].edep_num_electrons;
        //     fTempParticleTree.edge_start = fParticleTreeList[k].edge_start;
        //     fTempParticleTree.edge_end = fParticleTreeList[k].edge_end;
        //     fParticleTree->Fill();
        // }

        fTempParticleParentList.tracks = particleParentList.tracks;
        fTempParticleParentList.mothers = particleParentList.mothers;
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