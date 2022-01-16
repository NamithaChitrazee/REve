//Author: S Middleton
//Date: 2021
//Purpose: REveEventDisplay Driving module

#include "TRACE/trace.h"
//Art:
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/ParameterSetRegistry.h"

#include "artdaq-core/Data/ContainerFragment.hh"
#include "artdaq-core/Data/Fragment.hh"
//#include "artdaq/DAQdata/Globals.hh"

#include "cetlib_except/exception.h"

//ROOT:
//#include "art_root_io/TFileService.h" 
#include <TApplication.h>
#include <TSystem.h>
#include <TList.h>
#include <TObjArray.h>
#include <Rtypes.h>
#include <TFile.h>
#include <TTree.h>

//EVE-7
#include <ROOT/RWebWindow.hxx>
#include <ROOT/RWebWindowsManager.hxx>
#include <ROOT/REveManager.hxx>

#include <condition_variable>
#include <iostream>
#include <mutex>
#include <thread>

//Offline:
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wtype-limits" 
#pragma GCC diagnostic ignored "-Wpedantic"
#pragma GCC diagnostic pop

//EveMu2e
#include "REve/inc/REveMainWindow.hh"
#include "REve/inc/EventDisplayManager.hh"
#include "REve/inc/CollectionFiller.hh"
#include "REve/inc/DataCollections.hh"

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"


using namespace std;
using namespace mu2e;
namespace mu2e
{
    class XThreadTimer : public TTimer {
        std::function<void()> foo_;
        public:
            XThreadTimer(std::function<void()> f) : foo_(f)
            {
                SetTime(0);
                R__LOCKGUARD2(gSystemMutex);
                gSystem->AddTimer(this);
            }
            Bool_t Notify() override
            {
                foo_();
                gSystem->RemoveTimer(this);
                return kTRUE;
            }
    };

    class REveEventDisplay : public art::EDAnalyzer {
      public:
        struct Config{
          using Name=fhicl::Name;
          using Comment=fhicl::Comment;
          fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("for info"),0};
          fhicl::Atom<bool> showCRV{Name("showCRV"), Comment("set false if you just want to see DS"),false};   
          fhicl::Atom<bool> show2D{Name("show2D"), Comment(""),true};   
          fhicl::Table<CollectionFiller::Config> filler{Name("filler"),Comment("fill collections")};
          fhicl::Sequence<int>particles{Name("particles"),Comment("PDGcodes to plot")};
        };

        typedef art::EDAnalyzer::Table<Config> Parameters;
        explicit REveEventDisplay(const Parameters& conf);
        virtual ~REveEventDisplay();
        virtual void beginJob() override;
        virtual void beginRun(const art::Run& run) override;
        virtual void analyze(const art::Event& e);
        virtual void endJob() override;
        void signalAppStart();
    private:

        art::ServiceHandle<art::TFileService> tfs;
        Config _conf;
        bool showCRV_;
        
        void setup_eve();
        void run_application();
        void process_single_event();
        void printOpts();

        // Application control
        TApplication application_{"none", nullptr, nullptr};
        std::thread appThread_{};

        // Display control
        art::EventID displayedEventID_{};
        REX::REveManager* eveMng_{nullptr};
        EventDisplayManager* eventMgr_{nullptr};

        // Control between the main thread and event-display thread
        std::condition_variable cv_{};
        std::mutex m_{};

        TDirectory*   directory_ = nullptr;   
        CollectionFiller filler_;
        REveMainWindow *frame_;
        DataCollections data;
        bool firstLoop_ = true; 
        std::vector<int> particles_;
        
    };


  REveEventDisplay::REveEventDisplay(const Parameters& conf)  :
    art::EDAnalyzer(conf),
    showCRV_(conf().showCRV()),
    filler_(conf().filler()),
    particles_(conf().particles())
    {}

  REveEventDisplay::~REveEventDisplay() {}

  void REveEventDisplay::signalAppStart()
  {
      std::unique_lock lock{m_};
      cv_.notify_all();
  }



  void REveEventDisplay::beginJob(){
      std::cout<<"[REveEventDisplay : beginJob()] -- starting ..."<<std::endl;
      {
      std::unique_lock lock{m_};

      appThread_ = std::thread{[this] { run_application(); }};

      // Wait for app init to finish ... this will process pending timer events.
      XThreadTimer sut([this]{ signalAppStart(); });
      std::cout<<"[REveEventDisplay : beginJob()] -- starting wait on app start"<<std::endl;
      cv_.wait(lock);
      std::cout<<"[REveEventDisplay : beginJob()] -- app start signal received, starting eve init"<<std::endl;
      XThreadTimer suet([this]{ setup_eve(); });
      std::cout<<"[REveEventDisplay : beginJob()] -- starting wait on eve setup"<<std::endl;
      cv_.wait(lock);
      std::cout<<"[REveEventDisplay : beginJob()] -- eve setup apparently complete"<<std::endl;
      }
  }


  void REveEventDisplay::beginRun(const art::Run&){


  }

  void REveEventDisplay::printOpts(){
    std::cout<<"*********** REve Mu2e **************"
    <<" User Options: "
    <<" addHits : "<<filler_.addHits_
    <<" addTimeClusters : "<<filler_.addTimeClusters_
    <<" addCRVpulses : "<<filler_.addCrvHits_
    <<" addClusters : "<<filler_.addClusters_
    <<" addTracks : "<<filler_.addKalSeeds_
    <<" addCosmicTrackSeeds : "<<filler_.addCosmicTrackSeeds_<<std::endl;
  }
  
  
  
  void REveEventDisplay::analyze(art::Event const& event){
      std::cout<<"starting analyze"<<std::endl;
      //remove previous event objects
      data.Reset();
      
      // Update state relevant for displaying new event.
      displayedEventID_ = event.id();

      // Hand off control to display thread
      std::unique_lock lock{m_};
      std::cout<<"[REveEventDisplay : analyze()] -- Fill collections "<<std::endl;
      if(filler_.addClusters_)  filler_.FillRecoCollections(event, data, CaloClusters);
      if(filler_.addHits_)  filler_.FillRecoCollections(event, data, ComboHits);
      if(filler_.addCrvHits_) filler_.FillRecoCollections(event, data, CRVRecoPulses);
      if(filler_.addTimeClusters_) filler_.FillRecoCollections(event, data, TimeClusters);
      if(filler_.addTrkHits_) filler_.FillRecoCollections(event, data, TrkHits); 
      if(filler_.addKalSeeds_)  filler_.FillRecoCollections(event, data, KalSeeds);
      if(filler_.addCosmicTrackSeeds_)  filler_.FillRecoCollections(event, data, CosmicTrackSeeds);
      if(filler_.addMCTraj_)  filler_.FillMCCollections(event, data, MCTrajectories);
     
      std::cout<<"[REveEventDisplay : analyze()] -- Event processing started "<<std::endl;
      XThreadTimer proc_timer([this]{ process_single_event(); });
      std::cout<<"[REveEventDisplay : analyze()] -- transferring to TApplication thread "<<std::endl;
      cv_.wait(lock);
      std::cout<<"[REveEventDisplay : analyze()] -- TApplication thread returning control "<<std::endl;
      std::cout<<"[REveEventDisplay : analyze()] Ended Event "<<std::endl; 
  }

    void REveEventDisplay::endJob()
    {
      std::cout<<"[REveEventDisplay : EndJob] Start "<<std::endl;
      application_.Terminate(0);

      if (appThread_.joinable()) {
        appThread_.join();
      }
      std::cout<<"[REveEventDisplay : EndJob] End "<<std::endl;
    }



  // Functions invoked by the threads explicitly created by this module
  void REveEventDisplay::run_application()
  {
      // Without this the startup timer might not get invoked.
      gSystem->ProcessEvents();
      application_.Run(true);
  }


  void REveEventDisplay::setup_eve()
  {
      REX::RWebWindowsManager::AssignMainThrd();
      eveMng_ = REX::REveManager::Create();

      eventMgr_ = new EventDisplayManager{eveMng_, cv_, m_};
      auto world = eveMng_->GetWorld();
      assert(world);
      world->AddElement(eventMgr_);
      world->AddCommand("QuitRoot",  "sap-icon://log",  eventMgr_, "QuitRoot()");
      world->AddCommand("NextEvent", "sap-icon://step", eventMgr_, "NextEvent()");
      frame_ = new REveMainWindow();
      frame_->makeGeometryScene(eveMng_,showCRV_);

      std::unique_lock lock{m_};
      cv_.notify_all();
 
  }

  // Actually interesting function responsible for drawing the current event
  void REveEventDisplay::process_single_event()
    { 
      std::cout<<"[REveEventDisplay : process_single_event] Start "<<std::endl;
      eveMng_->DisableRedraw();
      eveMng_->GetWorld()->BeginAcceptingChanges();
      eveMng_->GetScenes()->AcceptChanges(true);
      std::cout<<"[REveEventDisplay : process_single_event] -- extract event scene "<<std::endl;
      REX::REveElement* scene = eveMng_->GetEventScene();

      std::cout<<"[REveEventDisplay : process_single_event] -- calls to data interface "<<std::endl;
      DrawOptions drawOpts(false, filler_.addCosmicTrackSeeds_, filler_.addKalSeeds_, filler_.addClusters_, filler_.addHits_, filler_.addCrvHits_, filler_.addTimeClusters_, filler_.addTrkHits_, filler_.addMCTraj_);
      frame_->showEvents(eveMng_, scene, firstLoop_, data, drawOpts, particles_);

      std::cout<<"[REveEventDisplay : process_single_event] -- cluster added to scene "<<std::endl;
      firstLoop_ = false;
      eveMng_->GetScenes()->AcceptChanges(false);
      eveMng_->GetWorld()->EndAcceptingChanges();
      eveMng_->EnableRedraw();
      
      std::cout<<"[REveEventDisplay : process_single_event] End "<<std::endl;
    }
  }
  using mu2e::REveEventDisplay;
  DEFINE_ART_MODULE(REveEventDisplay)
