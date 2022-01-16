#ifndef CollectionFiller_hh
#define CollectionFiller_hh

//Art:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include <TObject.h>
#include <TROOT.h>
#include <TGComboBox.h>
#include <TGListBox.h>

#include "REve/inc/DataCollections.hh"

namespace mu2e{

  enum RecoDataProductName {ComboHits, CRVRecoPulses, TimeClusters, CaloClusters, KalSeeds, CosmicTrackSeeds, TrkHits};
  enum MCDataProductName {MCTrajectories};
	class CollectionFiller
	{
      public:
        struct Config{
          using Name=fhicl::Name;
          using Comment=fhicl::Comment;
          fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("for info"),0};
          fhicl::Sequence<art::InputTag>chTag{Name("ComboHitCollection"),Comment("chTag")};
          fhicl::Atom<art::InputTag>tcTag{Name("TimeClusterCollection"),Comment("ttcTag"), "makePH"};
		fhicl::Sequence<art::InputTag>crvdigiTag{Name("CrvRecoPulseCollection"),Comment("crvTag")};
          fhicl::Sequence<art::InputTag>cluTag{Name("CaloClusterCollection"),Comment("cluTag")};
		 fhicl::Sequence<art::InputTag>cryHitTag{Name("CaloHitCollection"),Comment("cryHitTag")};
          fhicl::Sequence<art::InputTag>kalSeedTag{Name("KalSeedCollection"),Comment("kalseedTag")};
          fhicl::Atom<art::InputTag>cosmicTrackSeedTag{Name("CosmicTrackSeedCollection"),Comment("cosmicTrackSeedTag")};
          fhicl::Sequence<art::InputTag>MCTrajTag{Name("MCTrajectoryCollection"),Comment("MCTrajTag")};
		
          fhicl::Atom<bool> addHits{Name("addHits"), Comment("set to add the hits"),false};
          fhicl::Atom<bool> addCrvHits{Name("addCrvHits"), Comment("set to add crv hits"),false};
		fhicl::Atom<bool> addTimeClusters{Name("addTimeClusters"), Comment("set to add the TC hits"),false};
	        fhicl::Atom<bool> addTrkHits{Name("addTrkHits"), Comment("set to add the Trk hits"),false};
          fhicl::Atom<bool> addClusters{Name("addClusters"), Comment("set to add caloclusters"),false};
          fhicl::Atom<bool> addKalSeeds{Name("addKalSeeds"), Comment("set to add kalseeds"),false};
          fhicl::Atom<bool> addCosmicTrackSeeds{Name("addCosmicTrackSeeds"), Comment("set to add cosmic track seeds"),false};
          fhicl::Atom<bool> addMCTraj{Name("addMCTraj"), Comment("set to add add MC information"),false};
          fhicl::Atom<bool> FillAll{Name("FillAll"), Comment("to see all available products"), false};
        };

        explicit CollectionFiller(const Config& conf);
        CollectionFiller(const CollectionFiller &);
        CollectionFiller& operator=(const CollectionFiller &);

        std::vector<art::InputTag> chTag_;
        art::InputTag tcTag_; //TODO - make in same style as others
	std::vector<art::InputTag> crvcoinTag_;
        std::vector<art::InputTag> cluTag_;
	std::vector<art::InputTag> cryHitTag_;
        std::vector<art::InputTag> kalSeedTag_;
        art::InputTag cosmicTrackSeedTag_;
        std::vector<art::InputTag> MCTrajTag_;
        art::Event *_event;
        art::Run *_run;
        bool addHits_,  addCrvHits_, addTimeClusters_, addTrkHits_, addClusters_, addKalSeeds_, addCosmicTrackSeeds_, addMCTraj_, FillAll_;
        void FillRecoCollections(const art::Event& evt, DataCollections &data, RecoDataProductName code);
        void FillMCCollections(const art::Event& evt, DataCollections &data, MCDataProductName code);
        virtual ~CollectionFiller(){};

      private:
        Config _conf;
      
};

}

#endif 
