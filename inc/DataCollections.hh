#ifndef _DataCollections_hh
#define _DataCollections_hh
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryPoint.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
//Art/FCL:
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include <TObject.h>
#include <TROOT.h>
#include <TGComboBox.h>
#include <TGListBox.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <string>

namespace mu2e{

	class DataCollections
	{
    public:
      explicit DataCollections(){};
      DataCollections(const DataCollections &){};
      DataCollections& operator=(const DataCollections &);

      //DataProducts:
      const mu2e::ComboHitCollection* chcol = 0;   
		const mu2e::CrvRecoPulseCollection* crvRecoPulse = 0;
      const TimeClusterCollection *tccol = 0;
      const mu2e::CaloClusterCollection* clustercol = 0;
      const mu2e::KalSeedCollection* kalSeedcol = 0;
      const mu2e::CosmicTrackSeedCollection* CosmicTrackSeedcol = 0;
      const mu2e::MCTrajectoryCollection *mctrajcol = 0;
      //lists:
      std::vector<const mu2e::KalSeedCollection*> track_list;
      std::vector<const mu2e::CaloClusterCollection*> calocluster_list;
      std::vector<const mu2e::ComboHitCollection*> combohit_list;
      std::vector<const mu2e::CRVRecoPulseCollection*> crvpulse_list;
      std::vector<const mu2e::MCTrajectoryCollection*> mctrack_list;
      //Input Tag Labels:
      std::vector<std::string> track_labels;
      std::vector<std::string> calocluster_labels;
      std::vector<std::string> mctrack_labels;
      std::vector<std::string> combohit_labels;
		std::vector<std::string> crvpulse_labels;
      //Link Labels and Lists:
      std::tuple<std::vector<std::string>, std::vector<const mu2e::KalSeedCollection*>> track_tuple;
      std::tuple<std::vector<std::string>, std::vector<const mu2e::CaloClusterCollection*>> calocluster_tuple;
      std::tuple<std::vector<std::string>, std::vector<const mu2e::ComboHitCollection*>> combohit_tuple;
      std::tuple<std::vector<std::string>, std::vector<const mu2e::CRVRecoPulseCollection*>> crvpulse_tuple;
      std::tuple<std::vector<std::string>, std::vector<const mu2e::MCTrajectoryCollection*>> mctrack_tuple;
      
      void Reset(){
        this->track_list.clear();
        this->calocluster_list.clear();
        this->combohit_list.clear();
	this->crvpulse_list.clear();
        this->mctrack_list.clear();
        this->track_labels.clear();
        this->calocluster_labels.clear();
        this->mctrack_labels.clear();
        this->combohit_labels.clear();
      }
      
      virtual ~DataCollections(){};

	};
}
#endif 
