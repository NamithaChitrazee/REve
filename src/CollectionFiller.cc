#include <TObject.h>
#include <TSystem.h>
#include <TFile.h>
#include "REve/inc/CollectionFiller.hh"
#include "art/Framework/Principal/SubRun.h"

using namespace mu2e;
namespace mu2e{

    CollectionFiller::CollectionFiller(const Config& conf) :
        chTag_(conf.chTag()),
        cluTag_(conf.cluTag()),
        kalSeedTag_(conf.kalSeedTag()),
        cosmicTrackSeedTag_(conf.cosmicTrackSeedTag()),
        MCTrajTag_(conf.MCTrajTag()),
        addHits_(conf.addHits()),
        addClusters_(conf.addClusters()),
        addKalSeeds_(conf.addKalSeeds()),
        addCosmicTrackSeeds_(conf.addCosmicTrackSeeds()),
        addMCTraj_(conf.addMCTraj()),
        FillAll_(conf.FillAll())
        {}

    template <typename T> std::string TurnNameToString( const T& value )
    {
        std::ostringstream ss;
        ss << value;
        return ss.str();
    }
  

    void CollectionFiller::FillRecoCollections(const art::Event& evt, DataCollections &data, RecoDataProductName CollectionName){
        if(FillAll_  or (CollectionName == ComboHits)){ 
            auto chH = evt.getValidHandle<mu2e::ComboHitCollection>(chTag_);
            data.chcol = chH.product();
        }
        if(FillAll_  or (CollectionName == CaloClusters)){
            auto chH = evt.getValidHandle<mu2e::CaloClusterCollection>(cluTag_);
            data.clustercol = chH.product();
        }
        if(FillAll_ or (CollectionName==KalSeeds)){
          
            for(const auto &tag : kalSeedTag_){
                auto chH = evt.getValidHandle<mu2e::KalSeedCollection>(tag);
                data.kalSeedcol = chH.product();
                data.track_list.push_back(data.kalSeedcol);

                std::string name = TurnNameToString(tag);
                std::cout<<"Plotting KalSeed Instance: "<<name<<std::endl;
                data.track_labels.push_back(name);

            }
            data.track_tuple = std::make_tuple(data.track_labels,data.track_list);
            std::cout<<"size "<<data.track_list.size()<<std::endl;
        }
        if(FillAll_  or (CollectionName == CosmicTrackSeeds)){
            auto chH = evt.getValidHandle<mu2e::CosmicTrackSeedCollection>(cosmicTrackSeedTag_);
            data.CosmicTrackSeedcol = chH.product();
        }
    }
    
    void CollectionFiller::FillMCCollections(const art::Event& evt, DataCollections &data, MCDataProductName CollectionName){

        if(FillAll_ or (CollectionName==MCTrajectories)){
          
            for(const auto &tag : kalSeedTag_){
                auto chH = evt.getValidHandle<mu2e::MCTrajectoryCollection>(tag);
                data.mctrajcol = chH.product();
                data.mctrack_list.push_back(data.mctrajcol);

                std::string name = TurnNameToString(tag);
                std::cout<<"Plotting MCTrajectory Instance: "<<name<<std::endl;
                data.mctrack_labels.push_back(name);

            }
            data.mctrack_tuple = std::make_tuple(data.mctrack_labels,data.mctrack_list);
            std::cout<<"size "<<data.track_list.size()<<std::endl;
        }
        
    }

}


