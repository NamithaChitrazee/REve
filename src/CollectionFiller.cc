#include <TObject.h>
#include <TSystem.h>
#include <TFile.h>
#include "REve/inc/CollectionFiller.hh"
#include "art/Framework/Principal/SubRun.h"

using namespace mu2e;
namespace mu2e{

    CollectionFiller::CollectionFiller(const Config& conf) :
        chTag_(conf.chTag()),
        tcTag_(conf.tcTag()),
        crvcoinTag_(conf.crvdigiTag()),
        cluTag_(conf.cluTag()),
        cryHitTag_(conf.cryHitTag()),
        kalSeedTag_(conf.kalSeedTag()),
        cosmicTrackSeedTag_(conf.cosmicTrackSeedTag()),
        MCTrajTag_(conf.MCTrajTag()),
        addHits_(conf.addHits()),
        addCrvHits_(conf.addCrvHits()),
        addTimeClusters_(conf.addTimeClusters()),
        addTrkHits_(conf.addTrkHits()),
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
            for(const auto &tag : chTag_){
                auto chH = evt.getValidHandle<mu2e::ComboHitCollection>(tag);
                data.chcol = chH.product();
                data.combohit_list.push_back(data.chcol);
                std::string name = TurnNameToString(tag);
                std::cout<<"Plotting ComboHit Instance: "<<name<<std::endl;
                data.combohit_labels.push_back(name);
            }
            data.combohit_tuple = std::make_tuple(data.combohit_labels,data.combohit_list);
        }
         if(FillAll_ or (addCrvHits_ and CollectionName==CRVRecoPulses)){
           auto chH = evt.getValidHandle<mu2e::CrvRecoPulseCollection>(crvcoinTag_);
           data.crvcoincol = chH.product();
         }
        
        if(FillAll_  or (CollectionName == TimeClusters)){ 
           auto chH = evt.getValidHandle<mu2e::TimeClusterCollection>(tcTag_);
           data.tccol = chH.product();
        }
        if(FillAll_ or (addTrkHits_ and CollectionName == TrkHits)){ 

           for(const auto &tag : chTag_){
                auto chH = evt.getValidHandle<mu2e::ComboHitCollection>(tag);
                data.chcol = chH.product();
                data.combohit_list.push_back(data.chcol);
                std::string name = TurnNameToString(tag);
                std::cout<<"Plotting ComboHit Instance: "<<name<<std::endl;
                data.combohit_labels.push_back(name);
            }
            data.combohit_tuple = std::make_tuple(data.combohit_labels,data.combohit_list);
        }
        if(FillAll_  or (CollectionName == CaloClusters)){
            for(const auto &tag : cluTag_){
                auto chH = evt.getValidHandle<mu2e::CaloClusterCollection>(tag);
                data.clustercol = chH.product();
                data.calocluster_list.push_back(data.clustercol);
                std::string name = TurnNameToString(tag);
                std::cout<<"Plotting CaloCluster Instance: "<<name<<std::endl;
                data.calocluster_labels.push_back(name);
            }
            data.calocluster_tuple = std::make_tuple(data.calocluster_labels,data.calocluster_list);
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
            
        }
        if(FillAll_  or (CollectionName == CosmicTrackSeeds)){
            auto chH = evt.getValidHandle<mu2e::CosmicTrackSeedCollection>(cosmicTrackSeedTag_);
            data.CosmicTrackSeedcol = chH.product();
        }
    }
    
    void CollectionFiller::FillMCCollections(const art::Event& evt, DataCollections &data, MCDataProductName CollectionName){

        if(FillAll_ or (CollectionName==MCTrajectories)){
    
                for(const auto &tag : MCTrajTag_){
                  auto chH = evt.getValidHandle<mu2e::MCTrajectoryCollection>(tag);
                  data.mctrajcol = chH.product();
                  data.mctrack_list.push_back(data.mctrajcol);

                  std::string name = TurnNameToString(tag);
                  std::cout<<"Plotting MCTrajectory Instance: "<<name<<std::endl;
                  data.mctrack_labels.push_back(name);

            }
            data.mctrack_tuple = std::make_tuple(data.mctrack_labels,data.mctrack_list);
         
                }
        
    }

}
