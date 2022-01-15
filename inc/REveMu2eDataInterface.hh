#ifndef REveMu2eDataInterface_hh
#define REveMu2eDataInterface_hh
#include <ROOT/REveElement.hxx>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveLine.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveTrackPropagator.hxx>
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/CosmicRayShieldGeom/inc/CosmicRayShield.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include <TApplication.h>
#include <TEvePad.h>
#include <TObject.h>
#include <TSystem.h>
#include <limits>
#include <vector>
#include <tuple>
#include <algorithm>

using namespace mu2e;
namespace REX = ROOT::Experimental;
namespace mu2e{
    class REveMu2eDataInterface {
        public:
          static int const mstyle = 4;
          static int const msize=6;
          explicit REveMu2eDataInterface(){};
          explicit REveMu2eDataInterface(const REveMu2eDataInterface &);
          REveMu2eDataInterface& operator=(const REveMu2eDataInterface &);
          virtual ~REveMu2eDataInterface() = default;
          #ifndef __CINT__
      
          void AddComboHits(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const ComboHitCollection*>> combohit_tuple, REX::REveElement* &scene);
          //TODO - place these following functions in same style as others.
          void AddTimeClusters(REX::REveManager *&eveMng, bool firstloop, const mu2e::TimeClusterCollection *tccol, REX::REveElement* &scene);
          void AddTrkHits(REX::REveManager *&eveMng, bool firstLoop_, const mu2e::ComboHitCollection *chcol, std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, REX::REveElement* &scene);
          void AddCaloClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CaloClusterCollection*>> calocluster_tuple, REX::REveElement* &scene);
          void AddCRVInfo(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CrvRecoPulseCollection*>> crvpulse_tuple, REX::REveElement* &scene);
          void AddKalSeedCollection(REX::REveManager *&eveMng,bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const KalSeedCollection*>> track_tuple, REX::REveElement* &scene);
          void AddCosmicTrackFit(REX::REveManager *&eveMng, bool firstLoop_, const mu2e::CosmicTrackSeedCollection *cosmiccol, REX::REveElement* &scene);
          unsigned int  maxHits = 10000;
          unsigned int  minHits = 0;
          #endif
          ClassDef(REveMu2eDataInterface, 0);
      }; 
}


#endif 


