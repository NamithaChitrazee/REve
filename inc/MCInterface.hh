#ifndef REveMu2eMCInterface_hh
#define REveMu2eMCInterface_hh
#include <ROOT/REveElement.hxx>
#include <ROOT/REvePointSet.hxx>
#include <ROOT/REveManager.hxx>
#include <ROOT/REveLine.hxx>
#include <ROOT/REveScene.hxx>
#include <ROOT/REveCompound.hxx>
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryPoint.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include <TApplication.h>
#include <TEvePad.h>
#include <TObject.h>
#include <TSystem.h>
#include <vector>
#include <tuple>
#include <string>


using namespace mu2e;
namespace REX = ROOT::Experimental;
namespace mu2e{
    class REveMu2eMCInterface {
        public:
          static int const mstyle = 20 ;
          static int const msize = 5;
          explicit REveMu2eMCInterface(){};
          explicit REveMu2eMCInterface(const REveMu2eMCInterface &);
          REveMu2eMCInterface& operator=(const REveMu2eMCInterface &);
          virtual ~REveMu2eMCInterface() = default;
          #ifndef __CINT__
          inline constexpr double pointmmTocm(double mm){ return mm/10; };
          int Contains(std::vector<int> v, int x);
          const char* GetParticleName(int PDGCode);
          void SetLineColorPID(int PDGCode,REX::REveLine *line);
          void toExtracted(CLHEP::Hep3Vector& Pos);
          void AddMCTrajectoryCollection(REX::REveManager *&eveMng,bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const MCTrajectoryCollection*>> track_tuple, REX::REveElement* &scene, std::vector<int> particles, bool extracted );
          void AddSurfaceStepCollection(REX::REveManager *&eveMng,bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const SurfaceStepCollection*>> track_tuple, REX::REveElement* &scene, std::vector<int> particles, bool extracted );
          #endif
          ClassDef(REveMu2eMCInterface, 0);
      };
}


#endif
