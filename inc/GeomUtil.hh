#ifndef _REveMu2eGeomUtil_hh
#define _REveMu2eGeomUtil_hh

#include <ROOT/REveElement.hxx>
#include <ROOT/REveManager.hxx>

#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

namespace REX = ROOT::Experimental;

namespace mu2e {
    class REveMu2eGeomUtil  : public REX::REveElement {

        public :
            explicit REveMu2eGeomUtil() { SetErrorHandler(DefaultErrorHandler); }
            virtual ~REveMu2eGeomUtil() {}
            #ifndef __CINT__
            double FindStoppingTarget_z();
            double GetStoppingTarget_z(){ return StoppingTarget_z; }
            double StoppingTarget_z;

            #endif
            ClassDef(REveMu2eGeomUtil, 0);



      };

}
#endif
