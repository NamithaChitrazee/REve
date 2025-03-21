#ifndef _GeomUtil_hh
#define _GeomUtil_hh

#include <ROOT/REveElement.hxx>
#include <ROOT/REveManager.hxx>

#include "Offline/StoppingTargetGeom/inc/StoppingTarget.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"

namespace REX = ROOT::Experimental;

namespace mu2e {
    class GeomUtil  : public REX::REveElement {

        public :
            explicit GeomUtil() { SetErrorHandler(DefaultErrorHandler); }
            virtual ~GeomUtil() {}
            #ifndef __CINT__
            double FindStoppingTarget_z();
            double GetStoppingTarget_z(){ return StoppingTarget_z; }
            double StoppingTarget_z;

            #endif
            ClassDef(GeomUtil, 0);



      };

}
#endif
