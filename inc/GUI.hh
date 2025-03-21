#ifndef _GUI_hh
#define _GUI_hh

#include <ROOT/REveElement.hxx>
#include "nlohmann/json.hpp"
namespace REX = ROOT::Experimental;
using namespace ROOT::Experimental;

namespace mu2e {
class GUI : public ROOT::Experimental::REveElement
{
  public:
     int feventid{0};
     int fsubrunid{0};
     int frunid{0};
     int WriteCoreJson(nlohmann::json &j, int rnr_offset) override;
     void setRun(int run);
     void setEvent(int event);
     void getRunEvent();
     void PrintEventInfo();
     int runNumber{0};
     int eventNumber{0};

};
}

#endif
