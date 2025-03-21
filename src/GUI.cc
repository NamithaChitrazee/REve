#include "EventDisplay/inc/GUI.hh"

using namespace mu2e;


int GUI::WriteCoreJson(nlohmann::json &j, int rnr_offset)
{
  j["path"] = "Event/SubRun/Run"; 
  j["eventid"] = feventid;
  j["subrunid"] = fsubrunid;
  j["runid"] = frunid;
  j["UT_PostStream"] = "UT_refresh_event_info";
  return ROOT::Experimental::REveElement::WriteCoreJson(j, 0);
}
void GUI::setEvent(int event){ eventNumber = event; }
void GUI::setRun(int run){ runNumber = run; }
void GUI::getRunEvent(){ std::cout<<"GUI SETTTTT"<<runNumber<<" "<<eventNumber<<std::endl;}


  
