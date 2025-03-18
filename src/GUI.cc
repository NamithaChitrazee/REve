#include "Mu2eEventDisplay/inc/REveMu2eGUI.hh"

using namespace mu2e;


int REveMu2eGUI::WriteCoreJson(nlohmann::json &j, int rnr_offset)
{
  j["path"] = "Event/SubRun/Run"; 
  j["eventid"] = feventid;
  j["subrunid"] = fsubrunid;
  j["runid"] = frunid;
  j["UT_PostStream"] = "UT_refresh_event_info";
  return ROOT::Experimental::REveElement::WriteCoreJson(j, 0);
}
void REveMu2eGUI::setEvent(int event){ eventNumber = event; }
void REveMu2eGUI::setRun(int run){ runNumber = run; }
void REveMu2eGUI::getRunEvent(){ std::cout<<"GUI SETTTTT"<<runNumber<<" "<<eventNumber<<std::endl;}


  
