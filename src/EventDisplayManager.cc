#include "EventDisplay/inc/EventDisplayManager.hh"


namespace mu2e {

  EventDisplayManager::EventDisplayManager(
    ROOT::Experimental::REveManager* eveMgr,
    std::condition_variable& cv,
    std::mutex& m,
    GUI *fGui,
    TextSelect *fText)
    : REveElement{"EventManager"}, eveMng_{eveMgr}, cv_{&cv}, m_{&m}, fGui_(fGui), fText_(fText)
  {}

  void
  EventDisplayManager::NextEvent()
  {
    std::unique_lock lock{*m_};
    cv_->notify_all();
  }

  void
  EventDisplayManager::QuitRoot()
  {
    std::cout<<"Exit Signal 15, leaving REve Display "<<std::endl;
    exit(15);
  }
  void EventDisplayManager::autoplay(int  x)
   {
      std::cout << "EventManger autoplay() ....... " << x << std::endl;
   }
   int EventDisplayManager::getR(){
    return this->run;
  }
   void EventDisplayManager::setR(int runId){
    std::cout<<"[EventDisplayManager::setR] taking run number"<<this->run<<" and passing in "<<runId<<std::endl;
    this->run = runId;
    std::cout<<"[EventDisplayManager::setR] run number now "<<this->run<<std::endl;
  }
  
   void EventDisplayManager::goToRunEvent(int runId, int eventId)
   {
      std::cout << "[EventDisplayManager::goToRunEvent] received: " << runId<<" "<<eventId << std::endl;
      setR(runId);
      std::cout<<"[EventDisplayManager::goToRunEvent] value is now"<<this->run<<std::endl;
      fText_->set(runId,eventId);
   }

}
