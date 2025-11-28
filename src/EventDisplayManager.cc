#include "EventDisplay/inc/EventDisplayManager.hh"

namespace mu2e {

  EventDisplayManager::EventDisplayManager(
    ROOT::Experimental::REveManager* eveMgr,
    std::condition_variable& cv,
    std::mutex& m,
    GUI *fGui) // <<< TextSelect *fText IS REMOVED
    : 
    REveElement{"EventManager"}, 
    eveMng_{eveMgr}, 
    cv_{&cv}, 
    m_{&m}, 
    fGui_(fGui) // <<< fText_(fText) IS REMOVED
{}

// Implement the setter method
void EventDisplayManager::setTextSelect(TextSelect* fText) {
    fText_ = fText;
    if (fText_ == nullptr) {
        std::cerr << "INTERNAL ERROR: setTextSelect called with a NULL pointer." << std::endl;
    } else {
        std::cout << "[EventDisplayManager::setTextSelect] fText_ successfully set." << std::endl;
    }
}
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
  // Add a function to store the ID
void EventDisplayManager::setTextSelectId(std::uint32_t textId) {
    fTextId_ = textId;
    std::cout << "[EventDisplayManager::setTextSelectId] fTextId_ set to: " << fTextId_ << std::endl;
}
void mu2e::EventDisplayManager::goToRunEvent(int runId, int eventId)
{
    std::cout << "[EventDisplayManager::goToRunEvent] received: " << runId<<" "<<eventId << std::endl;
    
    TextSelect* fText_obj = nullptr;
    
    if (ROOT::Experimental::gEve != nullptr) {
        
        ROOT::Experimental::REveElement* element = ROOT::Experimental::gEve->FindElementById(4285);//fTextId_); //FIXME note the hardcoding
        
        if (element != nullptr) {
            fText_obj = dynamic_cast<TextSelect*>(element);
        }
    }
    
    if (fText_obj == nullptr) {
        std::cerr << "CRITICAL ERROR: TextSelect object not found via gEve->FindElementWithId(" 
                  << fTextId_ << "). Cannot set Run/Event." << std::endl; 
        return;  
    }

    fText_obj->set(runId, eventId); 
}
}

