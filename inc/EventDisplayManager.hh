#ifndef _EventDisplayManager_hh
#define _EventDisplayManager_hh

#include <ROOT/REveElement.hxx>
#include <ROOT/REveScene.hxx>
#include <condition_variable>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <iostream>
#include "EventDisplay/inc/GUI.hh"
#include "EventDisplay/inc/TextSelect.hh"
#include "nlohmann/json.hpp"
#include <ROOT/REveManager.hxx>

namespace ROOT::Experimental {
  class REveManager;
}

namespace mu2e {

  constexpr auto invalid_event = std::numeric_limits<unsigned>::max();

  class EventDisplayManager : public ROOT::Experimental::REveElement {
    public:
      EventDisplayManager() = default; // ROOT needs a dictionary

      explicit EventDisplayManager(ROOT::Experimental::REveManager* eveMgr,
                                       std::condition_variable& cv,
                                       std::mutex& m,
                                       GUI *fGui);

    void NextEvent();
    void QuitRoot();
    void autoplay(int x);
    int getR();
    void setR(int runId);
    void goToRunEvent(int runId, int eventId);
    int run{0};
    void setTextSelect(TextSelect* fText);
    std::uint32_t fTextId_{0};
    void setTextSelectId(std::uint32_t textId);
    private:
      ROOT::Experimental::REveManager* eveMng_{nullptr};
      std::condition_variable* cv_{nullptr};
      std::mutex* m_{nullptr};
      bool doneProcessingEvents_{false};
      GUI *fGui_{nullptr};
      TextSelect *fText_{nullptr};
    };
}

#endif /* EventDisplayManager_h */
