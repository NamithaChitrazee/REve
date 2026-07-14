#include "EventDisplay/inc/TextSelect.hh"
#include <iostream>

namespace mu2e {

// Setter: Called by the REve thread (EventDisplayManager)
void TextSelect::set(int run, int subrun, int event) {
    std::lock_guard<std::mutex> lock(_mutex);
    runN    = run;
    subrunN = subrun;
    eventN  = event;
    std::cout << "[TextSelect::set] Run/Subrun/Event set to: "
              << runN << "/" << subrunN << "/" << eventN << std::endl;
}

// Setter: Called by the REve thread (EventDisplayManager)
void TextSelect::setAutoplay(int x) {
    // Lock the mutex for writing
    std::lock_guard<std::mutex> lock(_mutex);
    autoplay_ = x;
    std::cout << "[TextSelect::setAutoplay] set to: " << autoplay_ << std::endl;
}

int TextSelect::getAutoplay(){
  return autoplay_;
}

// Getter: Called by the art module thread (Mu2eEventDisplay::analyze)
std::tuple<int, int, int> TextSelect::getRunEvent() {
    std::lock_guard<std::mutex> lock(_mutex);
    return {runN, subrunN, eventN};
}

}
