#include "EventDisplay/inc/TextSelect.hh"
#include <iostream>

namespace mu2e {

// Setter: Called by the REve thread (EventDisplayManager)
void TextSelect::set(int run, int event) {
    // Lock the mutex for writing
    std::lock_guard<std::mutex> lock(_mutex);
    runN = run;
    eventN = event;
    std::cout<<"user has set "<<runN << " " << eventN <<std::endl;
    std::cout << "[TextSelect::set] Run/Event set to: " << runN << "/" << eventN << std::endl;
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
std::pair<int, int> TextSelect::getRunEvent() {
    // Lock the mutex for reading
    std::lock_guard<std::mutex> lock(_mutex);
    std::cout<<"user has received "<<runN << " " << eventN <<std::endl;
    return {runN, eventN};
}

}
