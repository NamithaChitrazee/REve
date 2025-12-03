#ifndef _TextSelect_hh
#define _TextSelect_hh

#include <ROOT/REveElement.hxx>
#include <mutex>

namespace REX = ROOT::Experimental;
using namespace ROOT::Experimental;

namespace mu2e {
class TextSelect : public ROOT::Experimental::REveElement
{
  public:
    // Constructor remains implicit or default
    TextSelect() : REveElement{"TextSelect"} {} 
      
    // Setter (Called by the GUI thread/EventDisplayManager)
    void set(int run, int event);

    // Getter (Called by the art module thread)
    // We now retrieve both run and event in one call for atomic safety
    std::pair<int, int> getRunEvent(); // <<< Modified getter signature

    // setAutoplay
    void setAutoplay(int x);
    int getAutoplay();
    private:
     int runN = 0;
     int eventN = 0;
     int autoplay_=0;
     std::mutex _mutex; // <<< The synchronization tool
};
}

#endif
