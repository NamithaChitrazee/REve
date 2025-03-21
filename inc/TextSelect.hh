#ifndef _TextSelect_hh
#define _TextSelect_hh

#include <ROOT/REveElement.hxx>

namespace REX = ROOT::Experimental;
using namespace ROOT::Experimental;

namespace mu2e {
class TextSelect : public ROOT::Experimental::REveElement
{
  public:
    //TextSelect(int &test_) : test(test_){};
    /*TextSelect() = default; // ROOT needs a dictionary
    TextSelect(){};

    explicit TextSelect(
        int runn,
        int eventn)
        : REveElement{"TextSelect"}, runN(runn), eventN(eventn)
      {}*/
      
       void set(int run, int event);
       int get();
       int setRun(int run);
       //int test{0};
   private:
     int runN = 0;
     int eventN = 0;
     
     
     
};
}

#endif
