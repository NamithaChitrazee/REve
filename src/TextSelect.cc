#include "Mu2eEventDisplay/inc/TextSelect.hh"

using namespace mu2e;

/*TextSelect::TextSelect(){};

TextSelect::TextSelect(
    int runn,
    int eventn)
    : REveElement{"TextSelect"}, runN(runn), eventN(eventn)
  {};*/
  
  
int TextSelect::get(){
  return runN;
}

int TextSelect::setRun(int run){
        int N = run;
        std::cout<<"set run "<<N<<std::endl;
        return N;
      }

void TextSelect::set(int run, int event){
  int N = setRun(run);
  runN = N;
  
  std::cout<<" [TextSelect::set] set to "<<N<<std::endl;
  std::cout<<" [TextSelect::set] run number in this class "<<runN<<std::endl;
}
   

  
