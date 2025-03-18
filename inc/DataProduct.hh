#ifndef REveMu2eDataProduct_HH_
#define REveMu2eDataProduct_HH_

#include <string>

namespace mu2e {

  class REveMu2eDataProduct {
   public:
    REveMu2eDataProduct(std::string name) : _name(name) {};
    REveMu2eDataProduct(){};
    std::string& name() { return _name; }

   private:
    std::string _name;

  };
}  // namespace mu2e

#endif
