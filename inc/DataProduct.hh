#ifndef DataProduct_HH_
#define DataProduct_HH_

#include <string>

namespace mu2e {

  class DataProduct {
   public:
    DataProduct(std::string name) : _name(name) {};
    DataProduct(){};
    std::string& name() { return _name; }

   private:
    std::string _name;

  };
}  // namespace mu2e

#endif
