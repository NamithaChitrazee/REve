#ifndef DataProduct_HH_
#define DataProduct_HH_

#include <string>
namespace mu2e {

    /**
     * @brief A simple utility class to hold a string name for a data product.
     * * In the context of a large experiment like Mu2e, this class might be used
     * to consistently track the human-readable label (e.g., the Art InputTag)
     * of a collection being visualized.
     */
    class DataProduct {
    public:
        
        /**
         * @brief Constructor that initializes the name member.
         * @param name The string name to assign to this data product.
         */
        DataProduct(std::string name) : _name(name) {};
        
        /**
         * @brief Default constructor. Initializes name to an empty string.
         */
        DataProduct(){}; 
        
        /**
         * @brief Accessor method for the private name member.
         * Returns a reference, allowing the name to be modified (read/write access).
         * @return std::string& A reference to the internal name string.
         */
        std::string& name() { 
            return _name; 
        }

    private:
        // The private member storing the product's identifying name/label.
        std::string _name; 

    };
}  // namespace mu2e

#endif
