#ifndef _GUI_hh
#define _GUI_hh

#include <ROOT/REveElement.hxx>
#include "nlohmann/json.hpp"
namespace REX = ROOT::Experimental;
using namespace ROOT::Experimental;

namespace mu2e {
/**
 * @brief Represents the custom Graphical User Interface component within the REve visualization.
 * * This class holds the event metadata and overrides methods to push that data to the web client.
 * Inherits from REveElement so it can be added to the REve scene and participate in state updates.
 */
class GUI : public ROOT::Experimental::REveElement
{
  public:
     // --- Event Display Variables (Set by the Art thread) ---
     
     // Current Event ID being displayed. Used when writing to JSON.
     int feventid{0};
     // Current Subrun ID being displayed.
     int fsubrunid{0};
     // Current Run ID being displayed.
     int frunid{0};

     // --- REve Overrides and Communication ---

     /**
      * @brief Overrides the base method to serialize custom class data into a JSON stream.
      * * This is the crucial function that sends feventid/frunid to the web browser client.
      * @param j The nlohmann::json object to populate.
      * @param rnr_offset The rendering offset (typically 0).
      * @return int Status code from the base class.
      */
     int WriteCoreJson(nlohmann::json &j, int rnr_offset) override;

     // Storage for the Run number requested by the user.
     int runNumber{0};
     // Storage for the Event number requested by the user.
     int eventNumber{0};
};
}
#endif
