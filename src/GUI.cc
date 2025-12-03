#include "EventDisplay/inc/GUI.hh"

using namespace mu2e;


/**
 * @brief Writes core event metadata to the JSON structure for the REve browser interface.
 * * This method overrides the base REveElement method to ensure custom event information 
 * is included in the state update sent to the client (browser).
 * * @param j The nlohmann::json object being populated with element properties.
 * @param rnr_offset The offset used for rendering elements (typically 0).
 * @return int The status code returned by the base class method (0 on success).
 */
int GUI::WriteCoreJson(nlohmann::json &j, int rnr_offset)
{
    // Sets the path displayed in the browser's REve status bar or structure view.
    j["path"] = "Event/SubRun/Run"; 

    // Embed the current event identification numbers into the JSON object.
    // These members (feventid, fsubrunid, frunid) should be updated by the analyze() thread.
    j["eventid"] = feventid;
    j["subrunid"] = fsubrunid;
    j["runid"] = frunid;

    // This is a REve-specific instruction. It tells the browser client to execute 
    // the JavaScript function 'UT_refresh_event_info' after processing the stream 
    // update. This function is likely implemented in your custom HTML/JS (eventDisplay.html) 
    // to update displayed text fields with the new Run/Event IDs.
    j["UT_PostStream"] = "UT_refresh_event_info"; 

    // Call the base class method to ensure all standard REve properties are included.
    return ROOT::Experimental::REveElement::WriteCoreJson(j, 0);
}
