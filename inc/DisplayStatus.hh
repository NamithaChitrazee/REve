#include <string>
#include <mutex>

namespace mu2e {
    struct DisplayStatus {
        // Holds the command entered by the user in the REve interface
        std::string pendingCommand;
        
        // Mutex to protect access to the pendingCommand from different threads 
        // (REve UI thread vs. art event-processing thread).
        std::mutex commandMutex; 

        // Add other shared flags/data as needed, e.g.,
        bool shouldAdvanceEvent = false;
    };
}
