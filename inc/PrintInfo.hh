#ifndef _PrintInfo_hh
#define _PrintInfo_hh

#include <ROOT/REveElement.hxx>
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryPoint.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

#include <iostream>
#include <vector>
#include <tuple>
#include <string>

namespace REX = ROOT::Experimental;
using namespace ROOT::Experimental;
namespace mu2e {

/**
 * @brief REve element responsible for holding event data pointers and printing diagnostic information.
 * * This element is added to the REve World and is the target for print commands issued from the browser.
 * * The member variables (tuples) are populated by the Mu2eEventDisplay::analyze() method.
 */
class PrintInfo : public ROOT::Experimental::REveElement
{
  public:
    PrintInfo() = default; // Required for the ROOT dictionary/serialization.

    // --- Data Storage for CaloClusterCollection ---
    const mu2e::CaloClusterCollection* clustercol = 0; // Single raw pointer (often for legacy use)
    std::vector<const mu2e::CaloClusterCollection*> calocluster_list; // List of collection pointers
    std::vector<std::string> calocluster_labels; // Labels corresponding to the collections
    // Primary storage: Tuple of (Labels, Collection Pointers) for all CaloCluster products in the event
    std::tuple<std::vector<std::string>, std::vector<const mu2e::CaloClusterCollection*>> fcalocluster_tuple;

    // --- Data Storage for MCTrajectoryCollection ---
    const mu2e::MCTrajectoryCollection *mctrajcol = 0; 
    std::vector<const mu2e::MCTrajectoryCollection*> mctrack_list; 
    std::vector<std::string> mctrack_labels;
    // Primary storage: Tuple for all Monte Carlo Trajectory collections
    std::tuple<std::vector<std::string>, std::vector<const mu2e::MCTrajectoryCollection*>> fmctrack_tuple;

    // --- Data Storage for SimParticleCollection (Monte Carlo) ---
    const mu2e::SimParticleCollection *simcol = 0; 
    std::vector<const mu2e::SimParticleCollection*> sim_list; 
    std::vector<std::string> sim_labels;
    // Primary storage: Tuple for all SimParticle collections
    std::tuple<std::vector<std::string>, std::vector<const mu2e::SimParticleCollection*>> fsim_tuple;

    // --- Data Storage for KalSeedPtrCollection (Reconstruction) ---
    const mu2e::KalSeedPtrCollection* kalSeedcol = 0; 
    std::vector<const mu2e::KalSeedPtrCollection*> track_list; 
    std::vector<std::string> track_labels;
    // Primary storage: Tuple for all Kalman Seed collections (reconstructed tracks)
    std::tuple<std::vector<std::string>, std::vector<const mu2e::KalSeedPtrCollection*>> ftrack_tuple;

    // --- Data Storage for CrvCoincidenceClusterCollection (Reconstruction) ---
    const mu2e::CrvCoincidenceClusterCollection* crvcoincol = 0; 
    std::vector<const mu2e::CrvCoincidenceClusterCollection*> crvcoin_list; 
    std::vector<std::string> crvcoin_labels;
    // Primary storage: Tuple for all CRV Coincidence Cluster collections
    std::tuple<std::vector<std::string>, std::vector<const mu2e::CrvCoincidenceClusterCollection*>> fcrvcoin_tuple;

    // --- Print Command Methods (Invoked by REve commands) ---

    /**
     * @brief Prints summary information for reconstructed data products 
     * (Tracks, CaloClusters, CRV Coincidence Clusters).
     */
    void PrintRecoInfo(); 

    /**
     * @brief Prints summary information for Monte Carlo data products 
     * (MCTrajectories and SimParticles).
     */
    void PrintMCInfo(); 

    // The following are more specialized print methods, likely called internally 
    // by PrintRecoInfo/PrintMCInfo or bound to separate commands.
    void PrintMCTrajInfo();
    void PrintSimPartInfo();
    void PrintKalInfo();
    void PrintCaloInfo();
    void PrintCRVInfo();
    };
}

#endif
