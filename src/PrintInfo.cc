#include "EventDisplay/inc/PrintInfo.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

using namespace mu2e;


/**
 * @brief Prints summary information for Monte Carlo (MC) event data.
 * * This function is typically invoked by the "PrintMCInfo" REve command button.
 * It currently delegates to the PrintSimInfo method.
 */
void PrintInfo::PrintMCInfo(){
    PrintMCTrajInfo();
    PrintSimPartInfo();
}

/**
 * @brief Prints summary information for all Reconstructed (Reco) event data.
 * * This function is typically invoked by the "PrintRecoInfo" REve command button.
 * It aggregates printing from tracking, calorimetry, and CRV components.
 */
void PrintInfo::PrintRecoInfo(){
    PrintKalInfo();
    PrintCaloInfo();
    PrintCRVInfo();
}

void PrintInfo::PrintMCTrajInfo(){

    std::vector<const MCTrajectoryCollection*> track_list = std::get<1>(fmctrack_tuple);
    
    if(track_list.size() > 0){
        
        // Set up stream for consistent formatting across all printouts
        std::cout << std::fixed << std::setprecision(3);
        
        // Loop over all available MCTrajectory collections (different Art module outputs).
        for(unsigned int j=0; j< track_list.size(); j++){
            const MCTrajectoryCollection* trajcol = track_list[j];
            
            if(trajcol !=0){
                
                std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
                
                std::cout << "\n======================================================================================================================================================" << std::endl;
                std::cout << "MC TRAJ (Primary) SIM PARTICLE INFORMATION: Collection " << j + 1 << std::endl;
                std::cout << "Number of SimParticles = " << trajcol->size() << std::endl;
                std::cout << "======================================================================================================================================================" << std::endl;

                // --- Table Header ---
                // Using setw() for fixed column width and right/left alignment
                std::cout << std::left << std::setw(6) << "ID"
                          << std::setw(10) << "PDGID"
                          << std::setw(9) << "ENERGY"
                          << std::setw(12) << "p0x" 
                          << std::setw(12) << "p0y" 
                          << std::setw(12) << "p0z" 
                          << std::setw(12) << "posx"
                          << std::setw(12) << "posy"
                          << std::setw(12) << "posz"
                          << std::setw(10) << "t0"
                          << std::setw(12) << "p1x" 
                          << std::setw(12) << "p1y" 
                          << std::setw(12) << "p1z" 
                          << std::setw(10) << "t1"
                          << std::endl;
                std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;


                // 3. Iterate over the SimParticle Ptrs and their associated MCTrajectories.
                for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++){
                    
                    // Access the SimParticle object
                    const mu2e::SimParticle* simPtr = trajectoryIter->first.get();
                    
                    // --- Initial (Start) Kinematics ---
                    int pdgID = simPtr->pdgId();
                    double energy = simPtr->startMomentum().e();
                    double t0 = simPtr->startGlobalTime();
                    
                    CLHEP::Hep3Vector p0 = simPtr->startMomentum().vect();
                    CLHEP::Hep3Vector pos0 = simPtr->startPosition();
                    
                    // --- Final (End) Kinematics ---
                    double t1 = simPtr->endGlobalTime();
                    CLHEP::Hep3Vector p1 = simPtr->endMomentum().vect();
                    
                    // Extract the unique SimParticle ID.
                    int id = simPtr->id().asInt();

                    // --- Print Data Row (using stream manipulators) ---
                    std::cout << std::left << std::setw(6) << id
                              << std::setw(10) << pdgID
                              << std::right << std::setw(9) << energy
                              << std::setw(12) << p0.x() 
                              << std::setw(12) << p0.y() 
                              << std::setw(12) << p0.z() 
                              << std::setw(12) << pos0.x()
                              << std::setw(12) << pos0.y()
                              << std::setw(12) << pos0.z()
                              << std::setw(10) << t0
                              << std::setw(12) << p1.x() 
                              << std::setw(12) << p1.y() 
                              << std::setw(12) << p1.z() 
                              << std::setw(10) << t1
                              << std::endl;
                }
                std::cout << "======================================================================================================================================================" << std::endl;
            }
        }
    }
    // Restore default stream formatting (optional but good practice)
    std::cout << std::scientific << std::setprecision(6);
}

void PrintInfo::PrintSimPartInfo(){
    std::vector<const SimParticleCollection*> sim_list = std::get<1>(fsim_tuple);
    std::vector<const MCTrajectoryCollection*> track_list = std::get<1>(fmctrack_tuple);
    if(track_list.size() > 0){
        
        // Set up stream for consistent formatting across all printouts
        std::cout << std::fixed << std::setprecision(3);
        
        // Loop over all available SimParticle collections.
        for(unsigned int j=0; j< sim_list.size(); j++){
            const SimParticleCollection* simcol = sim_list[j];
            
            // Safety check: ensure the pointer to the collection is valid.
            if(simcol !=0){
                
                std::cout << "\n======================================================================================================================================================" << std::endl;
                std::cout << "SIM PARTICLE INFORMATION (Collection " << j + 1 << ")" << std::endl;
                std::cout << "Number of SimParticles = " << simcol->size() << std::endl;
                std::cout << "======================================================================================================================================================" << std::endl;

                // --- Table Header ---
                std::cout << std::left << std::setw(6) << "ID"
                          << std::setw(10) << "PDGID"
                          << std::setw(9) << "ENERGY"
                          << std::setw(12) << "p0x" 
                          << std::setw(12) << "p0y" 
                          << std::setw(12) << "p0z" 
                          << std::setw(12) << "posx"
                          << std::setw(12) << "posy"
                          << std::setw(12) << "posz"
                          << std::setw(10) << "t0"
                          << std::setw(12) << "p1x" 
                          << std::setw(12) << "p1y" 
                          << std::setw(12) << "p1z" 
                          << std::setw(10) << "t1"
                          << std::endl;
                std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;


                // Iterate over the SimParticleCollection (map<SimParticle::key_type, SimParticle>).
                for( auto const& simpair : *simcol) {
                    const mu2e::SimParticle& simpart = simpair.second;
                    
                    // --- Extract Initial (Start) Kinematics ---
                    int pdgID = simpart.pdgId();
                    double energy = simpart.startMomentum().e();
                    double t0 = simpart.startGlobalTime();
                    
                    CLHEP::Hep3Vector p0 = simpart.startMomentum().vect();
                    CLHEP::Hep3Vector pos0 = simpart.startPosition();
                    
                    // --- Extract Final (End) Kinematics ---
                    double t1 = simpart.endGlobalTime();
                    CLHEP::Hep3Vector p1 = simpart.endMomentum().vect();
                    
                    // Extract the unique SimParticle ID.
                    int id = simpart.id().asInt();

                    // --- Print Data Row (using stream manipulators) ---
                    std::cout << std::left << std::setw(6) << id
                              << std::setw(10) << pdgID
                              << std::right << std::setw(9) << energy
                              << std::setw(12) << p0.x() 
                              << std::setw(12) << p0.y() 
                              << std::setw(12) << p0.z() 
                              << std::setw(12) << pos0.x()
                              << std::setw(12) << pos0.y()
                              << std::setw(12) << pos0.z()
                              << std::setw(10) << t0
                              << std::setw(12) << p1.x() 
                              << std::setw(12) << p1.y() 
                              << std::setw(12) << p1.z() 
                              << std::setw(10) << t1
                              << std::endl;
                }
                std::cout << "======================================================================================================================================================" << std::endl;
            }
        }
    }
    // Restore default stream formatting
    std::cout << std::scientific << std::setprecision(6);
}

void PrintInfo::PrintKalInfo(){
    // Access the particle data table for PDG IDs to names conversion
    auto const& ptable = GlobalConstantsHandle<ParticleDataList>(); 
    
    // Extract the vector of collection pointers and labels from the ftrack_tuple
    std::vector<const KalSeedPtrCollection*> ktrack_list = std::get<1>(ftrack_tuple);
    std::vector<std::string> names = std::get<0>(ftrack_tuple);
    
    if(ktrack_list.size() > 0){
        
        // Set up stream for consistent formatting (e.g., 3 decimal places)
        std::cout << std::fixed << std::setprecision(3);
        
        for(unsigned int j=0; j< ktrack_list.size(); j++){
            const KalSeedPtrCollection* seedcol = ktrack_list[j];
            
            if(seedcol->size() !=0){
                
                std::cout << "\n==============================================================================================================================" << std::endl;
                std::cout << "KALSEED INFORMATION | Collection: " << names[j] << " | Tracks: " << seedcol->size() << std::endl;
                std::cout << "==============================================================================================================================" << std::endl;

                // --- Table Header ---
                std::cout << std::left << std::setw(10) << "Particle"
                          << std::right << std::setw(10) << "Mom (MeV)"
                          << std::setw(10) << "Cos(Th)"
                          << std::setw(10) << "t0 (ns)"
                          << std::setw(10) << "d0 (mm)"
                          << std::setw(8) << "NActive"
                          << std::setw(12) << "Fit Cons"
                          << std::setw(10) << "E Calo"
                          << std::left << std::setw(15) << "Fit Status"
                          << std::endl;
                std::cout << "------------------------------------------------------------------------------------------------------------------------------" << std::endl;


                for(auto const& kseedptr : *seedcol) {
                    auto const& kseed = *kseedptr;
                    
                    // Taking the first segment is common for simple summary printing
                    auto const& kseg = kseed.segments()[0]; 
                    auto const& momvec = kseg.momentum3();
                    
                    double d0 = 0.0;
                    
                    // --- Extract d0 based on fit hypothesis ---
                    if(kseed.loopHelixFit()){
                        d0 = kseg.loopHelix().minAxisDist();
                    } else if(kseed.centralHelixFit()){
                        d0 = kseg.centralHelix().d0();
                    }else if(kseed.kinematicLineFit()){
                        d0 = kseg.kinematicLine().d0();
                    }
                    
                    // --- Count active hits ---
                    unsigned nactive = 0;
                    for (auto const& hit : kseed.hits()){
                        // Check if the hit state is active (>= inactive, assuming states are ordered)
                        if (hit.strawHitState() >= WireHitState::inactive) { 
                            ++nactive; 
                        }
                    }

                    // --- Extract Calo Energy ---
                    double caloE = kseed.hasCaloCluster() ? kseed.caloCluster()->energyDep() : 0.0;
                    
                    // --- Extract Particle Name ---
                    std::string particleName = ptable->particle(kseed.particle()).name();
                    
                    // --- Print Data Row ---
                    std::cout << std::left << std::setw(10) << particleName
                              << std::right << std::setw(10) << momvec.R()           // Momentum Magnitude
                              << std::setw(10) << cos(momvec.Theta())    // Cosine of Polar Angle
                              << std::setw(10) << kseed.t0Val()          // Track t0
                              << std::setw(10) << d0                     // d0 or minAxisDist
                              << std::setw(8) << nactive                // Number of active hits
                              << std::setw(12) << kseed.fitConsistency() // Fit Consistency
                              << std::setw(10) << caloE                  // Calo Cluster Energy
                              << std::left << std::setw(15) << kseed.status().stringRep() // Fit Status
                              << std::endl;
                }
                std::cout << "==============================================================================================================================" << std::endl;
            }
        }
    }
    // Restore default stream formatting
    std::cout << std::scientific << std::setprecision(6);
}

void PrintInfo::PrintCaloInfo(){
    // Extract the vector of collection pointers from the fcalocluster_tuple.
    std::vector<const mu2e::CaloClusterCollection*> calocluster_list = std::get<1>(fcalocluster_tuple);
    std::vector<std::string> names = std::get<0>(fcalocluster_tuple);

    if(!calocluster_list.empty()){
        
        // Set up stream for consistent formatting (e.g., 3 decimal places)
        std::cout << std::fixed << std::setprecision(3);
        
        // Loop over all CaloCluster collections
        for(unsigned int j = 0; j< calocluster_list.size(); j++){
            const mu2e::CaloClusterCollection* clustercol = calocluster_list[j];
            
            if(clustercol && !clustercol->empty()){
                
                std::cout << "\n=============================================================================" << std::endl;
                std::cout << "CALO CLUSTER INFORMATION | Collection: " << names[j] << " | Clusters: " << clustercol->size() << std::endl;
                std::cout << "=============================================================================" << std::endl;

                // --- Table Header ---
                std::cout << std::right << std::setw(10) << "Energy"
                          << std::setw(12) << "Time (ns)"
                          << std::setw(12) << "X COG (mm)"
                          << std::setw(12) << "Y COG (mm)"
                          << std::setw(12) << "Z COG (mm)" // Added Z coordinate for completeness
                          << std::endl;
                std::cout << "-----------------------------------------------------------------------------" << std::endl;

                // Iterate over the clusters in the current collection
                for(unsigned int i = 0; i < clustercol->size(); i++){
                    auto const& cluster= (*clustercol)[i];
                    
                    double energy = cluster.energyDep();
                    double time = cluster.time();
                    CLHEP::Hep3Vector cog = cluster.cog3Vector();
                    
                    // Print Data Row
                    std::cout << std::right << std::setw(10) << energy
                              << std::setw(12) << time
                              << std::setw(12) << cog.x()
                              << std::setw(12) << cog.y()
                              << std::setw(12) << cog.z()
                              << std::endl;
                }
                std::cout << "=============================================================================" << std::endl;
            }
        }
    }
    // Restore default stream formatting
    std::cout << std::scientific << std::setprecision(6);
}


void PrintInfo::PrintCRVInfo(){
    // Extract the vector of collection pointers and labels from the fcrvcoin_tuple.
    std::vector<const mu2e::CrvCoincidenceClusterCollection*> crvcoin_list = std::get<1>(fcrvcoin_tuple);
    std::vector<std::string> names = std::get<0>(fcrvcoin_tuple);

    if(!crvcoin_list.empty()){
        
        // Set up stream for consistent formatting (e.g., 3 decimal places)
        std::cout << std::fixed << std::setprecision(3);
        
        // Loop over all CrvCoincidenceCluster collections
        for(unsigned int j = 0; j< crvcoin_list.size(); j++){
            const mu2e::CrvCoincidenceClusterCollection* crvcoincol = crvcoin_list[j];
            
            if(crvcoincol && !crvcoincol->empty()){
                
                std::cout << "\n=============================================================================" << std::endl;
                std::cout << "CRV COINCIDENCE INFORMATION | Collection: " << names[j] << " | Clusters: " << crvcoincol->size() << std::endl;
                std::cout << "=============================================================================" << std::endl;

                // --- Table Header ---
                // Using right alignment for numerical data
                std::cout << std::right << std::setw(12) << "Avg Time (ns)"
                          << std::setw(12) << "Avg X (mm)"
                          << std::setw(12) << "Avg Y (mm)"
                          << std::setw(12) << "Avg Z (mm)"
                          << std::endl;
                std::cout << "-----------------------------------------------------------------------------" << std::endl;

                // Iterate over the coincidence clusters in the current collection
                for(unsigned int i = 0; i < crvcoincol->size(); i++){
                    mu2e::CrvCoincidenceCluster const &cluster = (*crvcoincol)[i];
                    
                    double avTime = cluster.GetAvgHitTime();
                    CLHEP::Hep3Vector avPos = cluster.GetAvgHitPos();

                    // Print Data Row
                    std::cout << std::right << std::setw(12) << avTime
                              << std::setw(12) << avPos.x()
                              << std::setw(12) << avPos.y()
                              << std::setw(12) << avPos.z()
                              << std::endl;
                }
                std::cout << "=============================================================================" << std::endl;
            }
        }
    }
    // Restore default stream formatting
    std::cout << std::scientific << std::setprecision(6);
}
