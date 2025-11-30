#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "EventDisplay/inc/MCInterface.hh"

using namespace mu2e;
namespace REX = ROOT::Experimental;

int MCInterface::Contains(std::vector<int> v, int x)
{
  return std::count(v.begin(), v.end(), abs(x));
}

/*------------Function to look up the particle name from its PDG Code----------*/
const char* MCInterface::GetParticleName(int PDGCode)
{
    // Initialize the pointer to the default (fallback) name.
    const char* pid = "other";
    
    // Check against common PDG codes and update the name.
    switch(PDGCode) {
        case PDGCode::e_minus:
            pid = "electron-";
            break;
        case PDGCode::e_plus:
            pid = "positron+";
            break;
        case PDGCode::mu_minus:
            pid = "muon-";
            break;
        case PDGCode::mu_plus:
            pid = "muon+";
            break;
        case PDGCode::pi_minus:
            pid = "pion-";
            break;
        case PDGCode::pi_plus:
            pid = "pion+";
            break;
        case PDGCode::pi0:
            pid = "pion0";
            break;
        case PDGCode::proton:
            pid = "proton";
            break;
        case PDGCode::n0:
            pid = "neutron";
            break;
        case PDGCode::gamma:
            pid = "gamma";
            break;
        // default is handled by the initial 'pid = "other";'
    }
    
    return pid;
}

/*------------Function to  set line color based on PDG ID ---------*/
void MCInterface::SetLineColorPID(int PDGCode, REX::REveLine *line)
{
    // Initialize color to the default/catch-all case ("other")
    Color_t color = kCyan - 4;

    // Assign color based on PDG Code
    switch(PDGCode) {
        case PDGCode::e_minus:
            color = kRed;
            break;
        case PDGCode::e_plus:
            color = kGreen;
            break;
        case PDGCode::mu_minus:
            color = kBlack;
            break;
        case PDGCode::mu_plus:
            color = kViolet;
            break;
        case PDGCode::pi_minus:
            color = kMagenta;
            break;
        case PDGCode::pi_plus:
            color = kRed - 7;
            break;
        case PDGCode::pi0:
            color = kGreen - 7;
            break;
        case PDGCode::proton:
            color = kBlue;
            break;
        case PDGCode::n0:
            color = kViolet - 2;
            break;
        case PDGCode::gamma:
            color = kOrange;
            break;
        // default case is handled by the initial assignment: color = kCyan - 4;
    }

    // Apply the determined color to the line
    line->SetLineColor(color);
}


void MCInterface::AddMCTrajectoryCollection(REX::REveManager *&eveMng, bool firstloop, 
                                         std::tuple<std::vector<std::string>, 
                                         std::vector<const MCTrajectoryCollection *>> mctrack_tuple, 
                                         REX::REveElement* &scene, 
                                         std::vector<int> particleIds, 
                                         bool extracted)
{
    std::cout << "[MCInterface::AddMCTrajectoryCollection()]" << std::endl;

    std::string drawfilename("EventDisplay/config/drawutils.txt");
    SimpleConfig drawconfig(drawfilename);

    const auto& track_list = std::get<1>(mctrack_tuple);

    mu2e::GeomHandle<mu2e::DetectorSystem> det;
    std::set<art::Ptr<mu2e::SimParticle>> drawn_particles;
    bool particle_type_found = false;

    // Loop over MCTrajectory Collections (J-loop)
    for(unsigned int j = 0; j < track_list.size(); j++){
        const mu2e::MCTrajectoryCollection* trajcol = track_list[j];

        if(trajcol == nullptr) continue;
        
        // Loop over Trajectories in the Collection (Map Iteration)
        // Iterate directly over the map elements (SimParticle Ptr -> MCTrajectory)
        for(const auto& pair : *trajcol)
        {
            const art::Ptr<mu2e::SimParticle>& sim_particle = pair.first;
            const mu2e::MCTrajectory& trajectory = pair.second;
            
            int pdg_id = sim_particle->pdgId();

            // Check if the SimParticle is in the user-defined list to plot
            if(Contains(particleIds, pdg_id) == 1){
                particle_type_found = true;
                
                // Duplicate Check
                if (drawn_particles.count(sim_particle)) {
                    continue; 
                }
                drawn_particles.insert(sim_particle);
                
                const std::vector<mu2e::MCTrajectoryPoint>& points = trajectory.points();
                if (points.empty()) continue;

                // Create Title String
                const char* particlename = GetParticleName(pdg_id);
                std::string energy = std::to_string(points[0].kineticEnergy());

                // Extract and convert Process Codes
                ProcessCode creationCode = sim_particle->creationCode();
                ProcessCode stoppingCode = sim_particle->stoppingCode();

                // Get the string names
                std::string creationName = creationCode.name();
                std::string stoppingName = stoppingCode.name();


                std::string mctitle = "MCTrajectory: " + std::string(particlename) + '\n'
                    + "Sim ID: " + std::to_string(sim_particle->id().asInt()) + '\n'
                    + "Parent ID: " + std::to_string(sim_particle->parentId().asInt()) + '\n'
                    + "Kinetic Energy: " + energy + " MeV" + '\n'
                    + "Creation Process: " + creationName + " (" + std::to_string(creationCode) + ")" + '\n'
                    + "Stopping Process: " + stoppingName + " (" + std::to_string(stoppingCode) + ")" + '\n'
                    + "End Global Time: " + std::to_string(sim_particle->endGlobalTime()) + " ns";
                // Create and Populate REveLine
                auto line = new REX::REveLine(mctitle.c_str(), mctitle.c_str(), points.size()); 
                
                // Add points to the line
                for(const auto& point : points){
                    // Positions are typically in Mu2e coordinates (mm)
                    CLHEP::Hep3Vector Pos(point.x(), point.y(), point.z());
                    // Convert to detector coordinates (if needed)
                    CLHEP::Hep3Vector HitPos = det->toDetector(Pos); 
                    
                    // Apply Z-range cut (CRV cut)
                    if(pointmmTocm(HitPos.y()) < 800.0){ 
                        line->SetNextPoint(pointmmTocm(HitPos.x()), 
                                           pointmmTocm(HitPos.y()), 
                                           pointmmTocm(HitPos.z()));
                    }
                }
                
                // Styling and Scene Addition
                if (line->GetSize() > 0) {
                    SetLineColorPID(pdg_id, line);
                    line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
                    scene->AddElement(line);
                }
            } 
        } // End of map iteration
    } // End of collection loop
    
    if (!particle_type_found && !particleIds.empty()) {
        std::cout << "Warning: No Particles of User-Specified Type In File." << std::endl;
    }
}

void MCInterface::AddSurfaceStepCollection(REX::REveManager *&eveMng, bool firstloop,
                                         std::tuple<std::vector<std::string>, 
                                         std::vector<const SurfaceStepCollection *>> surfstep_tuple, 
                                         REX::REveElement* &scene, 
                                         std::vector<int> particleIds, 
                                         bool extracted)
{
    std::cout << "[MCInterface::AddSurfaceStepCollection() ]" << std::endl;

    // --- Scene Setup and Cleanup ---
    // If this is not the first loop (i.e., not the first event being displayed), clear the previous elements.
    if (!firstloop) {
        scene->DestroyElements();
    }
    
    // Extract the list of SurfaceStep collections and their associated names/tags from the tuple.
    const std::vector<const mu2e::SurfaceStepCollection*>& ssteps_list = std::get<1>(surfstep_tuple);
    const std::vector<std::string>& names = std::get<0>(surfstep_tuple);
    
    // Load drawing configuration (assuming this is used for marker size/style).
    std::string drawfilename("EventDisplay/config/drawutils.txt");
    mu2e::SimpleConfig drawconfig(drawfilename);

    // --- Main Loop: Iterate over all SurfaceStep collections ---
    for(unsigned int j = 0; j < ssteps_list.size(); j++){
        const mu2e::SurfaceStepCollection* surfstepcol = ssteps_list[j];
        
        if(surfstepcol == nullptr) continue; // Skip if collection pointer is null

        // Make a compound object to store all SurfaceStep hits from this collection, 
        // ensuring they are grouped logically in the REve display tree.
        std::string comptitle = "SurfaceStepCollection: " + names[j];
        auto allpoints = new REX::REveCompound(comptitle.c_str(), comptitle.c_str(), 1);

        // Inner Loop: Iterate over individual SurfaceStep objects
        for(const auto& surfstep : *surfstepcol) {
            
            // Get the PDG ID of the SimParticle associated with this step.
            auto pdgid = surfstep.simParticle()->pdgId();
            
            // Check if the SimParticle's PDG ID is in the user-defined list to plot.
            int plot_this_particle = Contains(particleIds, pdgid);
            
            // Get the midpoint position of the step (in Mu2e coordinates, mm).
            auto midpos = surfstep.midPosition();
            
            if(plot_this_particle == 1){
                // Make a descriptive label for the mouseover tooltip.
                std::string momentum = std::to_string(surfstep.momentum().R());
                std::string edep = std::to_string(surfstep.energyDeposit());
                
                std::string mctitle = "SurfaceStep on " + surfstep.surfaceId().name() + '\n'
                    + " x " + std::to_string(midpos.X()) + '\n'
                    + " y " + std::to_string(midpos.Y()) + '\n'
                    + " z " + std::to_string(midpos.Z()) + '\n'
                    + " time: " + std::to_string(surfstep.time()) + '\n'
                    + " momentum " + momentum + " MeV/c, energy loss = " + edep + " MeV";
                
                // Create the REve marker point for this step.
                auto surfpoint = new REX::REvePointSet(mctitle.c_str(), mctitle.c_str(), 1);
                
                // Configure marker appearance.
                surfpoint->SetMarkerStyle(MCInterface::mstyle);
                surfpoint->SetMarkerSize(MCInterface::msize);
                surfpoint->SetMarkerColor(kBlack);
                
                // Set the point position, converting coordinates from mm to cm.
                surfpoint->SetNextPoint(pointmmTocm(midpos.X()), 
                                        pointmmTocm(midpos.Y()), 
                                        pointmmTocm(midpos.Z()));
                
                // Add the individual point to the compound group.
                allpoints->AddElement(surfpoint);
            }
        } // End of inner loop (individual SurfaceSteps)

        // Only add the compound group to the scene if it contains any plotted points.
        scene->AddElement(allpoints);
    } // End of main loop (SurfaceStep collections)
}

void MCInterface::AddSimParticleCollection(REX::REveManager *&eveMng, bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const SimParticleCollection *>> sim_tuple, REX::REveElement* &scene, std::vector<int> particleIds, bool extracted){
  std::cout<<"[MCInterface::AddSimParticleCollection() ]"<<std::endl;
  std::vector<const SimParticleCollection*> sim_list = std::get<1>(sim_tuple);
  std::vector<std::string> names = std::get<0>(sim_tuple);

  if(sim_list.size() !=0){
    for(unsigned int i=0; i < sim_list.size(); i++){
      std::string comptitle = "SimParticleCollection" + names[i];

      // make compund object to store hits
      std::string drawfilename("EventDisplay/config/drawutils.txt");
      SimpleConfig drawconfig(drawfilename);

      // eXtract the track and input tag:
      std::vector<const SimParticleCollection*> sim_list = std::get<1>(sim_tuple);
      std::vector<std::string> names = std::get<0>(sim_tuple);

      // Loop over SimParticle
      //for(unsigned int j=0; j< sim_list.size(); j++){
        const SimParticleCollection* simcol = sim_list[i];
      
        if(simcol!=0){
          auto SimCollection = new REX::REveCompound("SimParticles","SimParticles",1);
          for( auto const& simpair : *simcol) {
            // Check user defined list of particles to plot
            const mu2e::SimParticle& simpart = simpair.second;
            auto pdgid = simpart.pdgId();
            auto startCode = simpart.creationCode().name() ;
            auto stopCode = simpart.stoppingCode().name()  ;
            int x = Contains(particleIds,pdgid);
            GeomHandle<DetectorSystem> det;
            

            if(x == 1){
              // Make label
              //std::string momentum = 0;//std::to_string(simpart.startMomentum().R());
              std::string edep = std::to_string(simpart.endKineticEnergy());
              CLHEP::Hep3Vector StartPos = det->toDetector(simpart.startPosition());
              CLHEP::Hep3Vector EndPos = det->toDetector(simpart.endPosition());
              double momentum = sqrt(simpart.startMomentum().x()*simpart.startMomentum().x()+simpart.startMomentum().y()*simpart.startMomentum().y() + simpart.startMomentum().z()*simpart.startMomentum().z());
              double endmomentum = sqrt(simpart.endMomentum().x()*simpart.endMomentum().x()+simpart.endMomentum().y()*simpart.endMomentum().y() + simpart.endMomentum().z()*simpart.endMomentum().z());
              std::string mctitle_start = " SimParticle PDGid " + std::to_string(simpart.pdgId()) + '\n'
                + " Creation code " + (startCode) + " Stopping code " + (stopCode) + '\n'
                + " Start Position: " + '\n'
                + " x "  + std::to_string(StartPos.x())
                + " y " + std::to_string(StartPos.y())
                + " z " + std::to_string(StartPos.z())
                + " time :" + std::to_string(simpart.startGlobalTime()  )+  '\n'
                + " Start Momentum " + std::to_string(momentum) + " Start Energy " + std::to_string(simpart.startMomentum().e()) + '\n'
                + " End Position: " + '\n'
                + " x "  + std::to_string( EndPos.x())
                + " y " + std::to_string( EndPos.y())
                + " z " + std::to_string( EndPos.z())
                + " time :" + std::to_string(simpart.endGlobalTime()  )+ '\n'
                + " End Momentum " + std::to_string(endmomentum) + " End Energy " + std::to_string(simpart.endMomentum().e()) ;

              
              //auto simpoint_start = new REX::REvePointSet(mctitle_start,mctitle_start,1);
              // create line with the above label
              auto simpart_line = new REX::REveLine(mctitle_start,mctitle_start,1);
              simpart_line->SetNextPoint(pointmmTocm(StartPos.x()), pointmmTocm(StartPos.y()),pointmmTocm(StartPos.z()));
              simpart_line->SetNextPoint(pointmmTocm(EndPos.x()), pointmmTocm(EndPos.y()),pointmmTocm(EndPos.z()));
              // set line colour
              SetLineColorPID(pdgid, simpart_line );
              simpart_line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
              SimCollection->AddElement(simpart_line);

              
            }
         // }
        }
         scene->AddElement(SimCollection);
      }
    }
  }
}
