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
        
        // Loop over Trajectories in the Collection (Map Iteration) ---
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

                // --- 5. Create Title String ---
                const char* particlename = GetParticleName(pdg_id);
                std::string energy = std::to_string(points[0].kineticEnergy());

                // Extract and convert Process Codes
                ProcessCode creationCode = sim_particle->creationCode();
                ProcessCode stoppingCode = sim_particle->stoppingCode();

                // Get the string names (assuming GetProcessCodeName utility exists)
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
                
                // --- 7. Styling and Scene Addition ---
                if (line->GetSize() > 0) {
                    SetLineColorPID(pdg_id, line);
                    line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
                    scene->AddElement(line);
                }
            } 
        } // End of map iteration
    } // End of collection loop
    
    // --- 8. Final Warning ---
    if (!particle_type_found && !particleIds.empty()) {
        std::cout << "Warning: No Particles of User-Specified Type In File." << std::endl;
    }
}

void MCInterface::AddSurfaceStepCollection(REX::REveManager *&eveMng, bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const SurfaceStepCollection *>> surfstep_tuple, REX::REveElement* &scene, std::vector<int> particleIds, bool extracted){
  std::cout<<"[MCInterface::AddSurfaceStepCollection() ]"<<std::endl;
  std::vector<const SurfaceStepCollection*> ssteps_list = std::get<1>(surfstep_tuple);
  std::vector<std::string> names = std::get<0>(surfstep_tuple);
  
  if(ssteps_list.size() !=0){
    for(unsigned int i=0; i < ssteps_list.size(); i++){
      std::string comptitle = "SurfaceStepCollection" + names[i];
      
      // make compund object to store hits
      auto allpoints = new REX::REveCompound(comptitle,comptitle,1);
      std::string drawfilename("EventDisplay/config/drawutils.txt");
      SimpleConfig drawconfig(drawfilename);

      // eXtract the track and input tag:
      std::vector<const SurfaceStepCollection*> surfstep_list = std::get<1>(surfstep_tuple);
      std::vector<std::string> names = std::get<0>(surfstep_tuple);

      // Loop over surface steps
      for(unsigned int j=0; j< surfstep_list.size(); j++){
        const SurfaceStepCollection* surfstepcol = surfstep_list[j];
        if(surfstepcol!=0){
          for( auto const& surfstep : *surfstepcol) {
            // Check user defined list of particles to plot
            auto pdgid = surfstep.simParticle()->pdgId();
            int x = Contains(particleIds,pdgid);
            auto midpos = surfstep.midPosition();
            if(x == 1){
              // Make label
              std::string momentum = std::to_string(surfstep.momentum().R());
              std::string edep = std::to_string(surfstep.energyDeposit());
              std::string mctitle = " SurfaceStep on " +  surfstep.surfaceId().name() + '\n'
                + " x "  + std::to_string(midpos.X())
                + " y " + std::to_string(midpos.Y())
                + " z " + std::to_string(midpos.Z())
                + " time :" + std::to_string(surfstep.time())+  '\n'
                + " momentum " + momentum + " MeV/c, energy loss = " + edep + "MeV";
              // add point
              auto surfpoint = new REX::REvePointSet(mctitle,mctitle,1);
              surfpoint->SetMarkerStyle(MCInterface::mstyle);
              surfpoint->SetMarkerSize(MCInterface::msize);
              surfpoint->SetMarkerColor(kBlack);
              surfpoint->SetNextPoint(pointmmTocm(midpos.X()),pointmmTocm(midpos.Y()) ,pointmmTocm(midpos.Z()));
              //scene->AddElement(surfpoint);
              allpoints->AddElement(surfpoint);
            }
          }
        }
      }
      scene->AddElement(allpoints);
    }
  }
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
