#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "EventDisplay/inc/MCInterface.hh"
using namespace mu2e;
namespace REX = ROOT::Experimental;

int MCInterface::Contains(std::vector<int> v, int x)
{
  return std::count(v.begin(), v.end(), abs(x));
}

const char* MCInterface::GetParticleName(int PDGCode){
  const char* pid = "pid";
  switch(PDGCode) {
    case PDGCode::e_minus:
      pid = "electron -";
      break;
    case PDGCode::e_plus:
      pid = "positron +";
      break;
    case PDGCode::mu_minus:
      pid = "muon - ";
      break;
    case PDGCode::mu_plus:
      pid = "muon + ";
      break;
    case PDGCode::pi_minus:
      pid = "pion -";
      break;
    case PDGCode::pi_plus:
      pid = "pion +";
      break;
    case PDGCode::pi0:
      pid = "pion 0";
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
    default:
      pid = "other";
      break;
  }
  return pid;
}

void MCInterface::SetLineColorPID(int PDGCode,REX::REveLine *line){
  Color_t color;
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
      color = kRed-7;
      break;
    case PDGCode::pi0:
      color = kGreen-7;
      break;
    case PDGCode::proton:
      color = kBlue;
      break;
    case PDGCode::n0:
      color = kViolet-2;
      break;
    case PDGCode::gamma:
      color = kOrange;
      break;
    default:
      color = kCyan-4;
      break;
  }
  line->SetLineColor(color);
}


void MCInterface::AddMCTrajectoryCollection(REX::REveManager *&eveMng, bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const MCTrajectoryCollection *>> mctrack_tuple, REX::REveElement* &scene, std::vector<int> particleIds, bool extracted){
  std::cout<<"[ REveMCInterface::AddMCTrajectoryCollection() ]"<<std::endl;
  std::string drawfilename("EventDisplay/config/drawutils.txt");
  SimpleConfig drawconfig(drawfilename);

  // eEtract the track and input tag:
  std::vector<const MCTrajectoryCollection*> track_list = std::get<1>(mctrack_tuple);
  std::vector<std::string> names = std::get<0>(mctrack_tuple);

  // Loop over tracks:
  for(unsigned int j=0; j< track_list.size(); j++){
    // extract trajectory
    const MCTrajectoryCollection* trajcol = track_list[j];

    if(trajcol!=0){
      std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
      bool isSame = false;
      //create vector of art Ptr to particles:
      std::vector<art::Ptr<SimParticle> > allParts;
      for(unsigned int k = 0; k < trajcol->size(); k++){
        for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++)
        {
          // Check user defined list of particles to plot
          int x = Contains(particleIds,trajectoryIter->first->pdgId());
          // get particle name
          const char* particlename = GetParticleName(trajectoryIter->first->pdgId());
          if(x == 1){
            const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();

            std::string energy = std::to_string(points[0].kineticEnergy());
            // check for duplicates
            for(unsigned int ipart = 0; ipart < allParts.size() ; ipart++){
              MCRelationship checkrel(trajectoryIter->first,allParts.at(ipart));
              if(checkrel==MCRelationship::same) {
                  isSame = true;
                }
            }
            if(!isSame) { allParts.push_back(trajectoryIter->first); }

            std::string mctitle = " MCTrajectory tag : particle = " + std::string(particlename)  +  '\n'
              + " energy = " + energy + "MeV" +  '\n'
              + " creation code = " + std::to_string(trajectoryIter->first->creationCode()) +  '\n'
              + " stopping code = " + std::to_string(trajectoryIter->first->stoppingCode())  +  '\n'
              + " end global time = " + std::to_string(trajectoryIter->first->endGlobalTime())
              + " ns";
            // create line with the above label
            auto line = new REX::REveLine(mctitle,mctitle,1);
            // add points
            for(unsigned int i=0; i < points.size();i++){
              CLHEP::Hep3Vector Pos(points[i].x(), points[i].y(), points[i].z());
              GeomHandle<DetectorSystem> det;
              CLHEP::Hep3Vector HitPos = det->toDetector(Pos);
              if(pointmmTocm(HitPos.y()) < 800){ // a reasonable height above the CRV
                line->SetNextPoint(pointmmTocm((HitPos.x())),pointmmTocm((HitPos.y())),pointmmTocm(HitPos.z()));
              }
            }
            // set line colour
            SetLineColorPID(trajectoryIter->first->pdgId(),line);
            line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
            if(!isSame) scene->AddElement(line);
          } else std::cout<<"Warning: No Particles of User-Specified Type In File "<<std::endl;
        }
      }
    }
  }
}

void MCInterface::AddSurfaceStepCollection(REX::REveManager *&eveMng, bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const SurfaceStepCollection *>> surfstep_tuple, REX::REveElement* &scene, std::vector<int> particleIds, bool extracted){
  std::cout<<"[ REveMCInterface::AddSurfaceStepCollection() ]"<<std::endl;
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
  std::cout<<"[ REveMCInterface::AddSimParticleCollection() ]"<<std::endl;
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
