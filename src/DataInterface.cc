#include "EventDisplay/inc/DataInterface.hh"
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/RecoDataProducts/inc/TrkStrawHitSeed.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include <sstream>
#include <iomanip>

using namespace mu2e;
namespace REX = ROOT::Experimental;

std::string drawfilename("Mu2eEventDisplay/config/drawutils.txt");
SimpleConfig drawconfig(drawfilename);


void DataInterface::AddCaloDigis(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CaloDigiCollection*>> calodigi_tuple, REX::REveElement* &scene, bool addCrystalDraw){

  std::cout<<"[DataInterface] AddCaloDigi "<<std::endl;
  std::vector<const CaloDigiCollection*> calodigi_list = std::get<1>(calodigi_tuple);
  std::vector<std::string> names = std::get<0>(calodigi_tuple);
  for(unsigned int j = 0; j< calodigi_list.size(); j++){

    const CaloDigiCollection* calodigicol = calodigi_list[j];
    if(calodigicol->size() != 0){
      if(!firstLoop_){
        scene->DestroyElements();;
      }
      mu2e::Calorimeter const &cal = *(mu2e::GeomHandle<mu2e::Calorimeter>());
      for(unsigned int i = 0; i < calodigicol->size(); i++){
        mu2e::CaloDigi const  &digi= (*calodigicol)[i];
        int sipmID = digi.SiPMID();
        //bool isEven = sipmID % 2;
        int cryID = sipmID/2;
        //std::cout<<"crystal ID "<<cryID<<" "<<isEven<<" "<<digi.t0()<<std::endl;
        Crystal const &crystal = cal.crystal(cryID);
        double crystalXLen = pointmmTocm(crystal.size().x());
        double crystalYLen = pointmmTocm(crystal.size().y());
        double crystalZLen = pointmmTocm(crystal.size().z());


        GeomHandle<DetectorSystem> det;
        double zpos = 0;
        double diskID = 0;
        if(cryID < 674) {
          zpos = 186.53;
          diskID = 0;
        }
        if(cryID >= 674) {
          zpos = 256.53;
          diskID = 1;
        }
        CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDisk(diskID,crystal.position()) ;

        std::string crytitle =  "Crystal ID = " + std::to_string(cryID) +  '\n'
          + " Time = " + std::to_string(digi.t0())+" ns ";
        char const *crytitle_c = crytitle.c_str();

        // plot the crystals which are present in this event:
        auto b = new REX::REveBox("crystal",crytitle_c);
        b->SetMainColor(602);
        double width = crystalXLen/2;
        double height = crystalYLen/2;
        double thickness = crystalZLen/2;
        b->SetVertex(0, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height ,pointmmTocm(crystalPos.z())- thickness + 2*thickness + zpos);//---
        b->SetVertex(1, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height, pointmmTocm(crystalPos.z())- thickness + 2*thickness  +zpos);//-+-
        b->SetVertex(2, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height ,pointmmTocm(crystalPos.z())- thickness + 2*thickness  + zpos);//++-
        b->SetVertex(3, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, pointmmTocm(crystalPos.z())-thickness  + 2*thickness + zpos);//+--
        b->SetVertex(4, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height ,pointmmTocm(crystalPos.z())+ thickness + 2*thickness + zpos);//--+
        b->SetVertex(5, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height , pointmmTocm(crystalPos.z())+ thickness  + 2*thickness + zpos);//-++
        b->SetVertex(6, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height , pointmmTocm(crystalPos.z()) + thickness + 2*thickness+zpos); //+++
        b->SetVertex(7,pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, pointmmTocm(crystalPos.z())+ thickness + 2*thickness + zpos);//+-+
        scene->AddElement(b);
      }
    }
  }
}
void DataInterface::AddCaloClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CaloClusterCollection*>> calocluster_tuple, REX::REveElement* &scene, bool addCrystalDraw, double t1, double t2){

  std::cout<<"[DataInterface] AddCaloClusters "<<std::endl;
  std::vector<const CaloClusterCollection*> calocluster_list = std::get<1>(calocluster_tuple);
  std::vector<std::string> names = std::get<0>(calocluster_tuple);
  for(unsigned int j = 0; j< calocluster_list.size(); j++){

    // Extract cluster list
    const CaloClusterCollection* clustercol = calocluster_list[j];
    if(clustercol->size() != 0){
      if(!firstLoop_){
        scene->DestroyElements();;
      }
      Color_t color = kWhite;
      mu2e::Calorimeter const &cal = *(mu2e::GeomHandle<mu2e::Calorimeter>());
      GeomHandle<DetectorSystem> det;
      // Add crystals
      double maxE = 1e-6;
      double minE = 1000;
      using CCCPtr = const CaloCluster*;
      std::list<CCCPtr> cluList;
      //std::array<Color_t, 10> colorList {kViolet, kBlue, kBlue, kGreen, kGreen, kYellow, kYellow, kOrange, kOrange, kRed};
      //std::vector<Color_t> colors;

      for(unsigned int i = 0; i < clustercol->size(); i++){
        auto const& cluster= (*clustercol)[i];
        cluList.push_back(&cluster);
        if(cluster.energyDep() > maxE) maxE = cluster.energyDep();
        if(cluster.energyDep() < minE) minE = cluster.energyDep();

        // sort in energy
        cluList.sort([] (const CaloCluster* lhs, const CaloCluster* rhs) {return lhs->energyDep() > rhs->energyDep();} );
      }

      /*if(cluList.size() <= 2){
        Color_t color;
        for(auto const& cluster : cluList){
        if(cluster.energyDep() == minE){
        color = kViolet;
        colors.push_back(color);
        }
        if(cluster.energyDep() == maxE){
        color = kRed;
        colors.push_back(color);
        }
        }
        }
        unsigned int c = -1;
        if(cluList.size() > 2){
        for(auto const& cluster : cluList){
        c+=1;
        Color_t color;
        if(cluster.energyDep() == minE){
        color = kViolet;
        colors.push_back(color);
        }
        if(cluster.energyDep() == maxE){
        color = kRed;
        colors.push_back(color);
        }
        if (c > 1 and c < colorList.size()-1){
        color = colorList.at(c);
        colors.push_back(color);
        }
        if (c > 1 and c >= colorList.size()-1){
        color = kRed;
        colors.push_back(color);
        }
        }
        } */

      unsigned int i = 0;
      for(auto cptr : cluList){
        auto const& cluster = *cptr;
        //mu2e::CaloCluster const  &cluster= cluList.at(i);
        i += 1;

        // Info for label:
        std::string cluster_energy = std::to_string(cluster.energyDep());
        std::string cluster_time = std::to_string(cluster.time());
        std::string cluster_x = std::to_string(cluster.cog3Vector().x());
        std::string cluster_y = std::to_string(cluster.cog3Vector().y());
        // Extract center of gravity, convert to coord sys
        CLHEP::Hep3Vector COG(cluster.cog3Vector().x(),cluster.cog3Vector().y(), cluster.cog3Vector().z());
        CLHEP::Hep3Vector crystalPos   = cal.geomUtil().mu2eToDisk(cluster.diskID(),COG);
        CLHEP::Hep3Vector pointInMu2e = det->toMu2e(crystalPos);
        // Info for label
        std::string cluster_z = std::to_string(abs(pointInMu2e.z()));
        // Make label and REve objects
        std::string label =   " Instance = " + names[0] +  '\n'
          + " Energy Dep. = "+cluster_energy+" MeV "+   '\n'
          + " Time = "+cluster_time+" ns " +  '\n'
          + " Pos =  ("+cluster_x+","+cluster_y+","+cluster_z+") mm";
        std::string name = "disk" + std::to_string(cluster.diskID()) + label;
        auto ps1 = new REX::REvePointSet(name, "CaloClusters Disk 1: "+label,0);
        auto ps2 = new REX::REvePointSet(name, "CaloClusters Disk 2: "+label,0);


        // Set positions of clusters
        if(cluster.diskID() == 0) ps1->SetNextPoint(pointmmTocm(COG.x()), pointmmTocm(COG.y()) , abs(pointmmTocm(pointInMu2e.z())));
        if(cluster.diskID() == 1) ps2->SetNextPoint(pointmmTocm(COG.x()), pointmmTocm(COG.y()) , abs(pointmmTocm(pointInMu2e.z())));

        // Set draw options
        /*TColor color1;
          color1.SetPalette(1,0);
          Color_t color = color1.GetNumber();// colors[i];*/

        if(abs(cluster.time() - t2) < 20) color = kRed;
        if(abs(cluster.time() - t2) >=  20 and abs(cluster.time() - t2) < 100) color = kOrange;
        if(abs(cluster.time() - t2) >= 100 and abs(cluster.time() - t2) < 200) color = kYellow;
        if(abs(cluster.time() - t2) >= 200 and abs(cluster.time() - t2) < 400) color = kGreen;
        if(abs(cluster.time() - t2) >= 400 and abs(cluster.time() - t2) < 600) color = kBlue;
        if(abs(cluster.time() - t2) >= 600 ) color = kViolet;

        ps1->SetMarkerColor(color);
        ps1->SetMarkerStyle(DataInterface::mstyle);
        ps1->SetMarkerSize(DataInterface::msize);

        ps2->SetMarkerColor(color);
        ps2->SetMarkerStyle(DataInterface::mstyle);
        ps2->SetMarkerSize(DataInterface::msize);

        // Add to REve world
        scene->AddElement(ps1);
        scene->AddElement(ps2);

        // Add crystals
        if(addCrystalDraw){
          auto allcryhits = new REX::REveCompound("CrystalHits for Cluster "+std::to_string(i),"CrystalHits for Cluster "+std::to_string(i),1);
          for(unsigned h =0 ; h < cluster.caloHitsPtrVector().size();h++)     {
            art::Ptr<CaloHit>  crystalhit = cluster.caloHitsPtrVector()[h];
            int cryID = crystalhit->crystalID();

            //int diskID = cal.crystal(crystalhit->crystalID()).diskID();
            Crystal const &crystal = cal.crystal(cryID);
            double crystalXLen = pointmmTocm(crystal.size().x());
            double crystalYLen = pointmmTocm(crystal.size().y());
            double crystalZLen = pointmmTocm(crystal.size().z());

            GeomHandle<DetectorSystem> det;
            CLHEP::Hep3Vector crystalPos = cal.geomUtil().mu2eToDisk(cluster.diskID(),crystal.position()) ;

            //  make title
            std::string crytitle =   "disk"+std::to_string(cal.crystal(crystalhit->crystalID()).diskID()) + " Crystal Hit = " + std::to_string(cryID) +  '\n'
              + " Energy Dep. = "+std::to_string(crystalhit->energyDep())+" MeV "+   '\n'
              + " Time = "+std::to_string(crystalhit->time())+" ns ";
            char const *crytitle_c = crytitle.c_str();
            auto b = new REX::REveBox(crytitle_c,crytitle_c);

            // plot the crystals which are present in this event in lego:
            b->SetMainColor(color);//s[i]);


            double width = crystalXLen/2;
            double height = crystalYLen/2;
            double thickness = crystalhit->energyDep()/maxE * crystalZLen/2; //length proportional to energy

            b->SetVertex(0, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height ,pointmmTocm(crystalPos.z())- thickness  + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//---
            b->SetVertex(1, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height, pointmmTocm(crystalPos.z())- thickness   +abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//-+-
            b->SetVertex(2, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height ,pointmmTocm(crystalPos.z())- thickness   + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//++-
            b->SetVertex(3, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, pointmmTocm(crystalPos.z())-thickness   + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//+--
            b->SetVertex(4, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())- height ,pointmmTocm(crystalPos.z())+ thickness  + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//--+
            b->SetVertex(5, pointmmTocm(crystalPos.x()) - width, pointmmTocm(crystalPos.y())+ height , pointmmTocm(crystalPos.z())+ thickness   + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//-++
            b->SetVertex(6, pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())+ height , pointmmTocm(crystalPos.z()) + thickness +abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2); //+++
            b->SetVertex(7,pointmmTocm(crystalPos.x()) + width, pointmmTocm(crystalPos.y())- height, pointmmTocm(crystalPos.z())+ thickness + abs(pointmmTocm(pointInMu2e.z()))+crystalZLen/2);//+-+
            allcryhits->AddElement(b);
          }
          scene->AddElement(allcryhits);
        }

      }
    }
  }
}

void DataInterface::AddComboHits(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const ComboHitCollection*>> combohit_tuple, REX::REveElement* &scene, bool strawdisplay, bool AddErrorBar_){

  std::vector<const ComboHitCollection*> combohit_list = std::get<1>(combohit_tuple);
  std::vector<std::string> names = std::get<0>(combohit_tuple);

  // Loop over hit lists
  for(unsigned int j = 0; j < combohit_list.size(); j++){
    const ComboHitCollection* chcol = combohit_list[j];
    int colour = (j+3);
    if(chcol->size() !=0 ){
      // Loop over hits
      for(unsigned int i=0; i< chcol->size(); i++){
        mu2e::ComboHit const  &hit= (*chcol)[i];
        // Display hit straws if selected too in FCL file
        if(strawdisplay){
          mu2e::GeomHandle<mu2e::Tracker> tracker;
          const auto& allStraws = tracker->getStraws();
          int sid = hit._sid.asUint16();
          CLHEP::Hep3Vector sposi(0.0,0.0,0.0), sposf(0.0,0.0,0.0);
          const mu2e::Straw& s = allStraws[sid];
          const CLHEP::Hep3Vector& p = s.getMidPoint();
          const CLHEP::Hep3Vector& d = s.getDirection();
          double x = p.x();
          double y = p.y();
          double z = p.z();
          double l = s.halfLength();
          double st=sin(d.theta());
          double ct=cos(d.theta());
          double sp=sin(d.phi()+TMath::Pi()/2.0);
          double cp=cos(d.phi()+TMath::Pi()/2.0);
          if(sid < drawconfig.getInt("maxSID")){
            double x1=x+l*st*sp;
            double y1=y-l*st*cp;
            double z1=z+l*ct;
            double x2=x-l*st*sp;
            double y2=y+l*st*cp;
            double z2=z-l*ct;
            std::string strawtitle;
            int idStraw =  s.id().getStraw();
            int idPanel =  s.id().getPanel();
            int idPlane =  s.id().getPlane();
            int colorid = idPlane + idPanel;
            strawtitle =Form("Straw %i Panel %i  Plane %i",idStraw,idPanel,idPlane);
            sposi.set(x1, y1, z1);
            sposf.set(x2, y2, z2);
            if(sposi.x()!=0){
              auto strawline = new REX::REveLine("StrawHit ",strawtitle,2);
              strawline->SetPoint(0,pointmmTocm(sposi.x()),pointmmTocm(sposi.y()),pointmmTocm(sposi.z()));
              strawline->SetNextPoint(pointmmTocm(sposf.x()),pointmmTocm(sposf.y()),pointmmTocm(sposf.z()));
              strawline->SetLineWidth(1);
              strawline->SetLineColor(colorid);
              if(strawline->GetSize() !=0 ) scene->AddElement(strawline);
            }
          }
        }

        if(AddErrorBar_){
          //XY
          auto const& p = hit.pos();
          auto w = hit.uDir();
          auto const& s = hit.wireRes();
          double x1 = (p.x()+s*w.x());
          double x2 = (p.x()-s*w.x());
          double z1 = (p.z()+s*w.z());
          double z2 = (p.z()-s*w.z());
          double y1 = (p.y()+s*w.y());
          double y2 = (p.y()-s*w.y());

          std::string errorbar = Form("ErrorBar Length: %f, %f, %f, %f, %f, %f",x1,y1,z1,x2,y2,z2);
          auto error = new REX::REveLine("errors",errorbar.c_str(),2);
          error->SetPoint(0, pointmmTocm(x1),pointmmTocm(y1),pointmmTocm(z1));
          error->SetNextPoint(pointmmTocm(x2),pointmmTocm(y2),pointmmTocm(z2));
          //std::cout<<"points "<<p.x()<<" "<<p.y()<<" "<<p.z()<<std::endl;
          //std::cout<<"errors "<<x1<<" "<<y1<<" "<<z1<<std::endl;
          error->SetLineColor(kRed);
          error->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
          scene->AddElement(error);
        }
        CLHEP::Hep3Vector HitPos(pointmmTocm(hit.pos().x()), pointmmTocm(hit.pos().y()), pointmmTocm(hit.pos().z()));
        std::string chtitle = "ComboHits tag = "
          + (names[j])  +  '\n'
          + " position : x "  + std::to_string(hit.pos().x())  +  '\n'
          + " y " + std::to_string(hit.pos().y())  +  '\n'
          + " z " + std::to_string(hit.pos().z())  +  '\n'
          + " time :" + std::to_string(hit.time()) +  '\n'
          + " energy dep : "
          + std::to_string(hit.energyDep())  +
          + "MeV";
        auto ps1 = new REX::REvePointSet(chtitle, chtitle,0);
        ps1->SetNextPoint(HitPos.x(), HitPos.y() , HitPos.z());
        ps1->SetMarkerColor(colour);
        ps1->SetMarkerStyle(DataInterface::mstyle);
        ps1->SetMarkerSize(DataInterface::msize);
        if(ps1->GetSize() !=0 ) scene->AddElement(ps1);
      }
    }
  }
}

/*------------Function to add CRV information to the display:-------------*/
void DataInterface::AddCRVInfo(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CrvRecoPulseCollection*>>  crvpulse_tuple, REX::REveElement* &scene, bool extracted, bool addCRVBars){
  std::cout<<"[ REveDataInterface::AddCRVInfo() ]"<<std::endl;
  std::vector<const CrvRecoPulseCollection*> crvpulse_list = std::get<1>(crvpulse_tuple);
  std::vector<std::string> names = std::get<0>(crvpulse_tuple);
  GeomHandle<CosmicRayShield> CRS;
  GeomHandle<DetectorSystem> det;
  if(crvpulse_list.size() !=0){
    for(unsigned int i=0; i < crvpulse_list.size(); i++){
      const CrvRecoPulseCollection* crvRecoPulse = crvpulse_list[i];
      if(crvRecoPulse->size() !=0){
        std::string crvtitle = " CRV Bar Hit tag : " + names[i];
        auto ps1 = new REX::REvePointSet(crvtitle, crvtitle,0);
        auto allcrvbars = new REX::REveCompound("allcrvbars"+std::to_string(i),"allcrvbars"+std::to_string(i),1);
        for(unsigned int j=0; j< crvRecoPulse->size(); j++){

          mu2e::CrvRecoPulse const &crvpulse = (*crvRecoPulse)[j];

          const CRSScintillatorBarIndex &crvBarIndex = crvpulse.GetScintillatorBarIndex();
          const CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
          const CRSScintillatorBarDetail &barDetail = crvCounter.getBarDetail();
          CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
          CLHEP::Hep3Vector HitPos(crvCounterPos.x(), crvCounterPos.y(), crvCounterPos.z());

          CLHEP::Hep3Vector pointInMu2e = det-> toDetector(crvCounterPos);
          CLHEP::Hep3Vector sibardetails(barDetail.getHalfLengths()[0],barDetail.getHalfLengths()[1],barDetail.getHalfLengths()[2]);
          std::string pulsetitle = " CRV Bar Hit tag : "
          + names[i] +  '\n'
          + "CRVRecoPulse in Bar ID" +  '\n'
          + std::to_string(crvBarIndex.asInt());
          char const *pulsetitle_c = pulsetitle.c_str();
          if(addCRVBars){
            if(!extracted){
              auto b = new REX::REveBox(pulsetitle_c,pulsetitle_c);
              b->SetMainColor(drawconfig.getInt("CRVBarColor"));
              b-> SetMainTransparency(drawconfig.getInt("CRVtrans"));
              b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
              double  length = pointmmTocm(crvCounter.getHalfLength());
              double  width = pointmmTocm(crvCounter.getHalfWidth());
              double  height = pointmmTocm(crvCounter.getHalfThickness());

              if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 2 and barDetail.getLengthDirection() == 0){ //CRV D, CRV U -- WORKS
                b->SetVertex(0,  -1*length, pointmmTocm(pointInMu2e.y()) - width, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height); //---
                b->SetVertex(1,  -1*length, pointmmTocm(pointInMu2e.y()) - width, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height);//-+-
                b->SetVertex(2,  -1*length, pointmmTocm(pointInMu2e.y()) + width , pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//++-
                b->SetVertex(3,  -1*length, pointmmTocm(pointInMu2e.y()) + width  ,pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height);//+--
                b->SetVertex(4, length, pointmmTocm(pointInMu2e.y()) - width,   pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height  );//--+
                b->SetVertex(5,  length, pointmmTocm(pointInMu2e.y()) - width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//-++
                b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//+++
                b->SetVertex(7, length,pointmmTocm(pointInMu2e.y()) + width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height );//+-+

              }

              if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ //CRV TS -- WORKS
                b->SetVertex(0, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + -1*length ); //---
                b->SetVertex(1, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + -1*length);//-+-
                b->SetVertex(2, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + -1*length);//++-
                b->SetVertex(3, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + -1*length);//+--
                b->SetVertex(4, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + length);//--+
                b->SetVertex(5, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + length);//-++
                b->SetVertex(6, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + length);//+++
                b->SetVertex(7, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + length);//+-+
              }

              if(barDetail.getWidthDirection() == 2 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 0){ //CRV T -- orientation DONE
                b->SetVertex(0, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width); //---
                b->SetVertex(1, -1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width);//-+-
                b->SetVertex(2, -1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width);//++-
                b->SetVertex(3, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width);//+--
                b->SetVertex(4, length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width );//--+
                b->SetVertex(5, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width );//-++
                b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width );//+++
                b->SetVertex(7, length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width );//+-+
              }

              if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 0 and barDetail.getLengthDirection() == 2){ //CRV R, CRV L -- orientation DONE
                b->SetVertex(0, pointmmTocm(pointInMu2e.y()) - height, -1*length, pointmmTocm(pointInMu2e.x()) - width );//---
                b->SetVertex(1, pointmmTocm(pointInMu2e.y()) + height, -1*length,pointmmTocm(pointInMu2e.x()) - width);//-+-
                b->SetVertex(2, pointmmTocm(pointInMu2e.y()) + height, -1*length,pointmmTocm(pointInMu2e.x()) + width);//++-
                b->SetVertex(3, pointmmTocm(pointInMu2e.y()) - height, -1*length,pointmmTocm(pointInMu2e.x()) + width);//+--
                b->SetVertex(4, pointmmTocm(pointInMu2e.y()) - height, length, pointmmTocm(pointInMu2e.x()) - width);//--+
                b->SetVertex(5, pointmmTocm(pointInMu2e.y()) + height, length, pointmmTocm(pointInMu2e.x()) - width );//-++
                b->SetVertex(6, pointmmTocm(pointInMu2e.y()) + height, length,pointmmTocm(pointInMu2e.x()) + width );//+++
                b->SetVertex(7, pointmmTocm(pointInMu2e.y()) - height, length,pointmmTocm(pointInMu2e.x()) + width );//+-+
              }
              allcrvbars->AddElement(b);
            }
            if(extracted){ //TODO same for nominal geom

              // CRV hit scintillation bars highlighted
              // std::string const& base;
              // std::string bartitle = crvCounter.name( base );
              // char const *bartitle_c = base.c_str(); //TODO title

              // Draw "bars hit" in red:
              auto b = new REX::REveBox(pulsetitle_c,pulsetitle_c);
              b->SetMainColor(drawconfig.getInt("CRVBarColor"));
              b-> SetMainTransparency(drawconfig.getInt("CRVtrans"));
              b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
              double  length = pointmmTocm(crvCounter.getHalfLength());
              double  width = pointmmTocm(crvCounter.getHalfWidth());
              double  height = pointmmTocm(crvCounter.getHalfThickness());

              if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ //T1

                b->SetVertex(0, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, -1*length+ pointmmTocm(pointInMu2e.z())); //---
                b->SetVertex(1, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y())  + height, -1*length+ pointmmTocm(pointInMu2e.z()));//-+-
                b->SetVertex(2, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, -1*length+ pointmmTocm(pointInMu2e.z()));//++-
                b->SetVertex(3, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, -1*length+ pointmmTocm(pointInMu2e.z()));//+--
                b->SetVertex(4, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height,  length + pointmmTocm(pointInMu2e.z()));//--+
                b->SetVertex(5, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, length + pointmmTocm(pointInMu2e.z()));//-++
                b->SetVertex(6, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height,length + pointmmTocm(pointInMu2e.z()));//+++
                b->SetVertex(7, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, length+ pointmmTocm(pointInMu2e.z()) );//+-+
              } else { //EX, T2
                b->SetVertex(0, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()));//---
                b->SetVertex(1,-1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()));//-+-
                b->SetVertex(2, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.x()) - width  + pointmmTocm(pointInMu2e.z()));//++-
                b->SetVertex(3, length, pointmmTocm(pointInMu2e.y()) - height , pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()) );//+--
                b->SetVertex(4, -1*length , pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z())); //--+
                b->SetVertex(5,-1*length, pointmmTocm(pointInMu2e.y()) + height , pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()) ); // -++
                b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + height , pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()));//+++
                b->SetVertex(7,  length , pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()));//+-+

              }
              scene->AddElement(b);
            }
          }
          // Add Reco Pulse to collection
          ps1->SetNextPoint(pointmmTocm(pointInMu2e.x()), pointmmTocm(pointInMu2e.y()) , pointmmTocm(pointInMu2e.z()));

        }
        scene->AddElement(allcrvbars);
        // Draw reco pulse collection
        ps1->SetMarkerColor(i+3);
        ps1->SetMarkerStyle(DataInterface::mstyle);
        ps1->SetMarkerSize(DataInterface::msize);
        if(ps1->GetSize() !=0 ) scene->AddElement(ps1);
      }
    }
  }

}


/*------------Function to add CRV information to the display:-------------*/
void DataInterface::AddCRVClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const CrvCoincidenceClusterCollection*>>  crvpulse_tuple, REX::REveElement* &scene, bool extracted, bool addCRVBars){

  std::vector<const CrvCoincidenceClusterCollection*> crvpulse_list = std::get<1>(crvpulse_tuple);
  std::vector<std::string> names = std::get<0>(crvpulse_tuple);
  GeomHandle<CosmicRayShield> CRS;
  GeomHandle<DetectorSystem> det;
  if(crvpulse_list.size() !=0){
    for(unsigned int i=0; i < crvpulse_list.size(); i++){
      const CrvCoincidenceClusterCollection* crvClusters = crvpulse_list[i];
      if(crvClusters->size() !=0){

        for(unsigned int j=0; j< crvClusters->size(); j++){
          auto allcrvbars = new REX::REveCompound("Bars for CRVConicidenceCluster"+std::to_string(j),"Bars for CRVConicidenceCluster"+std::to_string(j),1);
          mu2e::CrvCoincidenceCluster const &crvclu = (*crvClusters)[j];
          std::string crvtitle = "CRVCoincidenceCluster" + std::to_string(j) + " tag : " + names[i] + '\n'
          + " averge hit time = " + std::to_string(crvclu.GetAvgHitTime())+" ns " + '\n'
          + " PEs = " + std::to_string(crvclu.GetPEs());
          auto ps1 = new REX::REvePointSet(crvtitle, crvtitle,0);
          CLHEP::Hep3Vector pointInMu2e = det-> toDetector(crvclu.GetAvgHitPos());
          ps1->SetNextPoint(pointmmTocm(pointInMu2e.x()), pointmmTocm(pointInMu2e.y()) , pointmmTocm(pointInMu2e.z()));
          for(unsigned h =0 ; h < crvclu.GetCrvRecoPulses().size();h++)     {

            art::Ptr<CrvRecoPulse>  crvpulse = crvclu.GetCrvRecoPulses()[h];
            const mu2e::CRSScintillatorBarIndex &crvBarIndex = crvpulse->GetScintillatorBarIndex();
            const mu2e::CRSScintillatorBar &crvCounter = CRS->getBar(crvBarIndex);
            const mu2e::CRSScintillatorBarDetail &barDetail = crvCounter.getBarDetail();
            CLHEP::Hep3Vector crvCounterPos = crvCounter.getPosition();
            CLHEP::Hep3Vector HitPos(crvCounterPos.x(), crvCounterPos.y(), crvCounterPos.z());

            CLHEP::Hep3Vector pointInMu2e = det-> toDetector(crvCounterPos);
            CLHEP::Hep3Vector sibardetails(barDetail.getHalfLengths()[0],barDetail.getHalfLengths()[1],barDetail.getHalfLengths()[2]);
            std::string pulsetitle = " CRV Bar Hit for  tag : "
            + names[i] +  '\n'
            + "Bar ID"
            + std::to_string(crvBarIndex.asInt())+  '\n'
            + "Coincidence start time " + std::to_string(crvclu.GetStartTime())+  '\n'
            + "Coincidence end time " + std::to_string(crvclu.GetEndTime());
            char const *pulsetitle_c = pulsetitle.c_str();
            if(addCRVBars){
              if(!extracted){
                auto b = new REX::REveBox(pulsetitle_c,pulsetitle_c);
                b->SetMainColor(drawconfig.getInt("CRVBarColor"));
                b-> SetMainTransparency(drawconfig.getInt("CRVtrans"));
                b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
                double  length = pointmmTocm(crvCounter.getHalfLength());
                double  width = pointmmTocm(crvCounter.getHalfWidth());
                double  height = pointmmTocm(crvCounter.getHalfThickness());

                if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 2 and barDetail.getLengthDirection() == 0){ //CRV D, CRV U -- WORKS
                  b->SetVertex(0,  -1*length, pointmmTocm(pointInMu2e.y()) - width, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height); //---
                  b->SetVertex(1,  -1*length, pointmmTocm(pointInMu2e.y()) - width, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height);//-+-
                  b->SetVertex(2,  -1*length, pointmmTocm(pointInMu2e.y()) + width , pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//++-
                  b->SetVertex(3,  -1*length, pointmmTocm(pointInMu2e.y()) + width  ,pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height);//+--
                  b->SetVertex(4, length, pointmmTocm(pointInMu2e.y()) - width,   pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height  );//--+
                  b->SetVertex(5,  length, pointmmTocm(pointInMu2e.y()) - width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//-++
                  b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + height );//+++
                  b->SetVertex(7, length,pointmmTocm(pointInMu2e.y()) + width,  pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - height );//+-+

                }

                if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ //CRV TS -- WORKS
                  b->SetVertex(0, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + -1*length ); //---
                  b->SetVertex(1, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + -1*length);//-+-
                  b->SetVertex(2, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + -1*length);//++-
                  b->SetVertex(3, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + -1*length);//+--
                  b->SetVertex(4, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + length);//--+
                  b->SetVertex(5, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + length);//-++
                  b->SetVertex(6, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + length);//+++
                  b->SetVertex(7, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + length);//+-+

                }

                if(barDetail.getWidthDirection() == 2 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 0){ //CRV T -- WORKS
                  b->SetVertex(0, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width); //---
                  b->SetVertex(1, -1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width);//-+-
                  b->SetVertex(2, -1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width);//++-
                  b->SetVertex(3, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width);//+--
                  b->SetVertex(4, length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width );//--+
                  b->SetVertex(5, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) - width );//-++
                  b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width );//+++
                  b->SetVertex(7, length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.z()) + pointmmTocm(pointInMu2e.x()) + width );//+-+

                }

                if(barDetail.getWidthDirection() == 1 and barDetail.getThicknessDirection() == 0 and barDetail.getLengthDirection() == 2){ //CRV R, CRV L -- orientation DONE
                  b->SetVertex(0, pointmmTocm(pointInMu2e.y()) - height, -1*length, pointmmTocm(pointInMu2e.x()) - width ); //---
                  b->SetVertex(1, pointmmTocm(pointInMu2e.y()) + height, -1*length,pointmmTocm(pointInMu2e.x()) - width);//-+-
                  b->SetVertex(2, pointmmTocm(pointInMu2e.y()) + height, -1*length,pointmmTocm(pointInMu2e.x()) + width);//++-
                  b->SetVertex(3, pointmmTocm(pointInMu2e.y()) - height, -1*length,pointmmTocm(pointInMu2e.x()) + width);//+--
                  b->SetVertex(4, pointmmTocm(pointInMu2e.y()) - height, length, pointmmTocm(pointInMu2e.x()) - width);//--+
                  b->SetVertex(5, pointmmTocm(pointInMu2e.y()) + height, length, pointmmTocm(pointInMu2e.x()) - width );//-++
                  b->SetVertex(6, pointmmTocm(pointInMu2e.y()) + height, length,pointmmTocm(pointInMu2e.x()) + width );//+++
                  b->SetVertex(7, pointmmTocm(pointInMu2e.y()) - height, length,pointmmTocm(pointInMu2e.x()) + width );//+-+

                }
                allcrvbars->AddElement(b);

              }
              if(extracted){

                // Draw "bars hit" in red:
                auto b = new REX::REveBox(pulsetitle_c,pulsetitle_c);
                b->SetMainColor(drawconfig.getInt("CRVBarColor"));
                b-> SetMainTransparency(drawconfig.getInt("CRVtrans"));
                b->SetLineWidth(drawconfig.getInt("GeomLineWidth"));
                double  length = pointmmTocm(crvCounter.getHalfLength());
                double  width = pointmmTocm(crvCounter.getHalfWidth());
                double  height = pointmmTocm(crvCounter.getHalfThickness());

                if(barDetail.getWidthDirection() == 0 and barDetail.getThicknessDirection() == 1 and barDetail.getLengthDirection() == 2){ //T1

                  b->SetVertex(0, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height, -1*length + pointmmTocm(pointInMu2e.z())); //---
                  b->SetVertex(1, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y())  + height, -1*length + pointmmTocm(pointInMu2e.z()));//-+-
                  b->SetVertex(2, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height, -1*length + pointmmTocm(pointInMu2e.z()));//++-
                  b->SetVertex(3, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, -1*length + pointmmTocm(pointInMu2e.z()));//+--
                  b->SetVertex(4, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) - height,  length + pointmmTocm(pointInMu2e.z()));//--+
                  b->SetVertex(5, pointmmTocm(pointInMu2e.x()) - width, pointmmTocm(pointInMu2e.y()) + height, length + pointmmTocm(pointInMu2e.z()));//-++
                  b->SetVertex(6, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) + height ,length + pointmmTocm(pointInMu2e.z()) );//+++
                  b->SetVertex(7, pointmmTocm(pointInMu2e.x()) + width, pointmmTocm(pointInMu2e.y()) - height, length + pointmmTocm(pointInMu2e.z()));//+-+
                } else { //EX, T2
                  b->SetVertex(0, -1*length, pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()));//---
                  b->SetVertex(1,-1*length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()));//-+-
                  b->SetVertex(2, length, pointmmTocm(pointInMu2e.y()) + height, pointmmTocm(pointInMu2e.x()) - width  + pointmmTocm(pointInMu2e.z()));//++-
                  b->SetVertex(3, length, pointmmTocm(pointInMu2e.y()) - height , pointmmTocm(pointInMu2e.x()) - width + pointmmTocm(pointInMu2e.z()) );//+--
                  b->SetVertex(4, -1*length , pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z())); //--+
                  b->SetVertex(5,-1*length, pointmmTocm(pointInMu2e.y()) + height , pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()) ); // -++
                  b->SetVertex(6, length, pointmmTocm(pointInMu2e.y()) + height , pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()));//+++
                  b->SetVertex(7,  length , pointmmTocm(pointInMu2e.y()) - height, pointmmTocm(pointInMu2e.x()) + width + pointmmTocm(pointInMu2e.z()));//+-+

                }
                scene->AddElement(b);
              }
            }
          }
        ps1->SetMarkerColor(drawconfig.getInt("CRVHitColor"));
        ps1->SetMarkerStyle(DataInterface::mstyle);
        ps1->SetMarkerSize(DataInterface::msize);
        if(ps1->GetSize() !=0 ) scene->AddElement(ps1);
        scene->AddElement(allcrvbars);
        }
      }
    }
  }

}

/*------------Function to add TimeCluster Collection in 3D and 2D displays:-------------*/
void DataInterface::AddTimeClusters(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const TimeClusterCollection*>>  timecluster_tuple, REX::REveElement* &scene){

  std::vector<const TimeClusterCollection*> timecluster_list = std::get<1>(timecluster_tuple);
  std::vector<std::string> names = std::get<0>(timecluster_tuple);

  if(timecluster_list.size() !=0){
    for(unsigned int i = 0; i < timecluster_list.size(); i++){
      const TimeClusterCollection* tccol = timecluster_list[i];
      if(tccol->size() != 0){
        if(!firstLoop_){
          scene->DestroyElements();;
        }


        for(size_t j=0; j<tccol->size();j++){
          mu2e::TimeCluster const  &tclust= (*tccol)[j];
          std::string tctitle = "Time Cluster tag: " + names[i] + '\n'
            + "t0 " + std::to_string(tclust.t0().t0()) + " +/- " + std::to_string(tclust.t0().t0Err()) + " ns " + '\n' ;
          auto ps1 = new REX::REvePointSet("TimeClusters", tctitle, 0);
          CLHEP::Hep3Vector HitPos(tclust._pos.x(), tclust._pos.y(), tclust._pos.z());
          ps1->SetNextPoint(pointmmTocm(HitPos.x()), pointmmTocm(HitPos.y()) , pointmmTocm(HitPos.z()));
          ps1->SetMarkerColor(i+6);
          ps1->SetMarkerStyle(kOpenCircle);
          ps1->SetMarkerSize(DataInterface::msize);
          if(ps1->GetSize() !=0 ) scene->AddElement(ps1);
        }
      }
    }
  }
}


void DataInterface::AddHelixSeedCollection(REX::REveManager *&eveMng,bool firstloop,  std::tuple<std::vector<std::string>, std::vector<const HelixSeedCollection*>> helix_tuple, REX::REveElement* &scene){
  std::cout<<"[DataInterface] AddHelices "<<std::endl;
  std::vector<const HelixSeedCollection*> helix_list = std::get<1>(helix_tuple);
  std::vector<std::string> names = std::get<0>(helix_tuple);
  for(unsigned int j=0; j< helix_list.size(); j++){
    const HelixSeedCollection* seedcol = helix_list[j];
    if(seedcol!=0){
      for(unsigned int k = 0; k < seedcol->size(); k++){
        mu2e::HelixSeed hseed = (*seedcol)[k];
        const ComboHitCollection& hhits = hseed.hits();
        unsigned int nhhits = hhits.size();
        auto line = new REX::REveLine(names[j], names[j],nhhits);
        // std::cout<<"Helix hits = "<<nhhits<<" radius = "<<hseed.helix()._radius<<" x0 = "<<hseed.helix()._rcent*cos(hseed.helix()._fcent)<<std::endl;
        float helrad = hseed.helix()._radius;
        float circphi = 0.0;
        float xc = hseed.helix()._rcent*cos(hseed.helix()._fcent);
        float yc = hseed.helix()._rcent*sin(hseed.helix()._fcent);
        if(hhits.size() !=0 ){
          for(int i=-1500; i < 1500; i +=100){ //TODO Remove hard-coded numbers
            float z = i;
            if(hseed.helix()._lambda !=0.0) circphi = hseed.helix()._fz0 + z/hseed.helix()._lambda;
            float x= xc + helrad*cos(circphi);
            float y = yc + helrad*sin(circphi);
            //std::cout<<"x = "<<x<<" y = "<<y<<" z = "<<z<<std::endl;
            CLHEP::Hep3Vector HelPos(x, y, z);
            if(circphi !=0.0)line->SetNextPoint(pointmmTocm(HelPos.x()),pointmmTocm(HelPos.y()) ,pointmmTocm(HelPos.z()));
          }
        }
        line->SetLineColor(drawconfig.getInt("RecoTrackColor"));
        line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
        scene->AddElement(line);
      }
    }
  }
}

void DataInterface::AddKalIntersection(KalSeed const& kalseed, REX::REveElement* &scene, REX::REveCompound *products){
  //Plot intersecitons:
  std::vector<mu2e::KalIntersection> const& inters = kalseed.intersections();

  for(auto const& inter : inters){
    KinKal::VEC3 posKI = inter.position3();
    std::string title = "KalIntersection position : x "  + std::to_string(posKI.x())  +  '\n'
      + " y " + std::to_string(posKI.y())  +  '\n'
      + " z " + std::to_string(posKI.z())  +  '\n'
      + " time :" + std::to_string(inter.time()) +  '\n'
      + " mom , dmom : "
      + std::to_string(inter.mom()) + " , " + std::to_string(inter.dMom())  +
      + "MeV/c " + '\n'
      + "Surface " +  inter.surfaceId().name();
    auto interpoint = new REX::REvePointSet(title, title,1);

    interpoint->SetMarkerStyle(DataInterface::mstyle);
    interpoint->SetMarkerSize(DataInterface::msize);
    if(fabs(inter.dMom()) > 0.0){
      interpoint->SetMarkerColor(kViolet);// color material intersections differently from virtual surface intersections
    } else {
      interpoint->SetMarkerColor(kYellow);
    }
    interpoint->SetNextPoint(pointmmTocm(posKI.x()),pointmmTocm(posKI.y()) ,pointmmTocm(posKI.z()));
    if(interpoint->GetSize() !=0 ) products->AddElement(interpoint);
  }

}

template<class KTRAJc> void DataInterface::AddTrkStrawHit(KalSeed const& kalseed, REX::REveElement* &scene,  std::unique_ptr<KTRAJc> &lhptr, REX::REveCompound *trackproducts){
  std::cout<<"[DataInterface::AddTrkStrawHit]"<<std::endl;
  //Plot trk straw hits
  mu2e::GeomHandle<mu2e::Tracker> tracker;
  std::vector<mu2e::TrkStrawHitSeed> const& hits = kalseed.hits();

  for(unsigned int i = 0; i < hits.size(); i++){
    const mu2e::TrkStrawHitSeed &tshs = hits.at(i);
    auto const& straw = tracker->straw(tshs.strawId());
    mu2e::WireHitState whs(mu2e::WireHitState::State(tshs._ambig),
        mu2e::StrawHitUpdaters::algorithm(tshs._algo),
        tshs._kkshflag);
    bool active = whs.active();
    bool usedrift = whs.driftConstraint();
    if(active){ // maybe draw inactive hits but with a different color? TODO
      // then find the position at the reference POCA: start with the position on the wire
      auto tshspos = XYZVectorF(straw.wirePosition(tshs._rupos));
      double herr(0.0);
      //find the direction along DOCA
      auto tdir = lhptr->direction(tshs._ptoca);
      auto wdir = XYZVectorF(straw.wireDirection(tshs._rupos));
      auto ddir = wdir.Cross(tdir).Unit()*whs.lrSign();
      if(usedrift){
        // move the position out on the signed drift direction.
        tshspos += tshs._rdrift*ddir;
        herr = tshs._sderr;
      } else {
        herr = tshs._uderr;
      }
      // set the line length to be N sigma.  1 may be too short to see, TODO
      double nsigma(2.0);
      auto end1 = tshspos + nsigma*herr*ddir;
      auto end2 = tshspos - nsigma*herr*ddir;
      std::string err_title = "+/-"+std::to_string(nsigma) +"sigma";

      //goes along that same line (ddir)
      std::string title = "TrkStrawHitSeed : x "  + std::to_string(tshspos.x())  +  '\n'
        + " y " + std::to_string(tshspos.y())  +  '\n'
        + " z " + std::to_string(tshspos.z())  +  '\n'
        + " time :" + std::to_string(tshs.hitTime())+  '\n'
        + " energyDep :" + std::to_string(tshs.energyDep())+ "MeV" + '\n'
        + " error : " + err_title;
      auto point_with_error = new REX::REveCompound("TrkStrawHitSeed "+std::to_string(i), "TrkStrawHitSeed",1);
      auto trkstrawpoint = new REX::REvePointSet(title, title,1);
      trkstrawpoint->SetMarkerStyle(DataInterface::mstyle);
      trkstrawpoint->SetMarkerSize(DataInterface::msize);
      trkstrawpoint->SetMarkerColor(drawconfig.getInt("TrkHitColor"));
      auto line = new REX::REveLine("TrkStrawHit Error"+err_title,err_title, 1);
      line->SetNextPoint(pointmmTocm(end1.x()),pointmmTocm(end1.y()) ,pointmmTocm(end1.z()));
      line->SetNextPoint(pointmmTocm(end2.x()),pointmmTocm(end2.y()) ,pointmmTocm(end2.z()));
      line->SetLineColor(drawconfig.getInt("TrkHitColor"));
      if(!usedrift)trkstrawpoint->SetMarkerColor(drawconfig.getInt("TrkNoHitColor"));
      trkstrawpoint->SetNextPoint(pointmmTocm(tshspos.x()),pointmmTocm(tshspos.y()) ,pointmmTocm(tshspos.z()));
      point_with_error->AddElement(trkstrawpoint);
      point_with_error->AddElement(line);
      trackproducts->AddElement(point_with_error);
    }
  }
}

void DataInterface::AddTrkCaloHit(KalSeed const& kalseed, REX::REveElement* &scene){
  /*mu2e::TrkCaloHitSeed caloseed = kalseed.caloHit();
    art::Ptr<mu2e::CaloCluster> cluster = caloseed.caloCluster(); // TODO - use this to access x,y,z,t
    */
}

/*------------Function to color code the Tracker hits -------------*/
void DataInterface::AddTrkHits(REX::REveManager *&eveMng, bool firstLoop_, std::tuple<std::vector<std::string>, std::vector<const ComboHitCollection*>> combohit_tuple,std::tuple<std::vector<std::string>, std::vector<const KalSeedPtrCollection*>> track_tuple, REX::REveElement* &scene){
  std::cout<<"[DataInterface] AddTrkHits  "<<std::endl;
  std::vector<const ComboHitCollection*> combohit_list = std::get<1>(combohit_tuple);
  std::vector<const KalSeedPtrCollection*> track_list = std::get<1>(track_tuple);

  GeomHandle<DetectorSystem> det;
  std::vector<StrawId> trksid(drawconfig.getInt("maxStrawID"));
  unsigned int trkhitsize=0;
  //Save the hit straw IDs of the KalSeed hits
  for(unsigned int j = 0; j< track_list.size(); j++){
    const KalSeedPtrCollection* seedcol = track_list[j];
    if(seedcol!=0){
      for(unsigned int k = 0; k < seedcol->size(); k++){
        auto const& kseed = *(*seedcol)[k];
        const std::vector<mu2e::TrkStrawHitSeed> &hits = kseed.hits();
        trkhitsize = hits.size();
        for(unsigned int i = 0; i <trkhitsize; i++){
          const mu2e::TrkStrawHitSeed &hit = hits.at(i);
          trksid[i] = hit._sid;
        }
        std::vector<StrawId> usedtrksid(trkhitsize);
        std::vector<unsigned int> usedid(trkhitsize);
        //Compare the straw IDs of the Kal seed hits with the hits in the ComboHit Collection
        for(unsigned int j=0; j< combohit_list.size(); j++){
          const ComboHitCollection* chcol = combohit_list[j]; //TODO shouldnt rely on ComboHits
          if(chcol!=0){
            for(unsigned int i=0; i<chcol->size();i++){
              ComboHit hit = (*chcol)[i];
              for(unsigned int q=0; q<trkhitsize; q++){
                if(hit._sid == trksid[q]){
                  usedtrksid[q]=hit._sid;//Save the Straw ID if the KalSeed and Combo hit ID matches
                  usedid[q]=q;
                  CLHEP::Hep3Vector HitPos(hit.pos().x(), hit.pos().y(), hit.pos().z());
                  std::string chtitle = "TrkSeedHit position : x "  + std::to_string(hit.pos().x())  +  '\n'
                    + " y " + std::to_string(hit.pos().y())  +  '\n'
                    + " z " + std::to_string(hit.pos().z())  +  '\n'
                    + " time :" + std::to_string(hit.time()) +  '\n'
                    + " energy dep : "
                    + std::to_string(hit.energyDep())  +
                    + "MeV";

                  auto trkhit = new REX::REvePointSet(chtitle, chtitle,0);
                  trkhit ->SetMarkerStyle(DataInterface::mstyle);
                  trkhit ->SetMarkerSize(DataInterface::msize);
                  // trkhit ->SetMarkerColor(drawconfig.getInt("RecoTrackColor")-4);
                  trkhit ->SetMarkerColor(drawconfig.getInt("TrkHitColor"));
                  trkhit ->SetNextPoint(pointmmTocm(HitPos.x()),pointmmTocm(HitPos.y()) ,pointmmTocm(HitPos.z()));
                  if(trkhit->GetSize() !=0 ) scene->AddElement(trkhit);
                  // std::cout<<"TrkHit = "<<HitPos.x()<<"  "<<HitPos.y()<<" "<<HitPos.z()<<std::endl;
                }
              }
            }
          }
        }
      }
    }
  }
}

using LHPT = KalSeed::LHPT;
using CHPT = KalSeed::CHPT;
using KLPT = KalSeed::KLPT;
template<class KTRAJ> void DataInterface::AddKinKalTrajectory( std::unique_ptr<KTRAJ> &trajectory, REX::REveElement* &scene, unsigned int j, std::string kaltitle, double& t1, double& t2){
  t1=trajectory->range().begin();
  t2=trajectory->range().end();

  double x1=trajectory->position3(t1).x();
  double y1=trajectory->position3(t1).y();
  double z1=trajectory->position3(t1).z();

  auto line = new REX::REveLine(kaltitle,kaltitle, 100);
  line->SetPoint(0,pointmmTocm(x1), pointmmTocm(y1) , pointmmTocm(z1));
  for(double t=t1; t<=t2; t+=0.1)
  {
    const auto &p = trajectory->position3(t);
    double xt=p.x();
    double yt=p.y();
    double zt=p.z();
    line->SetNextPoint(pointmmTocm(xt), pointmmTocm(yt) , pointmmTocm(zt));
  }
  line->SetLineColor(j+6);
  line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
  scene->AddElement(line);
}

void DataInterface::FillKinKalTrajectory(REX::REveManager *&eveMng, bool firstloop, REX::REveElement* &scene, std::tuple<std::vector<std::string>, std::vector<const KalSeedPtrCollection*>> track_tuple, bool plotKalIntersection, bool addTrkHits, double& t1, double& t2){
  std::cout<<"[DataInterface::FillKinKalTrajectory()]"<<std::endl;
  auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
  std::vector<const KalSeedPtrCollection*> track_list = std::get<1>(track_tuple);
  std::vector<std::string> names = std::get<0>(track_tuple);
  for(unsigned int j=0; j< track_list.size(); j++){
    const KalSeedPtrCollection* seedcol = track_list[j];
    if(seedcol!=0){
      for(auto const& kseedptr : *seedcol){
        auto const& kseed = *kseedptr;
//        unsigned nhits = kseed.hits().size();
        unsigned nactive =0;
        for (auto const& hit : kseed.hits()){ if (hit.strawHitState() > WireHitState::inactive) ++nactive; }
        // use t0 to define the reference segment
        double t0;
        kseed.t0Segment(t0);
        std::stringstream ksstream;
        REX::REveCompound *trackproducts = new REX::REveCompound("Track Products for KalSeed "+ ptable->particle(kseed.particle()).name() + " LoopHelix ","Track Products for KalSeed "+ ptable->particle(kseed.particle()).name() + " LoopHelix " ,1);
        if(kseed.loopHelixFit())
        {
          auto trajectory=kseed.loopHelixFitTrajectory();
          auto const& lh = trajectory->nearestPiece(t0);
          auto momvec = lh.momentum3(t0);
          ksstream << ptable->particle(kseed.particle()).name() <<  " LoopHelix "
            << std::setw(6) << std::setprecision(3)
            <<  momvec.R() << " MeV/c, cos(Theta) " << cos(momvec.Theta()) << std::endl
            << "t0 " << lh.t0() << " ns "
            << "lam "  << lh.lam() << " mm "
            << "rad "  << lh.rad() << " mm " << std::endl
            << "cx "  << lh.cx() << " mm "
            << "cy "  << lh.cy() << " mm "
            << "phi0 "  << lh.phi0() << " rad " << std::endl
            << "N hits " <<  nactive << " fit cons. " << kseed.fitConsistency() << std::endl
            << "Instance " << names[j] << std::endl;
          AddKinKalTrajectory<LHPT>(trajectory,scene,j, ksstream.str(), t1, t2);
          if(addTrkHits) {
            AddTrkStrawHit<LHPT>(kseed, scene, trajectory,trackproducts);
          }
        }
        else if(kseed.centralHelixFit())
        {
          auto trajectory=kseed.centralHelixFitTrajectory();
          auto const& ch = trajectory->nearestPiece(t0);
          auto momvec = ch.momentum3(t0);
          ksstream << ptable->particle(kseed.particle()).name() << " CentralHelix "
            << std::setw(6) << std::setprecision(3)
            <<  momvec.R() << " MeV/c, cos(Theta) " << cos(momvec.Theta()) << std::endl
            << "t0 " << ch.t0() << " ns "
            << "tandip " << ch.tanDip( )
            << "d0 " << ch.d0( ) << " mm " << std::endl
            << "z0 " << ch.z0( ) << " mm "
            << "phi0 " << ch.phi0( ) << " rad "
            << "omega " << ch.omega( )<< " mm^-1 " << std::endl
            << "track arrival time " << t1 << std::endl
            << "Instance " << names[j] << std::endl;
          AddKinKalTrajectory<CHPT>(trajectory,scene,j, ksstream.str(), t1, t2);
          if(addTrkHits) {
            AddTrkStrawHit<CHPT>(kseed, scene, trajectory,trackproducts);
            }
        }
        else if(kseed.kinematicLineFit())
        {
          auto trajectory=kseed.kinematicLineFitTrajectory();
          auto const& kl = trajectory->nearestPiece(t0);
          auto momvec = kl.momentum3(t0);
          ksstream << ptable->particle(kseed.particle()).name() << " KinematicLine "
            << std::setw(6) << std::setprecision(3)
            <<  momvec.R() << " MeV/c, cos(Theta) " << cos(momvec.Theta()) << std::endl
            << " t0 " << kl.t0() << " ns "
            << " d0 " << kl.d0( ) << " mm "
            << " z0 " << kl.z0( ) << " mm " << std::endl
            << " phi0 " << kl.phi0( )<< " rad "
            << " theta " << kl.theta( ) << " rad " <<  std::endl
            << " track arrival time " << t1 << std::endl
            << " Instance " << names[j] << std::endl;
          AddKinKalTrajectory<KLPT>(trajectory,scene,j, ksstream.str(), t1, t2);
          if(addTrkHits) {
            AddTrkStrawHit<KLPT>(kseed, scene, trajectory, trackproducts );
          }
        }
        if(plotKalIntersection) AddKalIntersection(kseed, scene,trackproducts );
        if(addTrkHits or plotKalIntersection) {
        scene->AddElement(trackproducts);
        }
      }
    }
  }
}

void DataInterface::AddCosmicTrackFit(REX::REveManager *&eveMng, bool firstLoop_, const mu2e::CosmicTrackSeedCollection *cosmiccol, REX::REveElement* &scene){

  std::cout<<"[DataInterface] AddCosmicTrackSeed "<<std::endl;
  if(cosmiccol!=0){
    auto line = new REX::REveLine("Cosmic","Cosmic",2);
    for(unsigned int i=0; i< cosmiccol->size(); i++){
      mu2e::CosmicTrackSeed const  &sts= (*cosmiccol)[i];
      mu2e::CosmicTrack st = sts._track;
      double ty1 = sts._straw_chits[0].pos().y();
      double ty2 = sts._straw_chits[sts._straw_chits.size()-1].pos().y();
      double tx1 = st.MinuitParams.A0  - st.MinuitParams.A1*ty1;
      double tx2 = st.MinuitParams.A0  - st.MinuitParams.A1*ty2;
      double tz1 = st.MinuitParams.B0  - st.MinuitParams.B1*ty1;
      double tz2 = st.MinuitParams.B0  - st.MinuitParams.B1*ty2;
      line->SetNextPoint(pointmmTocm(tx1), pointmmTocm(ty1) , pointmmTocm(tz1));
      line->SetNextPoint(pointmmTocm(tx2), pointmmTocm(ty2) , pointmmTocm(tz2));
    }

    line->SetLineColor(drawconfig.getInt("RecoTrackColor"));
    line->SetLineWidth(drawconfig.getInt("TrackLineWidth"));
    scene->AddElement(line);
  }

}

