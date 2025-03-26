#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "EventDisplay/inc/MainWindow.hh"

namespace REX = ROOT::Experimental;
using namespace std;
using namespace mu2e;

//positions of components for DS component offsets
double x_d0 = 0; double x_d1 = 0; double x_cal = 0;double x_trk = 0;double x_ds3 = 0;double x_hall = 0;double x_world = 0;double x_st = 0; double x_ds2 = 0; double x_ipa =0;
double y_d0 = 0; double y_d1 = 0; double y_cal  = 0;double y_trk = 0;double y_ds3 = 0;double y_hall = 0;double y_world = 0;double y_st = 0; double y_ds2 = 0; double y_ipa =0;
double z_d0 = 0; double z_d1 = 0; double z_cal = 0;double z_trk = 0;double z_ds3 = 0;double z_hall = 0;double z_world = 0;double z_st = 0; double z_ds2 = 0; double z_ipa =0;
double x_crv = 0 ; double y_crv = 0 ; double z_crv = 0;  double z_crv_u = 0;  double z_crv_d = 0;
double nCrystals = 674;
double FrontTracker_gdmltag;

// Get geom information from Event Display configs:
std::string drawoptfilename("EventDisplay/config/drawutils.txt");
SimpleConfig drawconfigf(drawoptfilename);

void MainWindow::makeEveGeoShape(TGeoNode* n, REX::REveTrans& trans, REX::REveElement* holder, int val, bool crystal1, bool crystal2, std::string name, int color)
{
  bool isMother = name.find("ProductionTargetMother") != string::npos;
  auto gss = n->GetVolume()->GetShape();
  auto b1s = new REX::REveGeoShape(n->GetName());

  b1s->InitMainTrans();
  b1s->RefMainTrans().SetFrom(trans.Array());
  b1s->SetShape(gss);

  int transpar = drawconfigf.getInt("BLtrans");
  if(name.find("CRS") != string::npos) transpar = drawconfigf.getInt("CRVtrans");
  if(name.find("Tracker") != string::npos) transpar = drawconfigf.getInt("TRKtrans") ;
  if(name.find("CaloDisk") != string::npos) transpar = drawconfigf.getInt("CALtrans") ;
  if(name.find("CaloWrapper") != string::npos ) transpar = drawconfigf.getInt("CRYtrans") ;

  b1s->SetMainTransparency(transpar);
  b1s->SetMainColor(color);
  b1s->SetEditMainColor(true);
  b1s-> SetEditMainTransparency(true);
  if(!isMother) holder->AddElement(b1s);

  // make 2D projections
  if( crystal1 ){ // add crystals to correct disk
    mngXYCaloDisk0->ImportElements(b1s, XYCaloDisk0GeomScene);
  }if( crystal2 ){
    mngXYCaloDisk1->ImportElements(b1s, XYCaloDisk1GeomScene);
  }

  //show only first tracker plane
  if( val == FrontTracker_gdmltag ){
    mngTrackerXY->ImportElements(b1s, TrackerXYGeomScene); //shows only one plane for 2D view for simplicity
  }

  bool isCRV = name.find("CRS") != string::npos;
  // remove CRV from the YZ view as it blocks tracker/calo view (will have its own view)
  if(!isCRV) { mngRhoZ->ImportElements(b1s, rhoZGeomScene); }
}

int j = 0;
int fp =0;

void MainWindow::showNodesByName(TGeoNode* n, const std::string& str, bool onOff, int _diagLevel, REX::REveTrans& trans,  REX::REveElement* holder, int maxlevel, int level, bool caloshift, bool crystal,  std::vector<double> shift, bool print, bool single, int color) {
  ++level;
  if (level > maxlevel){
    return;
  }
  std::string name(n->GetName());
  j++;
  //if(print) std::cout<<j<<" "<<name<<std::endl;
  bool cry1 = false; // these help us know which disk and to draw crystals
  bool cry2 = false;
  int ndau = n->GetNdaughters();
  if (name.find(str)!= std::string::npos and single){
    REX::REveTrans ctrans;
    ctrans.SetFrom(trans.Array());
    {
      TGeoMatrix     *gm = n->GetMatrix();
      const Double_t *rm = gm->GetRotationMatrix();
      const Double_t *tv = gm->GetTranslation();
      REX::REveTrans t;// NOTE: starts index at "1" whereas the matrices are indexed from 0
      t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
      t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
      t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
      t(1,4) = tv[0] + shift[0]; t(2,4) = tv[1]  + shift[1]; t(3,4) = tv[2] + shift[2];
      //std::cout<<name<<"  "<<tv[0] + shift[0]<<" "<<tv[1]  + shift[1] << " "<< tv[2] + shift[2]<<std::endl;
      ctrans *= t;
    }
    n->ls();
    makeEveGeoShape(n, ctrans, holder, j, cry1, cry2, name, color);
  }
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * pn = n->GetDaughter(i);
    if (name.find(str)!= std::string::npos and i==0 ){ //To make the geometry translucant we only add i==0
      REX::REveTrans ctrans;
      ctrans.SetFrom(trans.Array());
      {
        TGeoMatrix     *gm = n->GetMatrix();
        const Double_t *rm = gm->GetRotationMatrix();
        const Double_t *tv = gm->GetTranslation();
        REX::REveTrans t;// NOTE: starts index at "1" whereas the matrices are indexed from 0
        t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
        t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
        t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
        t(1,4) = tv[0] + shift[0]; t(2,4) = tv[1]  + shift[1]; t(3,4) = tv[2] + shift[2];
       // std::cout<<name<<"  "<<tv[0] + shift[0]<<" "<<tv[1]  + shift[1] << " "<< tv[2] + shift[2]<<std::endl;

        if(name.find("TrackerPlaneEnvelope_00") != string::npos) {
          FrontTracker_gdmltag = j;
        }
        if(caloshift){
          fp++;
          double d = 0;
          if(fp < nCrystals) { d = z_d0 - z_cal; }
          else { d = z_d1 - z_cal; }
          if(crystal) { t(1,4) = tv[0] + shift[0]; t(2,4) = tv[1] + shift[1] ; t(3,4) = tv[2] + d + shift[2]; }
          else { t(1,4) = tv[0] + shift[0]; t(2,4) = tv[1] + shift[1]; t(3,4) = tv[2]+ shift[2]; }
          if (fp < nCrystals and crystal) cry1 = true;
          if (fp >= nCrystals and crystal) cry2 = true;
         }
        ctrans *= t;
      }
      n->ls();
      makeEveGeoShape(n, ctrans, holder, j, cry1, cry2, name, color);
    }
    showNodesByName( pn, str, onOff, _diagLevel, trans,  holder, maxlevel, level, caloshift, crystal, shift, print, single, color);
  }
}

void MainWindow::getOffsets(TGeoNode* n,const std::string& str, REX::REveTrans& trans, int maxlevel, int level, std::vector<std::pair<std::string, std::vector<float>>> & offsets) {
  ++level;
  if (level > maxlevel){
    return;
  }
  std::string name(n->GetName());
  j++;

  REX::REveTrans ctrans;
  ctrans.SetFrom(trans.Array());
  TGeoMatrix     *gm = n->GetMatrix();
  const Double_t *rm = gm->GetRotationMatrix();
  const Double_t *tv = gm->GetTranslation();
  REX::REveTrans t;// NOTE: starts index at "1" whereas the matrices are indexed from 0
  t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
  t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
  t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
  t(1,4) = tv[0]  ; t(2,4) = tv[1]  ; t(3,4) = tv[2] ;

  ctrans *= t;
 // std::cout<<j<<" "<<name<<" located at "<<tv[0]<<" "<< tv[1]<<" "<<tv[2]<<std::endl;
  std::vector<float> pos;
  std::pair<std::string, std::vector<float>> offset;
  pos.push_back(tv[0]);
  pos.push_back(tv[1]);
  pos.push_back(tv[2]);
  offset = make_pair(name, pos);
  offsets.push_back(offset);
  n->ls();
  int ndau = n->GetNdaughters();
  for ( int i=0; i<ndau; ++i ){
    TGeoNode * pn = n->GetDaughter(i);
      REX::REveTrans ctrans;
      ctrans.SetFrom(trans.Array());
      {
        TGeoMatrix     *gm = n->GetMatrix();
        const Double_t *rm = gm->GetRotationMatrix();
        const Double_t *tv = gm->GetTranslation();
        REX::REveTrans t;// NOTE: starts index at "1" whereas the matrices are indexed from 0
        t(1,1) = rm[0]; t(1,2) = rm[1]; t(1,3) = rm[2];
        t(2,1) = rm[3]; t(2,2) = rm[4]; t(2,3) = rm[5];
        t(3,1) = rm[6]; t(3,2) = rm[7]; t(3,3) = rm[8];
        t(1,4) = tv[0]  ; t(2,4) = tv[1]  ; t(3,4) = tv[2] ;

        ctrans *= t;
       // std::cout<<j<<" "<<name<<" located at "<<tv[0]<<" "<< tv[1]<<" "<<tv[2]<<std::endl;
        std::vector<float> pos;
        std::pair<std::string, std::vector<float>> offset;
        pos.push_back(tv[0]);
        pos.push_back(tv[1]);
        pos.push_back(tv[2]);
        offset = make_pair(name, pos);
        offsets.push_back(offset);
        n->ls();

      }
    getOffsets( pn, str, trans,  maxlevel, level, offsets);
  }
}

void MainWindow::GeomDrawerSol(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* beamlineholder, int maxlevel, int level, GeomOptions geomOpt, std::vector<std::pair<std::string, std::vector<float>>>& offsets) {

  double x_sol = 0; double y_sol = 0; double z_sol = 0;
  double x_pt = 0; double y_pt = 0; double z_pt =0;
  for(unsigned int i = 0; i < offsets.size(); i++){
      if(offsets[i].first.find("DS3Vacuum") != string::npos)      {
        x_sol = offsets[i].second[0];
        y_sol  = offsets[i].second[1];
        z_sol = offsets[i].second[2];
      }
      if(offsets[i].first.find("ProductionTarget") != string::npos)      {
        x_pt = offsets[i].second[0];
        y_pt  =offsets[i].second[1];
        z_pt = offsets[i].second[2];
      }
      }
      std::vector<double> shift;
      shift.push_back(-1*x_sol);
      shift.push_back(-1*y_sol);
      shift.push_back(z_sol);
      // set tracker to be central in the frame (shift to 0,0,0)
      if(geomOpt.showPS){
        static std::vector <std::string> substring_ps  {"PSVacuum"};
        for(auto& i: substring_ps){
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
        }
        static std::vector <std::string> substring_pt  {"ProductionTarget"};
        for(auto& i: substring_pt){
          std::vector<double> shiftpt;
          shiftpt.push_back(x_pt);
          shiftpt.push_back(y_pt);
          shiftpt.push_back(z_pt);
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shiftpt, false, true, drawconfigf.getInt("BLColor"));
        }
      }
      if(geomOpt.showTS){
        static std::vector <std::string> substring_ts  {"TS1Vacuum","TS1CryoInsVac","TS2Vacuum","TS2CryoInsVac","TS3Vacuum","TS3CryoInsVac","TS4Vacuum","TS4CryoInsVac","TS5Vacuum","TS5CryoInsVac"};
        for(auto& i: substring_ts){
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
        }
      }
      if(geomOpt.showDS){
        static std::vector <std::string> substring_ds  {"DS*"};
        for(auto& i: substring_ds){
          showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
        }
      }
  }

void MainWindow::GeomDrawerNominal(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* beamlineholder, REX::REveElement* trackerholder, REX::REveElement* caloholder, REX::REveElement* crystalsholder, REX::REveElement* crvholder, REX::REveElement* targetholder, int maxlevel, int level, GeomOptions geomOpt, std::vector<std::pair<std::string, std::vector<float>>>& offsets){
    for(unsigned int i = 0; i < offsets.size(); i++){
      if(offsets[i].first.find("World") != string::npos){
        x_world = offsets[i].second[0];
        y_world = offsets[i].second[1];
        z_world = offsets[i].second[2];
      }
      if(offsets[i].first.find("HallAir") != string::npos)
      {
        x_hall = x_world + offsets[i].second[0];
        y_hall = y_world + offsets[i].second[1];
        z_hall = z_world + offsets[i].second[2];
      }
      if(offsets[i].first.find("DS3Vacuum") != string::npos)      {
        x_ds3 = x_hall + offsets[i].second[0];
        y_ds3 = y_hall + offsets[i].second[1];
        z_ds3 = z_hall + offsets[i].second[2];
        z_crv = -1*offsets[i].second[2];
        x_crv = -1*offsets[i].second[0];
        y_crv = -1*offsets[i].second[1];
      }
      if(offsets[i].first.find("TrackerMother") != string::npos)       {
        x_trk = x_ds3 + offsets[i].second[0];
        y_trk = y_ds3 + offsets[i].second[1];
        z_trk = z_ds3 + offsets[i].second[2];
      }
      if(offsets[i].first.find("CalorimeterMother") != string::npos)       {
        x_cal = x_ds3 + offsets[i].second[0];
        y_cal = y_ds3 + offsets[i].second[1];
        z_cal = z_ds3 + offsets[i].second[2];
      }
      if(offsets[i].first.find("CaloDisk_00") != string::npos){
        x_d0 = x_cal + offsets[i].second[0];
        y_d0 = y_cal + offsets[i].second[1];
        z_d0 = z_cal + offsets[i].second[2];
      }
      if(offsets[i].first.find("CaloDisk_10") != string::npos){
        x_d1 = x_cal + offsets[i].second[0];
        y_d1 = y_cal + offsets[i].second[1];
        z_d1 = z_cal + offsets[i].second[2];
      }
      if(offsets[i].first.find("DS2Vacuum") != string::npos){
        x_ds2 = x_hall + offsets[i].second[0];
        y_ds2 = y_hall + offsets[i].second[1];
        z_ds2 = z_hall + offsets[i].second[2];
      }
      if(offsets[i].first.find("StoppingTargetMother") != string::npos) {
        x_st = x_ds2 + offsets[i].second[0];
        y_st = y_ds2 + offsets[i].second[1];
        z_st = z_ds2 + offsets[i].second[2];
      }
      if(offsets[i].first.find("protonabs1") != string::npos){
        x_ipa = x_st+ offsets[i].second[0];
        y_ipa = y_st + offsets[i].second[1];
        z_ipa = z_st + offsets[i].second[2];
      }
    }
    std::vector<double> shift;

    shift.push_back(0);
    shift.push_back(0);
    shift.push_back(0);
    // set tracker to be central in the frame (shift to 0,0,0)
    if(geomOpt.showTracker){
      static std::vector <std::string> substring_tracker  {"TrackerPlaneEnvelope"};
      for(auto& i: substring_tracker){
        shift.at(0) = 0;
        shift.at(1) = 0;
        shift.at(2) = 0;
        showNodesByName(node,i,kFALSE, 0, trans, trackerholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("TRKColor"));
      }
    }
    // everything else needs to be shifted such that its relative to the tracker center at 0,0,0
    if(geomOpt.showST and !geomOpt.extracted){
      static std::vector <std::string> substrings_stoppingtarget  {"TargetFoil"};
      shift.at(0) = x_st - x_trk;
      shift.at(1) = y_st - y_trk;
      shift.at(2) = z_st - z_trk;
      for(auto& i: substrings_stoppingtarget){
        showNodesByName(node,i,kFALSE, 0, trans, targetholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
      }
    }
    if(geomOpt.showCalo){
        static std::vector <std::string> substrings_disk  {"CaloDisk"};
        for(auto& i: substrings_disk){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk;
          showNodesByName(node,i,kFALSE, 0, trans, caloholder, maxlevel, level, true, false, shift, true, false, drawconfigf.getInt("CALColor") );
        }
      if(geomOpt.showCaloCrystals){
        static std::vector <std::string> substrings_crystals  {"CaloWrapper"};
        for(auto& i: substrings_crystals){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk;
          showNodesByName(node,i,kFALSE, 0, trans, crystalsholder, maxlevel, level, true, true, shift, true, false, drawconfigf.getInt("CALColor"));
        }
      }
    }
    if(geomOpt.showST and !geomOpt.extracted){
      static std::vector <std::string> substrings_pa  {"protonabs1","protonabs3","protonabs4"};
      shift.at(0) = x_ipa - x_trk;
      shift.at(1) = y_ipa - y_trk;
      shift.at(2) = z_ds2 - z_trk;
      for(auto& i: substrings_pa){
        showNodesByName(node,i,kFALSE, 0, trans, beamlineholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("BLColor"));
      }
  }
  if(geomOpt.showCRV and !geomOpt.extracted){
    static std::vector <std::string> substrings_crv {"CRSmother_CRV"};
    shift.at(0) = x_crv;
    shift.at(1) = y_crv;
    shift.at(2) = z_crv - z_trk;
    for(auto& i: substrings_crv){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, false, drawconfigf.getInt("CRVColor"));
    }
  }
}

void MainWindow::GeomDrawerExtracted(TGeoNode* node, REX::REveTrans& trans, REX::REveElement* beamlineholder, REX::REveElement* trackerholder, REX::REveElement* caloholder, REX::REveElement* crystalsholder, REX::REveElement* crvholder, REX::REveElement* targetholder, int maxlevel, int level, GeomOptions geomOpt, std::vector<std::pair<std::string, std::vector<float>>>& offsets){
    double x_crvex = 0; double y_crvex = 0;  double z_crvex = 0; double x_crvt1 = 0; double y_crvt1 = 0;  double z_crvt1 = 0; double x_crvt2 = 0; double y_crvt2 = 0;  double z_crvt2 = 0;
    for(unsigned int i = 0; i < offsets.size(); i++){
      if(offsets[i].first.find("World") != string::npos){
        x_world = offsets[i].second[0];
        y_world = offsets[i].second[1];
        z_world = offsets[i].second[2];
      }
      if(offsets[i].first.find("HallAir") != string::npos)
      {
        x_hall = x_world + offsets[i].second[0];
        y_hall = y_world + offsets[i].second[1];
        z_hall = z_world + offsets[i].second[2];
      }
      if(offsets[i].first.find("garageFakeDS3Vacuum") != string::npos)      {
        x_ds3 = x_hall + offsets[i].second[0];
        y_ds3 = y_hall + offsets[i].second[1];
        z_ds3 = z_hall + offsets[i].second[2];
      }
      if(offsets[i].first.find("TrackerMother") != string::npos)       {
        x_trk = x_ds3 + offsets[i].second[0];
        y_trk = y_ds3 + offsets[i].second[1];
        z_trk = z_ds3 + offsets[i].second[2];
      }
      if(offsets[i].first.find("CalorimeterMother") != string::npos)       {
        x_cal = x_ds3 + offsets[i].second[0];
        y_cal = y_ds3 + offsets[i].second[1];
        z_cal = z_ds3 + offsets[i].second[2];
      }
      if(offsets[i].first.find("CaloDisk_00") != string::npos){
        x_d0 = x_cal + offsets[i].second[0];
        y_d0 = y_cal + offsets[i].second[1];
        z_d0 = z_cal + offsets[i].second[2];
      }
      if(offsets[i].first.find("CaloDisk_10") != string::npos){
        x_d1 = x_cal + offsets[i].second[0];
        y_d1 = y_cal + offsets[i].second[1];
        z_d1 = z_cal + offsets[i].second[2];
      }
      if(offsets[i].first.find("CRSmother_CRV_EX") != string::npos)       {
          x_crvex = offsets[i].second[0]  + x_hall ;
          y_crvex =  offsets[i].second[1] + y_hall;
          z_crvex = offsets[i].second[2]  + z_hall ;
      }
      if(offsets[i].first.find("CRSmother_CRV_T1") != string::npos)       {
          x_crvt1 = offsets[i].second[0] + x_hall ;
          y_crvt1 = offsets[i].second[1] + y_hall ;
          z_crvt1 = offsets[i].second[2] + z_hall ;
      }
      if(offsets[i].first.find("CRSmother_CRV_T2") != string::npos)       {
          x_crvt2 = offsets[i].second[0] + x_hall ;
          y_crvt2 = offsets[i].second[1] + y_hall ;
          z_crvt2 = offsets[i].second[2] + z_hall ;
        }
    }
    std::vector<double> shift;

    shift.push_back(0);
    shift.push_back(0);
    shift.push_back(0);
    // set tracker to be central in the frame (shift to 0,0,0)
    if(geomOpt.showTracker){
      static std::vector <std::string> substring_tracker  {"TrackerPlaneEnvelope"};
      for(auto& i: substring_tracker){
        shift.at(0) = 0;
        shift.at(1) = 0;
        shift.at(2) = 0;
        showNodesByName(node,i,kFALSE, 0, trans, trackerholder, maxlevel, level, false, false, shift, false, true, drawconfigf.getInt("TRKColor"));
      }
    }
    // everything else needs to be shifted such that its relative to the tracker center at 0,0,0

    if(geomOpt.showCalo){
        static std::vector <std::string> substrings_disk  {"CaloDisk"};
        for(auto& i: substrings_disk){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk;
          showNodesByName(node,i,kFALSE, 0, trans, caloholder, maxlevel, level, true, false, shift, false, false, drawconfigf.getInt("CALColor") );
        }
      if(geomOpt.showCaloCrystals){
        static std::vector <std::string> substrings_crystals  {"CaloCrystal","CaloCrystalCsI","CaloCrystalLYSO"};
        for(auto& i: substrings_crystals){
          shift.at(0) = x_cal - x_trk;
          shift.at(1) = y_cal - y_trk;
          shift.at(2) = z_cal - z_trk;
          showNodesByName(node,i,kFALSE, 0, trans, crystalsholder, maxlevel, level, true, true, shift, false, false, drawconfigf.getInt("CALColor"));
        }
      }
    }
  if(geomOpt.showCRV and geomOpt.extracted){

    static std::vector <std::string> substrings_ex {"CRSmotherLayer_CRV_EX"};
    shift.at(0) = x_crvex - x_trk;
    shift.at(1) = y_crvex - y_trk;
    shift.at(2) = z_crvex - z_trk;

    for(auto& i: substrings_ex){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CRVColor"));
    }

    static std::vector <std::string> substrings_t1  {"CRSmotherLayer_CRV_T1"};
    shift.at(0) = x_crvt1 - x_trk;
    shift.at(1) = y_crvt1 - y_trk;
    shift.at(2) = z_crvt1 - z_trk;

    for(auto& i: substrings_t1){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CRVColor"));
    }

    static std::vector <std::string> substrings_t2  {"CRSmotherLayer_CRV_T2"};
    shift.at(0) = x_crvt2 - x_trk;
    shift.at(1) = y_crvt2 - y_trk;
    shift.at(2) = z_crvt2 - z_trk;

    for(auto& i: substrings_t2){
      showNodesByName(node,i,kFALSE, 0, trans, crvholder, maxlevel, level,  false, false, shift, false, true, drawconfigf.getInt("CRVColor"));
    }
  }
}

void MainWindow::projectEvents(REX::REveManager *eveMng)
{
  for (auto &ie : eveMng->GetEventScene()->RefChildren())
  {
    TrackerXYView->SetCameraType(REX::REveViewer::kCameraOrthoXOY);
    XYCaloDisk0View->SetCameraType(REX::REveViewer::kCameraOrthoXOY);
    XYCaloDisk1View->SetCameraType(REX::REveViewer::kCameraOrthoXOY);
    rhoZView->SetCameraType(REX::REveViewer::kCameraOrthoXOY);
    mngTrackerXY->ImportElements(ie, TrackerXYEventScene);

    mngRhoZ  ->ImportElements(ie, rhoZEventScene);

    if(ie->GetName().find("disk0") != string::npos ) mngXYCaloDisk0->ImportElements(ie, XYCaloDisk0EventScene);
    if(ie->GetName().find("disk1") != string::npos ) mngXYCaloDisk1->ImportElements(ie, XYCaloDisk1EventScene);
  }
}

void MainWindow::createProjectionStuff(REX::REveManager *eveMng)
{
  // -------------------Tracker XY View ----------------------------------
  TrackerXYGeomScene  = eveMng->SpawnNewScene("TrackerXY Geometry","TrackerXY");
  TrackerXYEventScene = eveMng->SpawnNewScene("TrackerXY Event Data","TrackerXY");

  mngTrackerXY = new REX::REveProjectionManager(REX::REveProjection::kPT_RPhi);

  TrackerXYView = eveMng->SpawnNewViewer("TrackerXY View", "");
  TrackerXYView->AddScene(TrackerXYGeomScene);
  TrackerXYView->AddScene(TrackerXYEventScene);

  // --------------------Tracker + Calo YZ View ------------------------------

  rhoZGeomScene  = eveMng->SpawnNewScene("ZY Detector Geometry", "ZY");
  rhoZEventScene = eveMng->SpawnNewScene("ZY Event Data","YZ");

  mngRhoZ = new REX::REveProjectionManager(REX::REveProjection::kPT_ZY );

  rhoZView = eveMng->SpawnNewViewer("ZY Detector View", "");
  rhoZView->AddScene(rhoZGeomScene);
  rhoZView->AddScene(rhoZEventScene);

  // ---------------------Calo Disk 1 XY View ----------------------------

  XYCaloDisk0GeomScene  = eveMng->SpawnNewScene("XYCaloDisk0 Geometry", "XYCaloDisk0");
  XYCaloDisk0EventScene = eveMng->SpawnNewScene("XYCaloDisk0 Event Data","XYCaloDisk0");

  mngXYCaloDisk0 = new REX::REveProjectionManager(REX::REveProjection::kPT_RPhi);

  XYCaloDisk0View = eveMng->SpawnNewViewer("XYCaloDisk0 View", "");
  XYCaloDisk0View->AddScene(XYCaloDisk0GeomScene);
  XYCaloDisk0View->AddScene(XYCaloDisk0EventScene);

  // -------------------- Calo Disk 2 XY View ------------------------

  XYCaloDisk1GeomScene  = eveMng->SpawnNewScene("XYCaloDisk1 Geometry", "XYCaloDisk1");
  XYCaloDisk1EventScene = eveMng->SpawnNewScene("XYCaloDisk1 Event Data","XYCaloDisk1");

  mngXYCaloDisk1 = new REX::REveProjectionManager(REX::REveProjection::kPT_RPhi);

  XYCaloDisk1View = eveMng->SpawnNewViewer("XYCaloDisk1 View", "");
  XYCaloDisk1View->AddScene(XYCaloDisk1GeomScene);
  XYCaloDisk1View->AddScene(XYCaloDisk1EventScene);
}


void MainWindow::showEvents(REX::REveManager *eveMng, REX::REveElement* &eventScene, bool firstLoop, bool firstLoopCalo, DataCollections &data, DrawOptions drawOpts, std::vector<int> particleIds, bool strawdisplay, GeomOptions geomOpts, KinKalOptions KKOpts){
  if(!firstLoop){
    eventScene->DestroyElements();
  }
  //...addReco:
  // for start and end of track
  double t1 = 0.;
  double t2 = 1696.;
  std::vector<const KalSeedPtrCollection*> track_list = std::get<1>(data.track_tuple);
  if(drawOpts.addTracks and track_list.size() !=0) {
    pass_data->FillKinKalTrajectory(eveMng, firstLoop, eventScene, data.track_tuple, KKOpts.addKalInter,  KKOpts.addTrkStrawHits, t1, t2); 
  }
  if(drawOpts.addComboHits) {
    std::vector<const ComboHitCollection*> combohit_list = std::get<1>(data.combohit_tuple);
    if(combohit_list.size() !=0 ) pass_data->AddComboHits(eveMng, firstLoop, data.combohit_tuple, eventScene, strawdisplay, drawOpts.addTrkErrBar);
  }

  if(drawOpts.addCRVInfo){
    std::vector<const CrvRecoPulseCollection*> crvpulse_list = std::get<1>(data.crvpulse_tuple);
    if(crvpulse_list.size() !=0) pass_data->AddCRVInfo(eveMng, firstLoop, data.crvpulse_tuple, eventScene, geomOpts.extracted, drawOpts.addCRVBars);
  }

  if(drawOpts.addCRVClusters){
    std::vector<const CrvCoincidenceClusterCollection*> crvcoin_list = std::get<1>(data.crvcoin_tuple);
    if(crvcoin_list.size() !=0) pass_data->AddCRVClusters(eveMng, firstLoop, data.crvcoin_tuple, eventScene, geomOpts.extracted, drawOpts.addCRVBars);
  }

  if(drawOpts.addCaloDigis){
    std::vector<const CaloDigiCollection*> calodigi_list = std::get<1>(data.calodigi_tuple);
    if(calodigi_list.size() !=0 ) pass_data->AddCaloDigis(eveMng, firstLoop, data.calodigi_tuple, eventScene, drawOpts.addCrystalDraw);
  }

  if(drawOpts.addClusters){
    std::vector<const CaloClusterCollection*> calocluster_list = std::get<1>(data.calocluster_tuple);
    if(calocluster_list.size() !=0 ) pass_data->AddCaloClusters(eveMng, firstLoopCalo, data.calocluster_tuple, eventScene, drawOpts.addCrystalDraw, t1, t2);
  }

  std::vector<const HelixSeedCollection*> helix_list = std::get<1>(data.helix_tuple);
  if(drawOpts.addHelices and helix_list.size() !=0) {
    pass_data->AddHelixSeedCollection(eveMng, firstLoop, data.helix_tuple, eventScene);
  }


  if(drawOpts.addCosmicTracks){
    pass_data->AddCosmicTrackFit(eveMng, firstLoop, data.CosmicTrackSeedcol, eventScene);
  }
  if(drawOpts.addTimeClusters){
    std::vector<const TimeClusterCollection*> timecluster_list = std::get<1>(data.timecluster_tuple);
    if(timecluster_list.size() !=0) pass_data->AddTimeClusters(eveMng, firstLoop, data.timecluster_tuple, eventScene);
  }

  //... add MC:
  std::vector<const MCTrajectoryCollection*> mctrack_list = std::get<1>(data.mctrack_tuple);
  if(drawOpts.addMCTrajectories and mctrack_list.size() !=0){
    pass_mc->AddMCTrajectoryCollection(eveMng, firstLoop,  data.mctrack_tuple, eventScene, particleIds, geomOpts.extracted);
  }

  std::vector<const SurfaceStepCollection*> surfstep_list = std::get<1>(data.surfstep_tuple);
  if(drawOpts.addSurfaceSteps and surfstep_list.size() !=0){
    pass_mc->AddSurfaceStepCollection(eveMng, firstLoop,  data.surfstep_tuple, eventScene, particleIds, geomOpts.extracted);
  }

  // ... project these events onto 2D geometry:
  projectEvents(eveMng);
}

void MainWindow::makeGeometryScene(REX::REveManager *eveMng, GeomOptions geomOpt, std::string gdmlname)
{
  std::vector<std::pair<std::string, std::vector<float>>> offsets; // a pair of offsets relative to the world
  TGeoManager *geom = TGeoManager::Import(gdmlname.c_str());
  TGeoVolume* topvol = geom->GetTopVolume();
  gGeoManager->SetTopVolume(topvol);
  gGeoManager->SetTopVisible(kFALSE);
  TGeoNode* topnode = gGeoManager->GetTopNode();
  createProjectionStuff(eveMng);
  std::string name(topnode->GetName());
  {
    // make holders for different parts of geometry
    //auto holder = new REX::REveElement("World");
    auto trackerholder = new REX::REveElement("Tracker");
    auto caloholder = new REX::REveElement("Calorimeter");
    auto crystalsholder = new REX::REveElement("CalorimeterCrystals");
    auto crvholder = new REX::REveElement("CRV");
    auto targetholder = new REX::REveElement("TargetElements");
    auto beamlineholder = new REX::REveElement("BeamlineElements");
    eveMng->GetGlobalScene()->AddElement(trackerholder);
    eveMng->GetGlobalScene()->AddElement(caloholder);
    eveMng->GetGlobalScene()->AddElement(crvholder);
    eveMng->GetGlobalScene()->AddElement(targetholder);
    eveMng->GetGlobalScene()->AddElement(beamlineholder);
    eveMng->GetGlobalScene()->AddElement(crystalsholder);

    REX::REveTrans trans;
    static std::vector <std::string> substrings {"*"};
    for(auto& i: substrings){
      getOffsets(topnode, i, trans, drawconfigf.getInt("maxlevel"),drawconfigf.getInt("level"), offsets);
    }

    if(!geomOpt.extracted) GeomDrawerNominal(topnode, trans, beamlineholder, trackerholder, caloholder, crystalsholder, crvholder, targetholder, drawconfigf.getInt("maxlevel"),drawconfigf.getInt("level"), geomOpt, offsets);
    if(geomOpt.extracted) GeomDrawerExtracted(topnode, trans, beamlineholder, trackerholder, caloholder, crystalsholder, crvholder, targetholder, drawconfigf.getInt("maxlevel"),drawconfigf.getInt("level"), geomOpt, offsets);
    if(geomOpt.showPS or geomOpt.showTS or geomOpt.showDS){
      GeomDrawerSol(topnode, trans, beamlineholder, drawconfigf.getInt("maxlevel"),drawconfigf.getInt("level"), geomOpt, offsets);
    }
  }

  try {
    eveMng->Show();
  } catch (...) {
    std::cout<< "Blocking window" << std::endl;
  }

}

#ifdef __CINT__
ClassImp(MainWindow)
#endif
