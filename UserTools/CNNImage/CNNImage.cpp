#include "CNNImage.h"

CNNImage::CNNImage():Tool(){}

bool CNNImage::Initialise(std::string configfile, DataModel &data){

  if (verbosity >=2) std::cout <<"Initialising tool: CNNImage..."<<std::endl;

  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  //m_variables.Print();
  m_data= &data; //assigning transient data pointer

  //Read in configuration values
  includeTopBottom = false;
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("DetectorConf",detector_config);
  m_variables.Get("DataMode",data_mode);
  m_variables.Get("SaveMode",save_mode);
  m_variables.Get("DimensionX",dimensionX);
  m_variables.Get("DimensionY",dimensionY);
  m_variables.Get("IncludeTopBottom",includeTopBottom);
  m_variables.Get("OutputFile",cnn_outpath);

  //Set default configuration values
  if (data_mode != "Normal" && data_mode != "Charge-Weighted" && data_mode != "TimeEvolution") data_mode = "Normal";
  if (save_mode != "Geometric" && save_mode != "PMT-wise") save_mode = "Geometric";
  if (verbosity > 0) {
    std::cout <<"data_mode: "<<data_mode<<std::endl;
    std::cout <<"save_mode: "<<save_mode<<std::endl;
  }

  //Dimensions in PMT-wise mode
  npmtsX = 150;    //in every row, we have 150 PMTs
  if (!includeTopBottom) npmtsY = 51;     //we have 51 rows of PMTs (excluding top and bottom PMTs)
  else npmtsY = 101;		//for top & bottom PMTs, include 25 extra rows of PMTs at top and at the bottom

  ifstream phi_file("phi_positions.txt");
  double temp_phi;
  while (!phi_file.eof()){
    phi_file >> temp_phi;
    phi_positions.push_back(temp_phi);
    if (phi_file.eof()) break;
  }
  phi_file.close();

  //Get geometry 
  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
  tank_radius = geom->GetTankRadius();
  tank_height = geom->GetTankHalfheight();
  Position detector_center = geom->GetTankCentre();
  tank_center_x = detector_center.X();
  tank_center_y = detector_center.Y();
  tank_center_z = detector_center.Z();
  n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
  m_data->CStore.Get("channelkey_to_pmtid",channelkey_to_pmtid);
  std::map<std::string,std::map<unsigned long,Detector*> >* Detectors = geom->GetDetectors();

  std::cout <<"n_tank_pmts: "<<n_tank_pmts<<std::endl;
  cout <<"radius: "<<tank_radius<<", height: "<<tank_height<<", center x: "<<tank_center_x<<", center y: "<<tank_center_y<<", center z: "<<tank_center_z<<std::endl;

  //read in PMT positions
  max_z = -1000000.;
  min_z = 1000000.;
  std::cout <<"Tank Detectors loop start"<<std::endl;
  for (std::map<unsigned long,Detector*>::iterator it  = Detectors->at("Tank").begin();
                                                    it != Detectors->at("Tank").end();
                                                  ++it){
    Detector* apmt = it->second;
    unsigned long detkey = it->first;
    pmt_detkeys.push_back(detkey);
    unsigned long chankey = apmt->GetChannels()->begin()->first;
    pmt_chankeys.push_back(chankey);
    Position position_PMT = apmt->GetDetectorPosition();
   // if (verbosity > 2) std::cout <<"detkey: "<<detkey<<", chankey: "<<chankey<<std::endl;
    x_pmt.insert(std::pair<int,double>(detkey,position_PMT.X()-tank_center_x));
    y_pmt.insert(std::pair<int,double>(detkey,position_PMT.Y()-tank_center_y));
    z_pmt.insert(std::pair<int,double>(detkey,position_PMT.Z()-tank_center_z));
    //if (verbosity > 2) std::cout <<"Detector ID: "<<detkey<<", position: ("<<position_PMT.X()<<","<<position_PMT.Y()<<","<<position_PMT.Z()<<")"<<std::endl;
    //if (verbosity > 2) std::cout <<"Rho PMT "<<detkey<<": "<<sqrt(x_pmt.at(detkey)*x_pmt.at(detkey)+z_pmt.at(detkey)*z_pmt.at(detkey))<<std::endl;
    //if (verbosity > 2) std::cout <<"Y PMT: "<<y_pmt.at(detkey)<<std::endl;
    //std::cout <<"tank location: "<<apmt->GetTankLocation()<<", y: "<<y_pmt.at(detkey)<<std::endl;
    if (z_pmt[detkey]>max_z && apmt->GetTankLocation()!="OD") max_z = z_pmt.at(detkey);
    if (z_pmt[detkey]<min_z && apmt->GetTankLocation()!="OD") min_z = z_pmt.at(detkey);

  }

  std::cout <<"Tank Detectors loop end"<<std::endl;
  std::cout <<"max_z = "<<max_z<<std::endl;
  std::cout <<"min_z = "<<min_z<<std::endl;

  std::vector<double> vector_y_top, vector_y_bottom, vector_y_barrel;

  //order the PMT 2D positions 
  for (unsigned int i_pmt = 0; i_pmt < y_pmt.size(); i_pmt++){

    double x,y;
    unsigned long detkey = pmt_detkeys[i_pmt];
    Position pmt_pos(x_pmt[detkey],y_pmt[detkey],z_pmt[detkey]);
    unsigned long chankey = pmt_chankeys[i_pmt];
    Detector *apmt = geom->ChannelToDetector(chankey);
    if (apmt->GetTankLocation()=="OD") continue;  //don't include OD PMTs
    if ((z_pmt[detkey] >= max_z-0.001 || z_pmt[detkey] <= min_z+0.001) && !includeTopBottom) continue;     //don't include top/bottom PMTs if specified
 
    if (z_pmt[detkey] >= max_z-0.001) {
      ConvertPositionTo2D_Top(pmt_pos, x, y);
      //std::cout <<"top: "<<y<<std::endl;
      vector_y_top.push_back(y);
    }
    else if (z_pmt[detkey] <= min_z+0.001) {
	ConvertPositionTo2D_Bottom(pmt_pos, x, y);
        //std::cout <<"bottom: "<<y<<std::endl;
        vector_y_bottom.push_back(y);
    }
    else {
       ConvertPositionTo2D(pmt_pos, x, y);
       //std::cout <<"barrel: "<<y<<std::endl;
       vector_y_barrel.push_back(y);
    }
    x = (round(1000*x)/1000.);
    y = (round(1000*y)/1000.);
    if (z_pmt[detkey] >= max_z-0.001) vec_pmt2D_x_Top.push_back(x);
    else if (z_pmt[detkey] <= min_z+0.001) vec_pmt2D_x_Bottom.push_back(x);
    else vec_pmt2D_x.push_back(x);
    vec_pmt2D_y.push_back(y);
    //if (verbosity > 2) std::cout <<detkey<<"  "<<x<<"  "<<y<<std::endl;
  }

  if (verbosity > 2) std::cout <<"vec_pmt2D_* size: "<<vec_pmt2D_x.size()<<std::endl;
  std::sort(vec_pmt2D_x.begin(),vec_pmt2D_x.end());
  std::sort(vec_pmt2D_y.begin(),vec_pmt2D_y.end());
  std::sort(vec_pmt2D_x_Top.begin(),vec_pmt2D_x_Top.end());
  std::sort(vec_pmt2D_x_Bottom.begin(),vec_pmt2D_x_Bottom.end());
  std::sort(vector_y_top.begin(),vector_y_top.end());
  std::sort(vector_y_bottom.begin(),vector_y_bottom.end());
  std::sort(vector_y_barrel.begin(),vector_y_barrel.end());
  vec_pmt2D_x.erase(std::unique(vec_pmt2D_x.begin(),vec_pmt2D_x.end()),vec_pmt2D_x.end());
  vec_pmt2D_y.erase(std::unique(vec_pmt2D_y.begin(),vec_pmt2D_y.end()),vec_pmt2D_y.end());
  vec_pmt2D_x_Top.erase(std::unique(vec_pmt2D_x_Top.begin(),vec_pmt2D_x_Top.end()),vec_pmt2D_x_Top.end());
  vec_pmt2D_x_Bottom.erase(std::unique(vec_pmt2D_x_Bottom.begin(),vec_pmt2D_x_Bottom.end()),vec_pmt2D_x_Bottom.end());
  vector_y_top.erase(std::unique(vector_y_top.begin(),vector_y_top.end()),vector_y_top.end());
  vector_y_bottom.erase(std::unique(vector_y_bottom.begin(),vector_y_bottom.end()),vector_y_bottom.end());
  vector_y_barrel.erase(std::unique(vector_y_barrel.begin(),vector_y_barrel.end()),vector_y_barrel.end());
 
  /*if (verbosity > 1) std::cout <<"Sorted 2D position vectors: "<<std::endl;
  for (unsigned int i_x=0;i_x<vec_pmt2D_x.size();i_x++){
    if (verbosity > 1) std::cout <<"x vector "<<i_x<<": "<<vec_pmt2D_x.at(i_x)<<std::endl;
  }
  for (unsigned int i_y=0;i_y<vec_pmt2D_y.size();i_y++){
    if (verbosity > 1) std::cout <<"y vector "<<i_y<<": "<<vec_pmt2D_y.at(i_y)<<std::endl;
  }
  for (unsigned int i_x=0;i_x<vec_pmt2D_x_Top.size();i_x++){
    if (verbosity > 1) std::cout <<"x top vector "<<i_x<<": "<<vec_pmt2D_x_Top.at(i_x)<<std::endl;
  }
  for (unsigned int i_x=0;i_x<vec_pmt2D_x_Bottom.size();i_x++){
    if (verbosity > 1) std::cout <<"x bottom vector "<<i_x<<": "<<vec_pmt2D_x_Bottom.at(i_x)<<std::endl;
  }
  for (unsigned int i_y=0; i_y < vector_y_top.size(); i_y++){
    std::cout <<"y (top): "<<vector_y_top.at(i_y)<<std::endl;
  }
  for (unsigned int i_y=0; i_y < vector_y_bottom.size(); i_y++){
    std::cout <<"y (bottom): "<<vector_y_bottom.at(i_y)<<std::endl;
  }
  for (unsigned int i_y=0; i_y < vector_y_barrel.size(); i_y++){
    std::cout <<"y (barrel): "<<vector_y_barrel.at(i_y)<<std::endl;
  }*/

  //define root and csv files to save histograms (root-files temporarily, for cross-checks)

  std::string str_root = ".root";
  std::string str_csv = ".csv";
  std::string str_time = "_time";
  std::string str_charge = "_charge";
  std::string str_Rings = "_Rings";
  std::string rootfile_name = cnn_outpath + str_root;
  std::string csvfile_name = cnn_outpath + str_charge + str_csv;
  std::string csvfile_time_name = cnn_outpath + str_time + str_csv;
  std::string csvfile_Rings = cnn_outpath + str_Rings + str_csv;

  file = new TFile(rootfile_name.c_str(),"RECREATE");
  outfile.open(csvfile_name.c_str());
  outfile_time.open(csvfile_time_name.c_str());
  outfile_Rings.open(csvfile_Rings.c_str());

  return true;
	
}


bool CNNImage::Execute(){

  if (verbosity >=2) std::cout <<"Executing tool: CNNImage..."<<std::endl;

  //get ANNIEEvent store information
  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if(!annieeventexists){ cerr<<"no ANNIEEvent store!"<<endl;}/*return false;*/

  m_data->Stores["ANNIEEvent"]->Get("MCParticles",mcparticles); //needed to retrieve true vertex and direction
  m_data->Stores["ANNIEEvent"]->Get("MCHits", MCHits);
  m_data->Stores["ANNIEEvent"]->Get("EventNumber", evnum);
  m_data->Stores["ANNIEEvent"]->Get("RunNumber",runnumber);
  m_data->Stores["ANNIEEvent"]->Get("SubRunNumber",subrunnumber);
  m_data->Stores["ANNIEEvent"]->Get("EventTime",EventTime);
  m_data->Stores.at("RecoEvent")->Get("NRings",nrings);

  //clear variables & containers
  charge.clear();
  time.clear();
  hitpmt_detkeys.clear();
  for (unsigned int i_pmt=0; i_pmt<pmt_detkeys.size();i_pmt++){
    unsigned long detkey = pmt_detkeys[i_pmt];
    charge.emplace(detkey,0.);
    time.emplace(detkey,0.);
  }


  TH1F *h_time = new TH1F("h_time","Super-K PMT times",2000,0,2000);
  TH1F *h_charge = new TH1F("h_charge","Super-K PMT charge",2000,0,1000);

  //make basic selection cuts to only look at clear event signatures

  bool bool_primary=false;
  bool bool_geometry=false;
  bool bool_nhits=false;

  if (verbosity > 1) std::cout <<"Loop through MCParticles..."<<std::endl;
  for(unsigned int particlei=0; particlei<mcparticles->size(); particlei++){
    MCParticle aparticle = mcparticles->at(particlei);
    if (verbosity > 2) std::cout <<"particle "<<particlei<<std::endl;
    if (verbosity > 2) std::cout <<"Parent ID: "<<aparticle.GetParentPdg()<<std::endl;
    if (verbosity > 2) std::cout <<"PDG code: "<<aparticle.GetPdgCode()<<std::endl;
    if (verbosity > 2) std::cout <<"Flag: "<<aparticle.GetFlag()<<std::endl;
    if (aparticle.GetParentPdg() !=0 ) continue;
    if (aparticle.GetFlag() !=0 ) continue;
//    if (!(aparticle.GetPdgCode() == 11 || aparticle.GetPdgCode() == 13)) continue;    //primary particle for Cherenkov tracking should be muon or electron
//    else {
      truevtx = aparticle.GetStartVertex();
      truevtx_x = truevtx.X()-tank_center_x;
      truevtx_y = truevtx.Y()-tank_center_y;
      truevtx_z = truevtx.Z()-tank_center_z;
      double distInnerStr_Hor = tank_radius - sqrt(pow(truevtx_x,2)+pow(truevtx_z,2));
      double distInnerStr_Vert1 = max_z - truevtx_z;
      double distInnerStr_Vert2 = truevtx_z - min_z;
      double distInnerStr_Vert;
      if (distInnerStr_Vert1 > 0 && distInnerStr_Vert2 > 0) {
	if (distInnerStr_Vert1>distInnerStr_Vert2) distInnerStr_Vert=distInnerStr_Vert2;
	else distInnerStr_Vert=distInnerStr_Vert1;
      } else if (distInnerStr_Vert1 <=0) distInnerStr_Vert=distInnerStr_Vert2;
      else distInnerStr_Vert=distInnerStr_Vert1;
      bool_geometry = (distInnerStr_Vert>0.2 && distInnerStr_Hor>0.2);
      bool_primary = true;
//    }
  }

  //---------------------------------------------------------------
  //-------------------Iterate over MCHits ------------------------
  //---------------------------------------------------------------

  int vectsize = MCHits->size();
  if (verbosity > 1) std::cout <<"Tool CNNImage: MCHits size: "<<vectsize<<std::endl;
  total_hits_pmts=0;
  for(std::pair<unsigned long, std::vector<MCHit>>&& apair : *MCHits){
    unsigned long chankey = apair.first;
    if (verbosity > 1) std::cout <<"chankey: "<<chankey;
    Detector* thistube = geom->ChannelToDetector(chankey);
    unsigned long detkey = thistube->GetDetectorID();
    //if (verbosity > 3) std::cout <<"detkey: "<<detkey<<std::endl;
    if (thistube->GetDetectorElement()=="Tank"){
      if (thistube->GetTankLocation()=="OD") continue;
      hitpmt_detkeys.push_back(detkey);
      std::vector<MCHit>& Hits = apair.second;
      int hits_pmt = 0;
      for (MCHit &ahit : Hits){
	std::cout <<", time: "<<ahit.GetTime()<<", charge: "<<ahit.GetCharge()<<std::endl;
        h_time->Fill(ahit.GetTime());
	//if (ahit.GetTime()>-10. && ahit.GetTime()<40.){
	  charge[detkey] += ahit.GetCharge();
	  if (data_mode == "Normal") time[detkey] += ahit.GetTime();
	  else if (data_mode == "Charge-Weighted") time[detkey] += (ahit.GetTime()*ahit.GetCharge());
	  hits_pmt++;
	//}
      }
      h_charge->Fill(charge[detkey]);
      if (data_mode == "Normal" && hits_pmt>0) time[detkey]/=hits_pmt;         //use mean time of all hits on one PMT
      else if (data_mode == "Charge-Weighted" && charge[detkey]>0.) time[detkey] /= charge[detkey];
      total_hits_pmts++;
    }
  }

  std::cout <<"done with MCHits loop"<<std::endl;
  if (total_hits_pmts>=10) bool_nhits=true;


  //---------------------------------------------------------------
  //------------- Determine max+min values ------------------------
  //---------------------------------------------------------------

  maximum_pmts = 0;
  max_time_pmts = -999999;
  min_time_pmts = 999999.;
  total_charge_pmts = 0;
  for (unsigned int i_pmt=0;i_pmt<hitpmt_detkeys.size();i_pmt++){
    unsigned long detkey = hitpmt_detkeys[i_pmt];
    if (charge[detkey]>maximum_pmts) maximum_pmts = charge[detkey];
    total_charge_pmts+=charge[detkey];
    std::cout<< time[detkey]<< endl;
    if (time[detkey]>max_time_pmts) max_time_pmts = time[detkey];
    if (time[detkey]<min_time_pmts) min_time_pmts = time[detkey];
  }
  std::cout<<"Max Time and min time: " << max_time_pmts<<", " << min_time_pmts<< endl;
  
  double global_max_time = max_time_pmts;
  double global_max_charge = maximum_pmts;
  double global_min_charge = 0.;
  double global_min_time = min_time_pmts;

  if (fabs(global_max_time-global_min_time)<0.01) global_max_time = global_min_time+1;
  if (global_max_charge<0.001) global_max_charge=1;  
  
  //---------------------------------------------------------------
  //-------------- Create CNN images ------------------------------
  //---------------------------------------------------------------

  //define histogram as an intermediate step to the CNN
  std::stringstream ss_cnn, ss_title_cnn, ss_cnn_time, ss_title_cnn_time, ss_cnn_pmtwise, ss_title_cnn_pmtwise, ss_cnn_time_pmtwise, ss_title_cnn_time_pmtwise;
  ss_cnn<<"hist_cnn"<<evnum;
  ss_title_cnn<<"EventDisplay (CNN), Event "<<evnum;
  ss_cnn_time<<"hist_cnn_time"<<evnum;
  ss_title_cnn_time<<"EventDisplay Time (CNN), Event "<<evnum;
  ss_cnn_pmtwise<<"hist_cnn_pmtwise"<<evnum;
  ss_title_cnn_pmtwise<<"EventDisplay (CNN, pmt wise), Event "<<evnum;
  ss_cnn_time_pmtwise << "hist_cnn_time_pmtwise"<<evnum;
  ss_title_cnn_time_pmtwise <<"EventDisplay Time (CNN, pmt wise), Event "<<evnum;
  TH2F *hist_cnn = new TH2F(ss_cnn.str().c_str(),ss_title_cnn.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_time = new TH2F(ss_cnn_time.str().c_str(),ss_title_cnn_time.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_pmtwise = new TH2F(ss_cnn_pmtwise.str().c_str(),ss_title_cnn_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_time_pmtwise = new TH2F(ss_cnn_time_pmtwise.str().c_str(),ss_title_cnn_time_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);

  for (int i_pmt=0;i_pmt<n_tank_pmts;i_pmt++){
    unsigned long detkey = pmt_detkeys[i_pmt];

    double x,y;
    Position pmt_pos(x_pmt[detkey],y_pmt[detkey],z_pmt[detkey]);
    ConvertPositionTo2D(pmt_pos, x, y);
    
    int binx = hist_cnn->GetXaxis()->FindBin(x);
    int biny = hist_cnn->GetYaxis()->FindBin(y);
    if (verbosity > 2) std::cout <<"CNNImage tool: chankey: "<<detkey<<", binx: "<<binx<<", biny: "<<biny<<", charge fill: "<<charge[detkey]<<", time fill: "<<time[detkey]<<std::endl;

    if (maximum_pmts < 0.001) maximum_pmts = 1.;
    if (fabs(max_time_pmts) < 0.001) max_time_pmts = 1.;
   
    // double charge_fill = charge[detkey]/maximum_pmts;
    double charge_fill = charge[detkey]/global_max_charge;

    hist_cnn->SetBinContent(binx,biny,hist_cnn->GetBinContent(binx,biny)+charge_fill);
    //    double time_fill = (time[detkey]-min_time_pmts)/(max_time_pmts-min_time_pmts);
    double time_fill = 0.;
    if (charge_fill > 1e-10) time_fill = (time[detkey]-global_min_time)/(global_max_time-global_min_time);

    hist_cnn_time->SetBinContent(binx,biny,hist_cnn_time->GetBinContent(binx,biny)+time_fill);

  //  std::cout<<"Time Detector key: " << time[detkey] << ", Time fill: " <<   time_fill <<endl;

    //fill the pmt-wise histogram
    if ((z_pmt[detkey]>=max_z || z_pmt[detkey]<=min_z) && !includeTopBottom) continue;       //don't include endcaps in the pmt-wise histogram for now
    if (z_pmt[detkey]>=max_z) ConvertPositionTo2D_Top(pmt_pos,x,y);
    if (z_pmt[detkey]<=min_z) ConvertPositionTo2D_Bottom(pmt_pos,x,y);
    double xCorr, yCorr;
    xCorr = (round(1000*x)/1000.);
    yCorr = (round(1000*y)/1000.);
    std::vector<double>::iterator it_x, it_y;
    if (z_pmt[detkey]>=max_z){
      it_x = std::find(vec_pmt2D_x_Top.begin(),vec_pmt2D_x_Top.end(),xCorr);
    }
    else if (z_pmt[detkey]<=min_z){
      it_x = std::find(vec_pmt2D_x_Bottom.begin(),vec_pmt2D_x_Bottom.end(),xCorr);
    }
    else {
      it_x = std::find(vec_pmt2D_x.begin(),vec_pmt2D_x.end(),xCorr);
    }
    it_y = std::find(vec_pmt2D_y.begin(),vec_pmt2D_y.end(),yCorr);
    int index_x, index_y;
    if (z_pmt[detkey]>=max_z) index_x = std::distance(vec_pmt2D_x_Top.begin(),it_x);
    else if (z_pmt[detkey]<=min_z) index_x = std::distance(vec_pmt2D_x_Bottom.begin(),it_x);
    else index_x = std::distance(vec_pmt2D_x.begin(),it_x);
    index_y = std::distance(vec_pmt2D_y.begin(),it_y);
    //std::cout <<"index_x: "<<index_x<<", index_y: "<<index_y<<std::endl;
    hist_cnn_pmtwise->SetBinContent(index_x+1,index_y+1,charge_fill);
    hist_cnn_time_pmtwise->SetBinContent(index_x+1,index_y+1,time_fill);

  }

  //save information from histogram to csv file
  //(1 line corresponds to 1 event, histogram entries flattened out to a 1D array)
 
  bool passed_eventselection;
  m_data->Stores["RecoEvent"]->Get("EventCutStatus",passed_eventselection);
  std::cout <<"passed_eventselection: "<<passed_eventselection<<std::endl;



  //if (bool_primary && bool_geometry && bool_nhits) {
  if (passed_eventselection) {
    // safe Rings and MRD information
    outfile_Rings << nrings << endl;

    //save root histograms
    hist_cnn->Write();
    hist_cnn_time->Write();
    hist_cnn_pmtwise->Write();
    hist_cnn_time_pmtwise->Write();
    h_time->Write();
    h_charge->Write();
    if (save_mode == "Geometric"){
      for (int i_binY=0; i_binY < hist_cnn->GetNbinsY();i_binY++){
	for (int i_binX=0; i_binX < hist_cnn->GetNbinsX();i_binX++){
	  outfile << hist_cnn->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn->GetNbinsX()-1 || i_binY!=hist_cnn->GetNbinsY()-1) outfile<<",";
	  outfile_time << hist_cnn_time->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_time->GetNbinsX()-1 || i_binY!=hist_cnn_time->GetNbinsY()-1) outfile_time<<",";    
	}
      }
    } else if (save_mode == "PMT-wise"){
      //std::cout <<"Entering PMT-wise save mode"<<std::endl;
      for (int i_binY=0; i_binY < hist_cnn_pmtwise->GetNbinsY();i_binY++){
	for (int i_binX=0; i_binX < hist_cnn_pmtwise->GetNbinsX();i_binX++){
	  outfile << hist_cnn_pmtwise->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_pmtwise->GetNbinsY()-1) outfile<<",";
	  outfile_time << hist_cnn_time_pmtwise->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_time_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_time_pmtwise->GetNbinsY()-1) outfile_time<<",";    
	}
      }
    }
    outfile << std::endl;
    outfile_time << std::endl;
  }

  return true;
}


bool CNNImage::Finalise(){

  if (verbosity >=2 ) std::cout <<"Finalising tool: CNNImage..."<<std::endl;
  file->Close();
  outfile.close();
  outfile_time.close();
  outfile_Rings.close();
  return true;
}

void CNNImage::ConvertPositionTo2D(Position xyz_pos, double &x, double &y){

    if (fabs(xyz_pos.Z()-max_z)<0.01){
      //top PMTs
      x=0.5-size_top_drawing*xyz_pos.X()/tank_radius;
      y=0.5+((0.45*tank_height)/tank_radius+1)*size_top_drawing-size_top_drawing*xyz_pos.Y()/tank_radius;
    } else if (fabs(xyz_pos.Z()-min_z)<0.01){
      //bottom PMTs
      x=0.5-size_top_drawing*xyz_pos.X()/tank_radius;
      y=0.5-(0.45*tank_height/tank_radius+1)*size_top_drawing+size_top_drawing*xyz_pos.Y()/tank_radius;
    } else {
      //barrel PMTs
      double phi=0.;
      if (xyz_pos.Y()>0 && xyz_pos.X()>0) phi = atan(xyz_pos.X()/xyz_pos.Y())+TMath::Pi()/2;
      else if (xyz_pos.Y()>0 && xyz_pos.X()<0) phi = atan(xyz_pos.Y()/-xyz_pos.X());
      else if (xyz_pos.Y()<0 && xyz_pos.X()<0) phi = 3*TMath::Pi()/2+atan(xyz_pos.X()/xyz_pos.Y());
      else if (xyz_pos.Y()<0 && xyz_pos.X()>0) phi = TMath::Pi()+atan(-xyz_pos.Y()/xyz_pos.X());
      else if (fabs(xyz_pos.Y())<0.0001){
	if (xyz_pos.X()>0) phi = TMath::Pi();
	else if (xyz_pos.X()<0) phi = 2*TMath::Pi();
      }
      else if (fabs(xyz_pos.X())<0.0001){
	if (xyz_pos.Y()>0) phi = 0.5*TMath::Pi();
	else if (xyz_pos.Y()<0) phi = 3*TMath::Pi()/2;
      }
      else phi = 0.;
      if (phi>2*TMath::Pi()) phi-=(2*TMath::Pi());
      phi-=TMath::Pi();
      if (phi < - TMath::Pi()) phi = -TMath::Pi();
      if (phi<-TMath::Pi() || phi>TMath::Pi())  std::cout <<"Drawing Event: Phi out of bounds! "<<", x= "<<xyz_pos.X()<<", y="<<xyz_pos.Y()<<", z="<<xyz_pos.Z()<<std::endl;
      x=0.5+phi*size_top_drawing;
      y=0.5+xyz_pos.Z()/tank_height*tank_height/tank_radius*size_top_drawing;
      }
}

void CNNImage::ConvertPositionTo2D_Top(Position xyz_pos, double &x, double &y){

	double rho = sqrt(xyz_pos.X()*xyz_pos.X()+xyz_pos.Y()*xyz_pos.Y());
        double rho_slice = 0.6666666;
	int num_slices = 25;
	//std::cout <<"rho (top): "<<rho<<std::endl;
	for (int i_slice = num_slices; i_slice >0; i_slice--){
		if (rho > (i_slice-1)*rho_slice){
			y = (51+25+num_slices-i_slice)/double(npmtsY);
			break;
		}
	}
        //std::cout <<"y: "<<y<<std::endl;

	double phi = 0.;
	if (xyz_pos.Y()>0 && xyz_pos.X()>0) phi = atan(xyz_pos.X()/xyz_pos.Y())+TMath::Pi()/2;
	else if (xyz_pos.Y()>0 && xyz_pos.X()<0) phi = atan(xyz_pos.Y()/-xyz_pos.X());
	else if (xyz_pos.Y()<0 && xyz_pos.X()<0) phi = 3*TMath::Pi()/2+atan(xyz_pos.X()/xyz_pos.Y());
	else if (xyz_pos.Y()<0 && xyz_pos.X()>0) phi = TMath::Pi()+atan(-xyz_pos.Y()/xyz_pos.X());
	else if (fabs(xyz_pos.Y())<0.0001){
	  if (xyz_pos.X()>0) phi = TMath::Pi();
	  else if (xyz_pos.X()<0) phi = 2*TMath::Pi();
	}
	else if (fabs(xyz_pos.X())<0.0001){
	  if (xyz_pos.Y()>0) phi = 0.5*TMath::Pi();
	  else if (xyz_pos.Y()<0) phi = 3*TMath::Pi()/2;
	}
	else phi = 0.;
	if (phi>2*TMath::Pi()) phi-=(2*TMath::Pi());
	phi-=TMath::Pi();
	if (phi < - TMath::Pi()) phi = -TMath::Pi();
	if (phi<-TMath::Pi() || phi>TMath::Pi())  std::cout <<"Drawing Event: Phi out of bounds! "<<", x= "<<xyz_pos.X()<<", y="<<xyz_pos.Y()<<", z="<<xyz_pos.Z()<<std::endl;
	x=0.5+phi*size_top_drawing; 
        double diff = 100000.;
        double x_min=0.;
	for (int i_phi = 0; i_phi < (int) phi_positions.size(); i_phi++){
          double temp_diff = fabs(phi_positions.at(i_phi)-x);
          if (temp_diff < diff) {x_min = phi_positions.at(i_phi); diff = temp_diff;}
        }
        x = x_min;

}

void CNNImage::ConvertPositionTo2D_Bottom(Position xyz_pos, double &x, double &y){

        double rho_slice = 0.6666666;
	int num_slices = 25;
	double rho = sqrt(xyz_pos.X()*xyz_pos.X()+xyz_pos.Y()*xyz_pos.Y());
	for (int i_slice = num_slices; i_slice >0; i_slice--){
		if (rho > (i_slice-1)*rho_slice){
			y = (25-(num_slices-i_slice))/double(npmtsY);
			break;
		}
	}
 
	double phi=0.;  
	if (xyz_pos.Y()>0 && xyz_pos.X()>0) phi = atan(xyz_pos.X()/xyz_pos.Y())+TMath::Pi()/2;
	else if (xyz_pos.Y()>0 && xyz_pos.X()<0) phi = atan(xyz_pos.Y()/-xyz_pos.X());
	else if (xyz_pos.Y()<0 && xyz_pos.X()<0) phi = 3*TMath::Pi()/2+atan(xyz_pos.X()/xyz_pos.Y());
	else if (xyz_pos.Y()<0 && xyz_pos.X()>0) phi = TMath::Pi()+atan(-xyz_pos.Y()/xyz_pos.X());
	else if (fabs(xyz_pos.Y())<0.0001){
	  if (xyz_pos.X()>0) phi = TMath::Pi();
	  else if (xyz_pos.X()<0) phi = 2*TMath::Pi();
	}
	else if (fabs(xyz_pos.X())<0.0001){
	  if (xyz_pos.Y()>0) phi = 0.5*TMath::Pi();
	  else if (xyz_pos.Y()<0) phi = 3*TMath::Pi()/2;
	}
	else phi = 0.;
	if (phi>2*TMath::Pi()) phi-=(2*TMath::Pi());
	phi-=TMath::Pi();
	if (phi < - TMath::Pi()) phi = -TMath::Pi();
	if (phi<-TMath::Pi() || phi>TMath::Pi())  std::cout <<"Drawing Event: Phi out of bounds! "<<", x= "<<xyz_pos.X()<<", y="<<xyz_pos.Y()<<", z="<<xyz_pos.Z()<<std::endl;
	x=0.5+phi*size_top_drawing;

        double diff = 100000.;
        double x_min=0.;
	for (int i_phi = 0; i_phi < (int) phi_positions.size(); i_phi++){
          double temp_diff = fabs(phi_positions.at(i_phi)-x);
          if (temp_diff < diff) {x_min = phi_positions.at(i_phi); diff = temp_diff;}
        }
        x = x_min;

}
