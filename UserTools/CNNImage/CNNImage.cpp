#include "CNNImage.h"

CNNImage::CNNImage():Tool(){}

bool CNNImage::Initialise(std::string configfile, DataModel &data){

  if(configfile!="") m_variables.Initialise(configfile); // loading config file
  m_data= &data; //assigning transient data pointer

  //---------------------------------------------------------------
  //-------------------Read in configuration ----------------------
  //---------------------------------------------------------------
  
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

  Log("CNNImage tool: Initialising...",v_message,verbosity);
  
  //Set default configuration values
  if (data_mode != "Normal" && data_mode != "Charge-Weighted" && data_mode != "TimeEvolution") data_mode = "Normal";
  if (save_mode != "Geometric" && save_mode != "PMT-wise") save_mode = "Geometric";
  Log("CNNImage tool: Initialise: Data mode: "+data_mode+" [Normal/Charge-Weighted], Save Mode: "+save_mode+" [PMT-wise/Geometric]",v_message,verbosity);

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

  //---------------------------------------------------------------
  //-------------------Get geometry properties --------------------
  //---------------------------------------------------------------
  
  //Get general geometry properties
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

  Log("CNNImage tool:Initialise: Detector parameters: # PMTs: "+std::to_string(n_tank_pmts)+", radius: "+std::to_string(tank_radius)+", height: "+std::to_string(tank_height)+", center x: "+std::to_string(tank_center_x)+", center y: "+std::to_string(tank_center_y)+", center z: "+std::to_string(tank_center_z),v_message,verbosity);

  //Read in PMT positions
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

  Log("CNNImage tool: Loop over tank detectors finished. Max z = "+std::to_string(max_z)+", min z = "+std::to_string(min_z),v_message,verbosity);

  //---------------------------------------------------------------
  //-------------------Order the PMT positions --------------------
  //---------------------------------------------------------------
  
  std::vector<double> vector_y_top, vector_y_bottom, vector_y_barrel;
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
      vector_y_top.push_back(y);
    }
    else if (z_pmt[detkey] <= min_z+0.001) {
	ConvertPositionTo2D_Bottom(pmt_pos, x, y);
        vector_y_bottom.push_back(y);
    }
    else {
       ConvertPositionTo2D(pmt_pos, x, y);
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
 
  //Print out ordered PMT positions, just for debugging
  if (verbosity > 4) {
    std::cout <<"Sorted 2D position vectors: "<<std::endl;
    for (unsigned int i_x=0;i_x<vec_pmt2D_x.size();i_x++){
      std::cout <<"x vector "<<i_x<<": "<<vec_pmt2D_x.at(i_x)<<std::endl;
    }
    for (unsigned int i_y=0;i_y<vec_pmt2D_y.size();i_y++){
      std::cout <<"y vector "<<i_y<<": "<<vec_pmt2D_y.at(i_y)<<std::endl;
    }
    for (unsigned int i_x=0;i_x<vec_pmt2D_x_Top.size();i_x++){
      std::cout <<"x top vector "<<i_x<<": "<<vec_pmt2D_x_Top.at(i_x)<<std::endl;
    }
    for (unsigned int i_x=0;i_x<vec_pmt2D_x_Bottom.size();i_x++){
      std::cout <<"x bottom vector "<<i_x<<": "<<vec_pmt2D_x_Bottom.at(i_x)<<std::endl;
    }
    for (unsigned int i_y=0; i_y < vector_y_top.size(); i_y++){
      std::cout <<"y (top): "<<vector_y_top.at(i_y)<<std::endl;
    }
    for (unsigned int i_y=0; i_y < vector_y_bottom.size(); i_y++){
      std::cout <<"y (bottom): "<<vector_y_bottom.at(i_y)<<std::endl;
    }
    for (unsigned int i_y=0; i_y < vector_y_barrel.size(); i_y++){
      std::cout <<"y (barrel): "<<vector_y_barrel.at(i_y)<<std::endl;
    }
  }
  
  //---------------------------------------------------------------
  //-------------------Define output files ------------------------
  //---------------------------------------------------------------

  //define root and csv files to save histograms (root-files temporarily, for cross-checks)
  std::string str_root = ".root";
  std::string str_csv = ".csv";
  std::string str_time = "_time";
  std::string str_firsttime = "_firsttime";
  std::string str_charge = "_charge";
  std::string str_Rings = "_Rings";
  std::string str_abs = "_abs";
  std::string rootfile_name = cnn_outpath + str_root;
  std::string csvfile_name = cnn_outpath + str_charge + str_csv;
  std::string csvfile_time_name = cnn_outpath + str_time + str_csv;
  std::string csvfile_firsttime_name = cnn_outpath + str_firsttime + str_csv;
  std::string csvfile_Rings = cnn_outpath + str_Rings + str_csv;
  std::string csvfile_abs = cnn_outpath + str_charge + str_abs + str_csv;
  std::string csvfile_time_abs = cnn_outpath + str_time + str_abs + str_csv;
  std::string csvfile_firsttime_abs = cnn_outpath + str_firsttime + str_abs + str_csv;

  file = new TFile(rootfile_name.c_str(),"RECREATE");
  outfile.open(csvfile_name.c_str());
  outfile_time.open(csvfile_time_name.c_str());
  outfile_firsttime.open(csvfile_firsttime_name.c_str());
  outfile_abs.open(csvfile_abs.c_str());
  outfile_abs_time.open(csvfile_time_abs.c_str());
  outfile_abs_firsttime.open(csvfile_firsttime_abs.c_str());
  outfile_Rings.open(csvfile_Rings.c_str());
  h_trigger_neutron = new TH1F("h_trigger_neutron","Delayed trigger times (neutron)",200,0,70000);
  h_charge_prompt = new TH1F("h_charge_prompt","Prompt charge distribution",200,0,1000);
  h_charge_delayed = new TH1F("h_charge_delayed","Delayed charge distribution",200,0,1000);
  h_time_total = new TH1F("h_time_total","Overall hit time distribution",200,0,2000);
  h_ly_prompt = new TH1F("h_ly_prompt","Light yield prompt event",200,0,100);

  return true;
	
}


bool CNNImage::Execute(){

  if (verbosity >=2) std::cout <<"Executing tool: CNNImage..."<<std::endl;
  
  //---------------------------------------------------------------
  //-------------------Get ANNIEEvent information------------------
  //---------------------------------------------------------------

  int annieeventexists = m_data->Stores.count("ANNIEEvent");
  if(!annieeventexists){ cerr<<"no ANNIEEvent store!"<<endl;}/*return false;*/
  m_data->Stores["ANNIEEvent"]->Get("MCParticles",mcparticles); //needed to retrieve true vertex and direction
  m_data->Stores["ANNIEEvent"]->Get("MCHits", MCHits);
  m_data->Stores["ANNIEEvent"]->Get("EventNumber", evnum);
  m_data->Stores["ANNIEEvent"]->Get("RunNumber",runnumber);
  m_data->Stores["ANNIEEvent"]->Get("SubRunNumber",subrunnumber);
  m_data->Stores["ANNIEEvent"]->Get("EventTime",EventTime);
  m_data->Stores["ANNIEEvent"]->Get("TriggerData",TriggerData); 
  m_data->Stores["ANNIEEvent"]->Get("MCTriggernum",MCTriggernum); 
  m_data->Stores.at("RecoEvent")->Get("NRings",nrings);
  double positron_energy=-1;
  m_data->Stores.at("RecoEvent")->Get("TrueMuonEnergy", positron_energy);
  std::cout <<"positron_energy: "<<positron_energy<<std::endl;

  //General trigger information
  Log("CNNImage tool: Execute: MCTriggernum: "+std::to_string(MCTriggernum),v_message,verbosity);
  for (unsigned int i_trigger = 0; i_trigger < TriggerData->size(); i_trigger++){
    TriggerTime = TriggerData->at(i_trigger).GetTime();
    Log("CNNImage tool: Execute: Trigger #: "+std::to_string(i_trigger)+", TriggerTime: "+std::to_string(TriggerTime.GetNs()),v_message,verbosity);
    if (MCTriggernum == 1) h_trigger_neutron->Fill(TriggerTime.GetNs());
  }

  //Clear variables & containers
  charge.clear();
  time.clear();
  first_time.clear();
  hitpmt_detkeys.clear();
  for (unsigned int i_pmt=0; i_pmt<pmt_detkeys.size();i_pmt++){
    unsigned long detkey = pmt_detkeys[i_pmt];
    charge.emplace(detkey,0.);
    time.emplace(detkey,0.);
    first_time.emplace(detkey,0.);
  }

  std::stringstream ss_hist_time, ss_hist_time_title, ss_hist_charge, ss_hist_charge_title;
  ss_hist_time <<"h_time"<<evnum;
  ss_hist_time_title << "PMT hit times Event "<<evnum;
  ss_hist_charge <<"h_charge"<<evnum;
  ss_hist_charge_title << "Total charge Event "<<evnum;
  TH1F *h_time = new TH1F(ss_hist_time.str().c_str(),ss_hist_time_title.str().c_str(),2000,0,2000);
  TH1F *h_charge = new TH1F(ss_hist_charge.str().c_str(),ss_hist_charge_title.str().c_str(),2000,0,100);

  //---------------------------------------------------------------
  //-------------------Iterate over MCHits ------------------------
  //---------------------------------------------------------------

  int vectsize = MCHits->size();
  Log("CNNImage tool: MCHits size: "+std::to_string(vectsize),v_message,verbosity);
  total_hits_pmts=0;
  double total_charge=0.;
  for(std::pair<unsigned long, std::vector<MCHit>>&& apair : *MCHits){
    unsigned long chankey = apair.first;
    Detector* thistube = geom->ChannelToDetector(chankey);
    unsigned long detkey = thistube->GetDetectorID();
    //if (verbosity > 4) std::cout <<"CNNImage tool: detkey: "<<detkey<<std::endl;
    if (thistube->GetDetectorElement()=="Tank"){
      if (thistube->GetTankLocation()=="OD") continue;
      hitpmt_detkeys.push_back(detkey);
      std::vector<MCHit>& Hits = apair.second;
      int hits_pmt = 0;
      for (MCHit &ahit : Hits){
	if (verbosity > 4) std::cout <<"CNNImage tool: time: "<<ahit.GetTime()<<", charge: "<<ahit.GetCharge()<<std::endl;
        h_time->Fill(ahit.GetTime());
        h_time_total->Fill(ahit.GetTime());
	//TODO: Do we need a time cut here?
	//if (ahit.GetTime()>-10. && ahit.GetTime()<40.){
	  charge[detkey] += ahit.GetCharge();
	  if (data_mode == "Normal") time[detkey] += ahit.GetTime();
	  else if (data_mode == "Charge-Weighted") time[detkey] += (ahit.GetTime()*ahit.GetCharge());
          if (hits_pmt==0) first_time[detkey] = ahit.GetTime();
	  hits_pmt++;
	//}
      }
      h_charge->Fill(charge[detkey]);
      if (data_mode == "Normal" && hits_pmt>0) time[detkey]/=hits_pmt;         //use mean time of all hits on one PMT
      else if (data_mode == "Charge-Weighted" && charge[detkey]>0.) time[detkey] /= charge[detkey];
      total_hits_pmts++;
      total_charge+=charge[detkey];
    }
  }
  Log("CNNImage tool: MCHits loop finished.",v_message,verbosity);
  
  if (MCTriggernum == 0) {
    h_charge_prompt->Fill(total_charge);
    h_ly_prompt->Fill(total_charge/positron_energy);
  }
  else if (MCTriggernum == 1) h_charge_delayed->Fill(total_charge);

  //---------------------------------------------------------------
  //------------- Determine max+min values ------------------------
  //---------------------------------------------------------------

  maximum_pmts = 0;
  max_time_pmts = -999999;
  min_time_pmts = 999999.;
  max_firsttime_pmts = -999999.;
  min_firsttime_pmts = 9999999.;
  total_charge_pmts = 0;
  for (unsigned int i_pmt=0;i_pmt<hitpmt_detkeys.size();i_pmt++){
    unsigned long detkey = hitpmt_detkeys[i_pmt];
    if (charge[detkey]>maximum_pmts) maximum_pmts = charge[detkey];
    total_charge_pmts+=charge[detkey];
    //std::cout<< time[detkey]<< endl;
    if (time[detkey]>max_time_pmts) max_time_pmts = time[detkey];
    if (time[detkey]<min_time_pmts) min_time_pmts = time[detkey];
    if (first_time[detkey]>max_firsttime_pmts) max_firsttime_pmts = first_time[detkey];
    if (first_time[detkey]<min_firsttime_pmts) min_firsttime_pmts = first_time[detkey];
  }
  std::cout<<"Max Time and min time: " << max_time_pmts<<", " << min_time_pmts<<std::endl;
  std::cout <<"Max and min first-time: "<<max_firsttime_pmts<<", "<<min_firsttime_pmts<<std::endl;  

  double global_max_time = max_time_pmts;
  double global_max_charge = maximum_pmts;
  double global_min_charge = 0.;
  double global_min_time = min_time_pmts;

  if (fabs(global_max_time-global_min_time)<0.01) global_max_time = global_min_time+1;
  if (global_max_charge<0.001) global_max_charge=1;  
  if (fabs(max_firsttime_pmts-min_firsttime_pmts)<0.01) max_firsttime_pmts = min_firsttime_pmts+1;  

  //---------------------------------------------------------------
  //-------------- Create CNN images ------------------------------
  //---------------------------------------------------------------

  //define histogram as an intermediate step to the CNN
  std::stringstream ss_cnn, ss_title_cnn, ss_cnn_time, ss_title_cnn_time, ss_cnn_pmtwise, ss_title_cnn_pmtwise, ss_cnn_time_pmtwise, ss_title_cnn_time_pmtwise, ss_cnn_time_first_pmtwise, ss_title_cnn_time_first_pmtwise, ss_cnn_time_first, ss_title_cnn_time_first;
  std::stringstream ss_cnn_abs, ss_title_cnn_abs, ss_cnn_abs_time, ss_title_cnn_abs_time, ss_cnn_abs_pmtwise, ss_title_cnn_abs_pmtwise, ss_cnn_abs_time_pmtwise, ss_title_cnn_abs_time_pmtwise, ss_cnn_abs_time_first_pmtwise, ss_title_cnn_abs_time_first_pmtwise, ss_cnn_abs_time_first, ss_title_cnn_abs_time_first;
  ss_cnn<<"hist_cnn"<<evnum;
  ss_title_cnn<<"EventDisplay (CNN), Event "<<evnum;
  ss_cnn_time<<"hist_cnn_time"<<evnum;
  ss_title_cnn_time<<"EventDisplay Time (CNN), Event "<<evnum;
  ss_cnn_time_first<<"hist_cnn_time_first"<<evnum;
  ss_title_cnn_time_first<<"EventDisplay First HitTime (CNN), Event "<<evnum;
  ss_cnn_abs<<"hist_cnn_abs"<<evnum;
  ss_title_cnn_abs<<"EventDisplay Charge(CNN), Event "<<evnum;
  ss_cnn_abs_time<<"hist_cnn_abs_time"<<evnum;
  ss_title_cnn_abs_time<<"EventDisplay Absolute Time (CNN), Event "<<evnum;
  ss_cnn_abs_time_first<<"hist_cnn_abs_time_first"<<evnum;
  ss_title_cnn_abs_time_first<<"EventDisplay Absolute First HitTime (CNN), Event "<<evnum;
  ss_cnn_pmtwise<<"hist_cnn_pmtwise"<<evnum;
  ss_title_cnn_pmtwise<<"EventDisplay (CNN, pmt wise), Event "<<evnum;
  ss_cnn_time_pmtwise << "hist_cnn_time_pmtwise"<<evnum;
  ss_title_cnn_time_pmtwise <<"EventDisplay Time (CNN, pmt wise), Event "<<evnum;
  ss_cnn_time_first_pmtwise <<"hist_cnn_time_first_pmtwise"<<evnum;
  ss_title_cnn_time_first_pmtwise <<"EventDisplay First Hit Time (CNN, pmt wise), Event "<<evnum;
  ss_cnn_abs_pmtwise<<"hist_cnn_abs_pmtwise"<<evnum;
  ss_title_cnn_abs_pmtwise<<"EventDisplay Charge (CNN, pmt wise), Event "<<evnum;
  ss_cnn_abs_time_pmtwise << "hist_cnn_abs_time_pmtwise"<<evnum;
  ss_title_cnn_abs_time_pmtwise <<"EventDisplay Absolute Time (CNN, pmt wise), Event "<<evnum;
  ss_cnn_abs_time_first_pmtwise <<"hist_cnn_abs_time_first_pmtwise"<<evnum;
  ss_title_cnn_abs_time_first_pmtwise <<"EventDisplay Absolute First Hit Time (CNN, pmt wise), Event "<<evnum;
  TH2F *hist_cnn = new TH2F(ss_cnn.str().c_str(),ss_title_cnn.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_time = new TH2F(ss_cnn_time.str().c_str(),ss_title_cnn_time.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_time_first = new TH2F(ss_cnn_time_first.str().c_str(),ss_title_cnn_time_first.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_abs = new TH2F(ss_cnn_abs.str().c_str(),ss_title_cnn_abs.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_abs_time = new TH2F(ss_cnn_abs_time.str().c_str(),ss_title_cnn_abs_time.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_abs_time_first = new TH2F(ss_cnn_abs_time_first.str().c_str(),ss_title_cnn_abs_time_first.str().c_str(),dimensionX,0.5-TMath::Pi()*size_top_drawing,0.5+TMath::Pi()*size_top_drawing,dimensionY,0.5-(0.45*tank_height/tank_radius+2)*size_top_drawing, 0.5+(0.45*tank_height/tank_radius+2)*size_top_drawing);
  TH2F *hist_cnn_pmtwise = new TH2F(ss_cnn_pmtwise.str().c_str(),ss_title_cnn_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_time_pmtwise = new TH2F(ss_cnn_time_pmtwise.str().c_str(),ss_title_cnn_time_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_time_first_pmtwise = new TH2F(ss_cnn_time_first_pmtwise.str().c_str(),ss_title_cnn_time_first_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_abs_pmtwise = new TH2F(ss_cnn_abs_pmtwise.str().c_str(),ss_title_cnn_abs_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_abs_time_pmtwise = new TH2F(ss_cnn_abs_time_pmtwise.str().c_str(),ss_title_cnn_abs_time_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);
  TH2F *hist_cnn_abs_time_first_pmtwise = new TH2F(ss_cnn_abs_time_first_pmtwise.str().c_str(),ss_title_cnn_abs_time_first_pmtwise.str().c_str(),npmtsX,0,npmtsX,npmtsY,0,npmtsY);

  for (int i_pmt=0;i_pmt<n_tank_pmts;i_pmt++){

    //Convert PMT position to 2D hitmap location
    unsigned long detkey = pmt_detkeys[i_pmt];
    double x,y;
    Position pmt_pos(x_pmt[detkey],y_pmt[detkey],z_pmt[detkey]);
    ConvertPositionTo2D(pmt_pos, x, y);
    
    //Fill geometric 2D-hitmap
    int binx = hist_cnn->GetXaxis()->FindBin(x);
    int biny = hist_cnn->GetYaxis()->FindBin(y);
    Log("CNNImage tool: chankey: "+std::to_string(detkey)+", binx: "+std::to_string(binx)+", biny: "+std::to_string(biny)+", charge fill: "+std::to_string(charge[detkey])+", time fill: "+std::to_string(time[detkey]),vv_debug,verbosity);

    if (maximum_pmts < 0.001) maximum_pmts = 1.;
    double charge_fill = charge[detkey]/global_max_charge;
    hist_cnn->SetBinContent(binx,biny,hist_cnn->GetBinContent(binx,biny)+charge_fill);
    hist_cnn_abs->SetBinContent(binx,biny,hist_cnn_abs->GetBinContent(binx,biny)+charge[detkey]);
    if (fabs(max_time_pmts) < 0.001) max_time_pmts = 1.;
    double time_fill = 0.;
    double time_first_fill = 0.;
    if (charge_fill > 1e-10) {
      time_fill = (time[detkey]-global_min_time)/(global_max_time-global_min_time);
      time_first_fill = (first_time[detkey]-min_firsttime_pmts)/(max_firsttime_pmts-min_firsttime_pmts);
    }
    //For the time files, just accept newest entry as the new overall entry
    hist_cnn_time->SetBinContent(binx,biny,time_fill);
    hist_cnn_time_first->SetBinContent(binx,biny,time_first_fill);
    hist_cnn_abs_time->SetBinContent(binx,biny,time[detkey]);
    hist_cnn_abs_time_first->SetBinContent(binx,biny,first_time[detkey]);

    //Fill the pmt-wise histogram
    if ((z_pmt[detkey]>=max_z || z_pmt[detkey]<=min_z) && !includeTopBottom) continue;       //don't include endcaps in the pmt-wise histogram if specified
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
    hist_cnn_pmtwise->SetBinContent(index_x+1,index_y+1,charge_fill);
    hist_cnn_time_pmtwise->SetBinContent(index_x+1,index_y+1,time_fill);
    hist_cnn_time_first_pmtwise->SetBinContent(index_x+1,index_y+1,time_first_fill);
    hist_cnn_abs_pmtwise->SetBinContent(index_x+1,index_y+1,charge[detkey]);
    hist_cnn_abs_time_pmtwise->SetBinContent(index_x+1,index_y+1,time[detkey]);
    hist_cnn_abs_time_first_pmtwise->SetBinContent(index_x+1,index_y+1,first_time[detkey]);

  }

  //---------------------------------------------------------------
  //-------------- Check event selection --------------------------
  //---------------------------------------------------------------
 
  bool passed_eventselection;
  m_data->Stores["RecoEvent"]->Get("EventCutStatus",passed_eventselection);
  Log("CNNImage: Boolean passed_eventselection: "+std::to_string(passed_eventselection),v_message,verbosity);

  //---------------------------------------------------------------
  //-------------- Write to csv-file ------------------------------
  //---------------------------------------------------------------

  if (passed_eventselection) {
    // safe Ring information
    outfile_Rings << nrings << endl;

    //save root histograms
    hist_cnn->Write();
    hist_cnn_time->Write();
    hist_cnn_time_first->Write();
    hist_cnn_pmtwise->Write();
    hist_cnn_time_pmtwise->Write();
    hist_cnn_time_first_pmtwise->Write();
    hist_cnn_abs->Write();
    hist_cnn_abs_time->Write();
    hist_cnn_abs_time_first->Write();
    hist_cnn_abs_pmtwise->Write();
    hist_cnn_abs_time_pmtwise->Write();
    hist_cnn_abs_time_first_pmtwise->Write();
    h_time->Write();
    h_charge->Write();
    if (save_mode == "Geometric"){
      for (int i_binY=0; i_binY < hist_cnn->GetNbinsY();i_binY++){
	for (int i_binX=0; i_binX < hist_cnn->GetNbinsX();i_binX++){
	  outfile << hist_cnn->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn->GetNbinsX()-1 || i_binY!=hist_cnn->GetNbinsY()-1) outfile<<",";
	  outfile_time << hist_cnn_time->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_time->GetNbinsX()-1 || i_binY!=hist_cnn_time->GetNbinsY()-1) outfile_time<<",";    
	  outfile_firsttime << hist_cnn_time_first->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_time_first->GetNbinsX()-1 || i_binY!=hist_cnn_time_first->GetNbinsY()-1) outfile_firsttime<<",";    
	  outfile_abs << hist_cnn_abs->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_abs->GetNbinsX()-1 || i_binY!=hist_cnn_abs->GetNbinsY()-1) outfile_abs<<",";
	  outfile_abs_time << hist_cnn_abs_time->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_abs_time->GetNbinsX()-1 || i_binY!=hist_cnn_abs_time->GetNbinsY()-1) outfile_abs_time<<",";    
	  outfile_abs_firsttime << hist_cnn_abs_time_first->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_abs_time_first->GetNbinsX()-1 || i_binY!=hist_cnn_abs_time_first->GetNbinsY()-1) outfile_abs_firsttime<<",";    
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
	  outfile_firsttime << hist_cnn_time_first_pmtwise->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_time_first_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_time_first_pmtwise->GetNbinsY()-1) outfile_firsttime<<",";    
	  outfile_abs << hist_cnn_abs_pmtwise->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_abs_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_abs_pmtwise->GetNbinsY()-1) outfile_abs<<",";
	  outfile_abs_time << hist_cnn_abs_time_pmtwise->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_abs_time_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_abs_time_pmtwise->GetNbinsY()-1) outfile_abs_time<<",";    
	  outfile_abs_firsttime << hist_cnn_abs_time_first_pmtwise->GetBinContent(i_binX+1,i_binY+1);
	  if (i_binX != hist_cnn_abs_time_first_pmtwise->GetNbinsX()-1 || i_binY!=hist_cnn_abs_time_first_pmtwise->GetNbinsY()-1) outfile_abs_firsttime<<",";    
	}
      }
    }
    outfile << std::endl;
    outfile_time << std::endl;
    outfile_firsttime << std::endl;
    outfile_abs << std::endl;
    outfile_abs_time << std::endl;
    outfile_abs_firsttime << std::endl;
  }

  return true;
}


bool CNNImage::Finalise(){

  Log("CNNImage tool: Finalising...",v_message,verbosity);
  file->cd();
  h_trigger_neutron->Write();
  h_charge_prompt->Write();
  h_charge_delayed->Write();
  h_time_total->Write();
  h_ly_prompt->Write();
  file->Close();
  outfile.close();
  outfile_time.close();
  outfile_firsttime.close();
  outfile_abs.close();
  outfile_abs_time.close();
  outfile_abs_firsttime.close();
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
	for (int i_slice = num_slices; i_slice >0; i_slice--){
		if (rho > (i_slice-1)*rho_slice){
			y = (51+25+num_slices-i_slice)/double(npmtsY);
			break;
		}
	}

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
