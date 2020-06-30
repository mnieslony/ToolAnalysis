#include "CreatePIDPDF.h"

CreatePIDPDF::CreatePIDPDF():Tool(){}


bool CreatePIDPDF::Initialise(std::string configfile, DataModel &data){

	//---------------------------------------------------------------
	//----------------- Useful header -------------------------------
	//---------------------------------------------------------------

	if(configfile!="")  m_variables.Initialise(configfile); 
	//m_variables.Print();
	m_data= &data; 
	m_variables.Get("verbosity",verbosity);
	m_variables.Get("OutputFile_PDF",filename_pdf);
	m_variables.Get("OutputFile_CSV",filename_csv);
	m_variables.Get("OutputFile_CSVAverage",filename_csv_average);
	m_variables.Get("OutputFile_CSVData",filename_csv_data);
	m_variables.Get("OutputFile_ROOT",filename_root);
	m_variables.Get("OutputFile_ROOTOverview",filename_root_overview);
	m_variables.Get("UseParametric",parametric);

	std::cout <<"filename pdf: "<<filename_pdf<<std::endl;
	std::cout <<"filename root: "<<filename_root<<std::endl;
	std::cout <<"filename root (overview): "<<filename_root_overview<<std::endl;
	std::cout <<"filename csv: "<<filename_csv<<std::endl;
	std::cout <<"filename csv (average): "<<filename_csv_average<<std::endl;

	m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);
	Position detector_center = geom->GetTankCentre();
	tank_center_x = detector_center.X();
	tank_center_y = detector_center.Y();
	tank_center_z = detector_center.Z();
	tank_radius = geom->GetTankRadius();
	detector_h = geom->GetTankHalfheight();
	detector_h/=2;
	std::cout <<"Detector radius: "<<tank_radius<<", detector height: "<<detector_h<<std::endl;

	Position position_PMT;
	n_tank_pmts = geom->GetNumDetectorsInSet("Tank");
	n_lappds = geom->GetNumDetectorsInSet("LAPPD");
	n_mrd_pmts = geom->GetNumDetectorsInSet("MRD");
	n_veto_pmts = geom->GetNumDetectorsInSet("Veto");
	m_data->CStore.Get("detectorkey_to_lappdid",detectorkey_to_lappdid);
	m_data->CStore.Get("channelkey_to_pmtid",channelkey_to_pmtid);
	m_data->CStore.Get("channelkey_to_mrdpmtid",channelkey_to_mrdpmtid);
	m_data->CStore.Get("channelkey_to_faccpmtid",channelkey_to_faccpmtid);
	m_data->CStore.Get("pmt_tubeid_to_channelkey",pmt_tubeid_to_channelkey);

	//-----------------------------------------------------------------
	//---------Familiarization with new Detectorkey layout-------------
	//-----------------------------------------------------------------

	std::cout <<"DetectorKey-LAPPD size: "<<detectorkey_to_lappdid.size()<<", ChannelKey-PMT size: "<<channelkey_to_pmtid.size()<<", ChannelKey-MRD size: "<<channelkey_to_mrdpmtid.size()<<", ChannelKey-FACC size: "<<channelkey_to_faccpmtid.size()<<std::endl;
	std::map<unsigned long, int>::iterator it_detectorkey;
	for (it_detectorkey=detectorkey_to_lappdid.begin();it_detectorkey!=detectorkey_to_lappdid.end();it_detectorkey++){
		std::cout <<"DetectorKey LAPPD: "<<it_detectorkey->first<<std::endl;
	}
	std::map<unsigned long, int>::iterator it_channelkeyMRD;
	for (it_channelkeyMRD=channelkey_to_mrdpmtid.begin();it_channelkeyMRD!=channelkey_to_mrdpmtid.end();it_channelkeyMRD++){
		std::cout <<"ChannelKey MRD: "<<it_channelkeyMRD->first<<std::endl;
	}
	std::map<unsigned long, int>::iterator it_channelkeyFACC;
	for (it_channelkeyFACC=channelkey_to_faccpmtid.begin();it_channelkeyFACC!=channelkey_to_faccpmtid.end();it_channelkeyFACC++){
		std::cout <<"ChannelKey FACC: "<<it_channelkeyFACC->first<<std::endl;
	}
	std::map<unsigned long, int>::iterator it_channelkey;
	for (it_channelkey=channelkey_to_pmtid.begin();it_channelkey!=channelkey_to_pmtid.end();it_channelkey++){
		std::cout <<"ChannelKey PMT: "<<it_channelkey->first<<std::endl;
	}
	std::map<int,unsigned long>::iterator it_pmt_to_channel;
	for (it_pmt_to_channel = pmt_tubeid_to_channelkey.begin();it_pmt_to_channel!=pmt_tubeid_to_channelkey.end();it_pmt_to_channel++){
		std::cout <<"PMT ID: "<< it_pmt_to_channel->first <<", Channelkey: "<< it_pmt_to_channel->second <<std::endl;
	}
	std::cout <<"PID: Num Tank PMTs: "<<n_tank_pmts<<", num MRD PMTs: "<<n_mrd_pmts<<", num Veto PMTs: "<<n_veto_pmts<<", num LAPPDs: "<<n_lappds<<std::endl;
	std::vector<std::map<unsigned long,Detector>* >* Detectors = geom->GetDetectors();
	std::cout <<"Detectors size: "<<Detectors->size()<<std::endl;

	//----------------------------------------------------------------------
	//-----------read in PMT x/y/z positions into vectors-------------------
	//----------------------------------------------------------------------

	for (int i_detector = 0; i_detector < Detectors->size(); i_detector++){

		//std::cout <<"Detector: "<<i_detector<<std::endl;
		std::map<unsigned long, Detector>* temp_map = Detectors->at(i_detector);
		//std::map<unsigned long, Detector>::iterator it;

		for(auto&& apmt : (*temp_map)){
			unsigned long detkey = apmt.first;
			unsigned long chankey;
			//std::cout <<"i_detector: "<<i_detector;
			//std::cout <<", detkey: "<<detkey<<", n_lappds: "<<n_lappds<<std::endl;
			if (i_detector == 0) chankey = detkey*n_channels_per_lappd;
			else if (i_detector == 1) chankey = pmt_tubeid_to_channelkey[detkey-n_lappds+1];
			else chankey = n_channels_per_lappd*n_lappds+(detkey-n_lappds);
			Detector* thedetector = geom->ChannelToDetector(chankey);
			//std::cout <<"chankey: "<<chankey<<std::endl;
			if (i_detector == 1) {
				//std::cout <<"PMT id: "<<channelkey_to_pmtid[chankey];
				Position position_PMT = thedetector->GetDetectorPosition();
				//std::cout <<"chankey: "<<chankey<<std::endl;
				//std::cout <<"filling PMT IDs: ID "<<detkey-n_lappds<<", wcsim ID: "<< channelkey_to_pmtid[chankey]<<std::endl;
				pmts_x.insert(std::pair<int,double>(channelkey_to_pmtid[chankey],position_PMT.X()-tank_center_x));
				pmts_y.insert(std::pair<int,double>(channelkey_to_pmtid[chankey],position_PMT.Y()-tank_center_y));
				pmts_z.insert(std::pair<int,double>(channelkey_to_pmtid[chankey],position_PMT.Z()-tank_center_z));
				//std::cout <<"WCSim ID: "<<channelkey_to_pmtid[chankey]<<", position: ("<<position_PMT.X()<<","<<position_PMT.Y()<<","<<position_PMT.Z()<<")"<<std::endl;
				if ((position_PMT.Y()-tank_center_y) > y_max) y_max = position_PMT.Y()-tank_center_y;
				else if ((position_PMT.Y()-tank_center_y) < y_min) y_min = position_PMT.Y()-tank_center_y;
			}
		}
	}

	//-----------------------------------------------------------------------
	//--------- Histograms for likelihood calculation (pid) -----------------
	//-----------------------------------------------------------------------

	//x-axis: energy (100 MeV - 2000 MeV)
	//y-axis: distance PMT - vertex (0.0-5.0m)
	//z-axis: angle track-PMT: (0 ... Pi)

	out_file_pdf = new TFile(filename_pdf.c_str(),"RECREATE");
	electrons_pe = new TH3F("electrons_pe","Electrons photoelectrons",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	electrons_time = new TH3F("electrons_time","Electrons times",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	n_electrons = new TH3F("n_electrons","# of bin contents (electrons)",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	n_electrons_total = new TH3F("n_electrons_total","# of bin contents (total, electrons)",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	muons_pe = new TH3F("muons_pe","Muons photoelectrons",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	muons_time = new TH3F("muons_time","Muons times",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	n_muons = new TH3F("n_muons","# of bin contents (muons)",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	n_muons_total = new TH3F("n_muons_total","# of bin contents (total, muons)",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	electrons_prob_unhit = new TH3F("electrons_prob_unhit","Probability for no hits (electrons)",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	muons_prob_unhit = new TH3F("muons_prob_unhit","Probability for no hits (muons)",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	electrons_pe_lappd = new TH3F("electrons_pe_lappd","Electrons photoelectrons LAPPD",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	electrons_time_lappd = new TH3F("electrons_time_lappd","Electrons times LAPPD",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	n_electrons_lappd = new TH3F("n_electrons_lappd","# of bin contents (electrons, LAPPD)",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	muons_pe_lappd = new TH3F("muons_pe_lappd","Muons photoelectrons LAPPD",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	muons_time_lappd = new TH3F("muons_time_lappd","Muons times LAPPD",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	n_muons_lappd = new TH3F("n_muons_lappd","# of bin contents (muons, LAPPD)",num_distance,0,5.0,num_angle,0,TMath::Pi(),num_energy,100,2000);
	log_tubenr = new TH1F("log_tubenr","Tube # from LoadWCSim",200,0,200);

	std::cout <<"num_distance: "<<num_distance<<", num_angle: "<<num_angle<<", num_energy: "<<num_energy<<std::endl;
	for (int i_distance = 0; i_distance < num_distance; i_distance++){
		std::stringstream ss_dist;
		ss_dist<<i_distance;
		std::string str_dist = ss_dist.str();
		//std::cout <<"iD = "<<i_distance<<", ";
		for (int i_angle =0; i_angle < num_angle; i_angle++){
			std::stringstream ss_angle;
			ss_angle<<i_angle;
			std::string str_angle = ss_angle.str();
			//std::cout <<"iA = "<<i_angle<<", ";
			for (int i_energy = 0; i_energy < num_energy; i_energy++){
				std::stringstream ss_energy;
				ss_energy << i_energy;
				std::string str_energy = ss_energy.str();
				//std::cout <<"iE = "<<i_energy<<std::endl;
				std::string histname_charge_e = "charge_e_x"+str_dist+"_y"+str_angle+"_z"+str_energy;
				std::string histname_time_e = "time_e_x"+str_dist+"_y"+str_angle+"_z"+str_energy;
				std::string histname_charge_mu = "charge_mu_x"+str_dist+"_y"+str_angle+"_z"+str_energy;
				std::string histname_time_mu = "time_mu_x"+str_dist+"_y"+str_angle+"_z"+str_energy;
				//std::cout <<"histogram names: "<<histname_charge_e<<", "<<histname_time_e<<std::endl;
				charge_electrons[i_energy][i_distance][i_angle] = new TH1F(histname_charge_e.c_str(),"Charge histogram (e-) for bin",1000,0,300);
				time_electrons[i_energy][i_distance][i_angle] = new TH1F(histname_time_e.c_str(),"Time histogram (e-) for bin",2000,0,2000);
				charge_muons[i_energy][i_distance][i_angle] = new TH1F(histname_charge_mu.c_str(),"Charge histogram (mu-) for bin",1000,0,300);
				time_muons[i_energy][i_distance][i_angle] = new TH1F(histname_time_mu.c_str(),"Time histogram (mu-) for bin",2000,0,2000);
			}
		}
	}

	all_theta = new TH1F("all_theta","Theta overview",100,0,TMath::Pi());
  	all_npmts = new TH1F("all_npmts","Number PMTs overview",143,0,143);
  	all_pe = new TH1F("all_pe","P.E.s overview",500,0,500);
  	all_pe_total = new TH1F("all_pe_total","Total P.E.s overview",500,0,12000);
  	all_time = new TH1F("all_time","Hit times overview",2000,0,2000);
  	all_theta_bary = new TH1F("all_theta_bary","Charge barycenter theta overview",100,0,TMath::Pi());
	all_theta_rms = new TH1F("all_theta_rms","Theta rms overview",100,0,TMath::Pi());
	all_theta_var = new TH1F("all_theta_var","Theta variance overview",100,0,TMath::Pi());
	all_theta_skew = new TH1F("all_theta_skew","Theta skewness overview",100,-10,2);
	all_theta_kurt = new TH1F("all_theta_kurt","Theta kurtosis overview",100,-10,2); 
 	all_ratio_ring = new TH1F("all_ratio_ring","Ratio Charge Ring/Total",100,0,1);
 	all_ratio_ring_noweight = new TH1F("all_ratio_ring_noweight","Ratio Charge Ring/Total (not weighted)",100,0,1);
	all_ratio_downstream = new TH1F("all_ratio_downstream","Ratio Charge Downstream/Total",100,0,1);
	all_theta_lappd = new TH1F("all_theta_lappd","Theta overview (LAPPD)",100,0,TMath::Pi());
  	all_theta_bary_lappd = new TH1F("all_theta_bary_lappd","Charge barycenter theta overview (LAPPD)",100,0,TMath::Pi());
  	all_pe_total_lappd = new TH1F("all_pe_total_lappd","Total P.E.s overview (LAPPD)",1000,0,1000);
  	all_time_lappd = new TH1F("all_time_lappd","Hit times overview (LAPPD)",1000,0,1000);
	all_theta_rms_lappd = new TH1F("all_theta_rms_lappd","Theta rms overview (LAPPD)",100,0,TMath::Pi());
	all_theta_var_lappd = new TH1F("all_theta_var_lappd","Theta variance overview (LAPPD)",100,0,TMath::Pi());
	all_theta_skew_lappd = new TH1F("all_theta_skew_lappd","Theta skewness overview (LAPPD)",100,-10,2);
	all_theta_kurt_lappd = new TH1F("all_theta_kurt_lappd","Theta kurtosis overview (LAPPD)",100,-10,2);
 	all_ratio_ring_lappd = new TH1F("all_ratio_ring_lappd","Ratio Charge Ring/Total (LAPPD)",100,0,1);
	all_ratio_downstream_lappd = new TH1F("all_ratio_downstream_lappd","Ratio Charge Downstream/Total (LAPPD)",100,0,1);
	all_pmts_mrd = new TH1F("all_pmts_mrd","Number of hit PMTs (MRD)",20,0,20);

	energy_step = 1900./num_energy;
	distance_step = 5./num_distance;
	angle_step = TMath::Pi()/double(num_angle);

	Log("CreatePIDPDF Tool: Initializing",v_message,verbosity);

	return true;
}


bool CreatePIDPDF::Execute(){

	Log("CreatePIDPDF Tool: Executing",v_debug,verbosity);

	get_ok = m_data->Stores.count("ANNIEEvent");
	if(!get_ok){
		Log("CreatePIDPDF Tool: No ANNIEEvent store!",v_error,verbosity);
		return false;
	};

	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",MCTriggernum);
	if(not get_ok){ Log("CreatePIDPDF Tool: Error retrieving MCTriggernum from ANNIEEvent!",v_error,verbosity); return false; }
	if(MCTriggernum>0){ Log("CreatePIDPDF Tool: Skipping delayed trigger",v_debug,verbosity); return true;}
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCFile",MCFile);
	if(not get_ok){ Log("CreatePIDPDF Tool: Error retrieving MCFile from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",MCEventNum);
	if(not get_ok){ Log("CreatePIDPDF Tool: Error retrieving MCEventNum from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventNumber",EventNumber);
	if(not get_ok){ Log("CreatePIDPDF Tool: Error retrieving EventNumber from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCParticles",MCParticles);
	if(not get_ok){ Log("CreatePIDPDF Tool: Error retrieving MCParticles,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCHits",MCHits);
	if(not get_ok){ Log("CreatePIDPDF Tool: Error retrieving MCHits,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCLAPPDHits",MCLAPPDHits);
	if(not get_ok){ Log("CreatePIDPDF Tool: Error retrieving MCLAPPDHits,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);
	if(not get_ok){ Log("CreatePIDPDF Tool: Error retrieving TDCData,true from ANNIEEvent!",v_error,verbosity); return false; }
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventTime",EventTime);
	if(not get_ok){ Log("CreatePIDPDF Tool: Error retrieving EventTime,true from ANNIEEvent!",v_error,verbosity); return false; }

	logmessage = "CreatePIDPDF Tool: Processing truth tracks and digits for "+MCFile
                            +", MCEvent "+to_string(MCEventNum)+", MCTrigger "+to_string(MCTriggernum);
	Log(logmessage,v_debug,verbosity);


	//-----------------------------------------------------------------------
	//--------------------- Finding primary particle ------------------------
	//-----------------------------------------------------------------------

	//std::cout <<"Start of Event!"<<std::endl;
	//std::cout <<"Trigger Number: "<<MCTriggernum<<", Event Number: "<<EventNumber<<std::endl;
	MCParticle primarymuon, primaryelectron; 
    	bool mufound=false;
	bool efound=false;

    	if(MCParticles){
		Log("CreatePIDPDF Tool: Num MCParticles = "+to_string(MCParticles->size()),v_message,verbosity);
		for(int particlei=0; particlei<MCParticles->size(); particlei++){
			MCParticle aparticle = MCParticles->at(particlei);
			if(aparticle.GetParentPdg()!=0) continue; 
			if (aparticle.GetFlag()!=0) continue;			//needed since primary particle is not created after a neutrino interaction
			if(aparticle.GetPdgCode()==13){	
				primarymuon = aparticle;
				mufound=true;
				is_electron=0;
				break;
			}
			else if (aparticle.GetPdgCode()==11){
				primaryelectron = aparticle;
				efound=true;
				is_electron=1;
				break;
			}
			else continue;
		}
	}
	else {
		Log("CreatePIDPDF Tool: No MCParticles in this event!",v_error,verbosity);
	}
	if (not mufound && not efound){
		Log("CreatePIDPDF Tool: Neither electrons nor muons in this event!",v_error,verbosity);
		return true;
	}

	//-----------------------------------------------------------------------
	//--------------------- Get vertex & direction --------------------------
	//-----------------------------------------------------------------------

	if (mufound){
		const Position neutrinovtx = primarymuon.GetStartVertex();
		vtx_x = neutrinovtx.X()-tank_center_x;
		vtx_y = neutrinovtx.Y()-tank_center_y;
		vtx_z = neutrinovtx.Z()-tank_center_z;
		const Direction muondirection = primarymuon.GetStartDirection();
		dir_x = muondirection.X();
		dir_y = muondirection.Y();
		dir_z = muondirection.Z();

		distWallVert = 1.98 - fabs(vtx_y);
		distWallHor = 1.524 - sqrt(vtx_x*vtx_x+vtx_z*vtx_z);
		distWallVert /= 1.98;
		distWallHor /= 1.524;
		distInnerStrHor = inner_radius - sqrt(vtx_x*vtx_x+vtx_z*vtx_z);
		distInnerStrVert = (vtx_y > 0)? (y_max - vtx_y) : (vtx_y - y_min);

		time_true = primarymuon.GetStartTime();
		t_reco = time_true-EventTime->GetNs();
		muonenergy = primarymuon.GetStartEnergy();
		logmessage = "CreatePIDPDF Tool: Interaction Vertex is at ("+to_string(neutrinovtx.X())
		            +", "+to_string(neutrinovtx.Y())+", "+to_string(neutrinovtx.Z())+")\n"
		            +"CreatePIDPDF Tool: Primary muon has energy "+to_string(muonenergy)+"GeV and direction ("
		            +to_string(muondirection.X())+", "+to_string(muondirection.Y())+", "+to_string(muondirection.Z())+" )";
		Log(logmessage,v_debug,verbosity);

	}
	else if (efound){
		const Position neutrinovtx = primaryelectron.GetStartVertex();
		vtx_x = neutrinovtx.X()-tank_center_x;
		vtx_y = neutrinovtx.Y()-tank_center_y;
		vtx_z = neutrinovtx.Z()-tank_center_z;
		const Direction electrondirection = primaryelectron.GetStartDirection();
		dir_x = electrondirection.X();
		dir_y = electrondirection.Y();
		dir_z = electrondirection.Z();
		
		distWallVert = 1.98 - fabs(vtx_y);
		distWallHor = 1.524 - sqrt(vtx_x*vtx_x+vtx_z*vtx_z);
		distWallVert /= 1.98;
		distWallHor /= 1.524;
		distInnerStrHor = inner_radius - sqrt(vtx_x*vtx_x+vtx_z*vtx_z);
		distInnerStrVert = (vtx_y > 0)? (y_max - vtx_y) : (vtx_y - y_min);

		time_true = primaryelectron.GetStartTime();
		t_reco = time_true-EventTime->GetNs();
		electronenergy = primaryelectron.GetStartEnergy();
		logmessage = "CreatePIDPDF Tool: Electron Interaction Vertex is at ("+to_string(neutrinovtx.X())
		            +", "+to_string(neutrinovtx.Y())+", "+to_string(neutrinovtx.Z())+")\n"
		            +"CreatePIDPDFTool: Primary electron has energy "+to_string(electronenergy)+"GeV and direction ("
		            +to_string(electrondirection.X())+", "+to_string(electrondirection.Y())+", "+to_string(electrondirection.Z())+" )";
		Log(logmessage,v_debug,verbosity);
	}

	//-----------------------------------------------------------------------
	//--------------------- Read out Monte Carlo hits -----------------------
	//-----------------------------------------------------------------------

	std::cout <<"CreatePIDPDF Tool: t_reco = "<<t_reco<<std::endl;
	if (MCHits && MCLAPPDHits && TDCData && t_reco > -100. && t_reco < 100. && distInnerStrVert > 0.2 && distInnerStrHor > 0.2){
		double totalQ=0.;
		double highestQ=0.;
		
		std::vector<double> singlePMT_dist, singlePMT_ang, singlePMT_q, singlePMT_t, singlePMT_distVert, singlePMT_distHor, singlePMT_ang_center;
		std::vector<int> singlePMT_id;

	/*	singlePMT_dist.assign(n_tank_pmts,0);
		singlePMT_ang.assign(n_tank_pmts,0);
		singlePMT_ang_center.assign(n_tank_pmts,0);
		singlePMT_q.assign(n_tank_pmts,0);
		singlePMT_t.assign(n_tank_pmts,0);
	*/
		std::vector<int> PMT_Hit(n_tank_pmts,0);
		Log("CreatePIDPDF Tool: Num PMT Digits = "+to_string(MCHits->size()),v_message,verbosity);
		for(std::pair<unsigned long,std::vector<Hit>>&& apair : *MCHits){
			unsigned long chankey = apair.first;
			//std::cout <<"chankey: "<<chankey;
			int wcsim_pmt_id = channelkey_to_pmtid[chankey];
			log_tubenr->Fill(wcsim_pmt_id);
			//std::cout <<"WCSim ID: "<<wcsim_pmt_id;
			Detector* thistube = geom->ChannelToDetector(chankey);
			Position detector_pos = thistube->GetDetectorPosition();
			//std::cout <<"detector position: ("<<detector_pos.X()-tank_center_x<<","<<detector_pos.Y()-tank_center_y<<","<<detector_pos.Z()-tank_center_z<<")"<<std::endl;
			if (thistube->GetDetectorElement()=="Tank"){
				std::vector<Hit>& hits = apair.second;
				double charge_pmt=0.;
				double t_pmt;
				//if (wcsim_pmt_id>140) std::cout <<"getting x/y/z for PMT id "<<wcsim_pmt_id<<std::endl;
				double x_pmt = pmts_x.at(wcsim_pmt_id);
				double y_pmt = pmts_y.at(wcsim_pmt_id);
				double z_pmt = pmts_z.at(wcsim_pmt_id);
				//std::cout <<", (x,y,z) = ("<<x_pmt<<","<<y_pmt<<","<<z_pmt<<")";
				double x_dir_pmt = x_pmt - vtx_x;
				double y_dir_pmt = y_pmt - vtx_y;
				double z_dir_pmt = z_pmt - vtx_z;
				double distance_pmt = sqrt(pow(x_pmt-vtx_x,2)+pow(y_pmt-vtx_y,2)+pow(z_pmt-vtx_z,2));
				double angle_pmt = acos((x_dir_pmt*dir_x+y_dir_pmt*dir_y+z_dir_pmt*dir_z)/distance_pmt);
				double distance_pmt_center = sqrt(x_pmt*x_pmt+y_pmt*y_pmt+z_pmt*z_pmt);
				//std::cout <<"distance = "<<distance_pmt<<", angle =  "<<angle_pmt<<std::endl;
				//std::cout <<", wcsim_pmt_id: "<<wcsim_pmt_id<<std::endl;
/*
				singlePMT_dist.at(wcsim_pmt_id-1) = distance_pmt;
				singlePMT_ang.at(wcsim_pmt_id-1) = angle_pmt;
				singlePMT_distVert.at(wcsim_pmt_id-1) = distWallVert;
				singlePMT_distHor.at(wcsim_pmt_id-1) = distWallHor;
*/				
				singlePMT_id.push_back(wcsim_pmt_id);
				singlePMT_dist.push_back(distance_pmt);
				singlePMT_ang.push_back(angle_pmt);
				singlePMT_distVert.push_back(distWallVert);
				singlePMT_distHor.push_back(distWallHor);
				all_theta->Fill(angle_pmt);

				if (parametric){		//apply parametric model from DigitBuilder

					std::vector<double> hitTimes;
					std::vector<double> hitCharges;
					for(Hit& ahit : hits){
						PMT_Hit.at(wcsim_pmt_id-1)=1;
						hitTimes.push_back(ahit.GetTime()*1.0); 
						hitCharges.push_back(ahit.GetCharge());
					}
					// Do median and sum
					std::sort(hitTimes.begin(), hitTimes.end());
					size_t timesize = hitTimes.size();
					if (timesize % 2 == 0){
						t_pmt = (hitTimes.at(timesize/2 - 1) + hitTimes.at(timesize/2))/2;
					} else {
						t_pmt = hitTimes.at(timesize/2);
					}

					for(std::vector<double>::iterator it = hitCharges.begin(); it != hitCharges.end(); ++it){
						charge_pmt += *it;
					}

					//fill arrays for 3D histograms
					if (mufound){
						if (muonenergy>=2000.) {/*std::cout <<"muon energy > 2000 MeV!"<<std::endl;*/ continue;}
						pe_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=charge_pmt;
						n_pe_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
						n_pe_muons_total[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
						times_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=t_pmt;
						n_times_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=1;
						charge_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(charge_pmt);
						time_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(t_pmt);
					}
					
					if (efound){
						if (electronenergy>=2000.) {/*std::cout <<"electron energy > 2000 MeV!"<<std::endl;*/ continue;}
						pe_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=charge_pmt;
						n_pe_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
						n_pe_electrons_total[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
						times_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=t_pmt;
						n_times_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
						charge_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(charge_pmt);
						time_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(t_pmt);
					}
					singlePMT_q.push_back(charge_pmt);
					singlePMT_t.push_back(t_pmt);
					all_pe->Fill(charge_pmt);
					all_time->Fill(t_pmt);
					if (charge_pmt>highestQ) highestQ = charge_pmt;
					totalQ += charge_pmt;				
					/*singlePMT_q.at(wcsim_pmt_id-1) = charge_pmt;
					singlePMT_t.at(wcsim_pmt_id-1) = t_pmt;*/
				} else {					//old ansatz: take mctruth information of every PMT hit separately
					int loop = 0;
					for (Hit& ahit: hits){			
						PMT_Hit.at(wcsim_pmt_id-1)=1;
						double charge_pmt = ahit.GetCharge();
						double t_pmt = ahit.GetTime();
						std::cout << "charge: "<<charge_pmt<<", time: "<<t_pmt<<std::endl;
						if (mufound){
							if (muonenergy>=2000.) {/*std::cout <<"muon energy > 2000 MeV!"<<std::endl;*/ continue;}
							pe_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=charge_pmt;
							n_pe_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
							n_pe_muons_total[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
							times_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=t_pmt;
							n_times_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=1;
							charge_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(charge_pmt);
							time_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(t_pmt);
						}
						if (efound){
							if (electronenergy>=2000.) {/*std::cout <<"electron energy > 2000 MeV!"<<std::endl;*/ continue;}
							pe_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=charge_pmt;
							n_pe_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
							n_pe_electrons_total[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
							times_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=t_pmt;
							n_times_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
							charge_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(charge_pmt);
							time_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(t_pmt);
						}
					//
					//take only first hit for csv file if not using parametric model
					//
					if (loop ==0){
						singlePMT_q.push_back(charge_pmt);
						singlePMT_t.push_back(t_pmt);
						if (charge_pmt > highestQ) highestQ = charge_pmt;
						totalQ+=charge_pmt;
						/*singlePMT_q.at(wcsim_pmt_id-1) = charge_pmt;
						singlePMT_t.at(wcsim_pmt_id-1) = t_pmt;*/
					}
					loop++;
					}
				}
			}
		}
		

		total_charge=0.;
		average_time=0.;
		nhits = 0;
		double barycenter_x=0.;
		double barycenter_y=0.;
		double barycenter_z=0.;
		average_dist=0.;
		double rms_angle=0.;
		double charge_ring = 0;
		double charge_downstream = 0;
		int nhits_ring=0;
		double ratio_ring = 0.;
		double ratio_downstream = 0.;
		double ratio_ring_noweight = 0.;
		double ratio_highestQ = 0.;

		double charge_bary = 0.;
		double rms_bary = 0.;
		double ratio_charge_bary = 0.;

		for (int i_charge = 0; i_charge < singlePMT_q.size(); i_charge++){
			total_charge+=singlePMT_q.at(i_charge);
			average_time+=singlePMT_t.at(i_charge);
			nhits++;
			average_dist+=singlePMT_dist.at(i_charge);
			rms_angle+=(singlePMT_ang.at(i_charge)*singlePMT_ang.at(i_charge));
			barycenter_x += (singlePMT_q.at(i_charge)*pmts_x.at(singlePMT_id.at(i_charge)));
			barycenter_y += (singlePMT_q.at(i_charge)*pmts_y.at(singlePMT_id.at(i_charge)));
			barycenter_z += (singlePMT_q.at(i_charge)*pmts_z.at(singlePMT_id.at(i_charge)));
			if (singlePMT_ang.at(i_charge) < TMath::Pi()/2.) charge_downstream+=(singlePMT_q.at(i_charge));
			if (singlePMT_ang.at(i_charge) < cherenkov_angle) {
				charge_ring+=(singlePMT_q.at(i_charge));
				nhits_ring++;
			}
		}
		if (singlePMT_t.size() != 0.){
			average_time /= singlePMT_t.size();
			average_dist /= singlePMT_t.size();
			rms_angle /= singlePMT_t.size();
			rms_angle = sqrt(rms_angle);
			barycenter_x /= total_charge;
			barycenter_y /= total_charge;
			barycenter_z /= total_charge;
			ratio_ring = charge_ring/total_charge;
			ratio_downstream = charge_downstream/total_charge;
			ratio_ring_noweight = double(nhits_ring)/double(nhits);
			ratio_highestQ = double(highestQ)/total_charge;
		}
		double x_dir_bary = barycenter_x - vtx_x;
		double y_dir_bary = barycenter_y - vtx_y;
		double z_dir_bary = barycenter_z - vtx_z;
		double barycenter_distance = sqrt(pow(barycenter_x-vtx_x,2)+pow(barycenter_y-vtx_y,2)+pow(barycenter_z-vtx_z,2));
		double barycenter_angle = acos((x_dir_bary*dir_x+y_dir_bary*dir_y+z_dir_bary*dir_z)/barycenter_distance);
		double barycenter_distance_center = sqrt(pow(barycenter_x,2)+pow(barycenter_y,2)+pow(barycenter_z,2));	
		for (int i_charge = 0; i_charge < singlePMT_q.size(); i_charge++){
			
			double x_pmt = pmts_x.at(singlePMT_id.at(i_charge));
			double y_pmt = pmts_y.at(singlePMT_id.at(i_charge));
			double z_pmt = pmts_z.at(singlePMT_id.at(i_charge));
			double distance_pmt = sqrt(x_pmt*x_pmt+y_pmt*y_pmt+z_pmt*z_pmt);		
			double angle_bary = acos((x_pmt*barycenter_x+y_pmt*barycenter_y+z_pmt*barycenter_z)/(distance_pmt*barycenter_distance_center));
			//std::cout <<"angle bary: "<<angle_bary/TMath::Pi()<<" Pi"<<std::endl;
			singlePMT_ang_center.push_back(angle_bary);
			if (angle_bary < TMath::Pi()/2.) charge_bary += singlePMT_q.at(i_charge);
			rms_bary += (angle_bary*angle_bary);	
		}
		if (singlePMT_q.size()!=0) rms_bary /= singlePMT_q.size();
		rms_bary = sqrt(rms_bary);
		ratio_charge_bary = charge_bary/total_charge;

		//std::cout <<"rms bary: "<<rms_bary/TMath::Pi()<<" Pi"<<std::endl;
		//std::cout <<"charge bary: "<<charge_bary<<", ratio Charge bary: "<<ratio_charge_bary<<std::endl;

		double variance = 0.;
		double skewness = 0.;
		double kurtosis = 0.;
		for (int i_charge = 0; i_charge < singlePMT_q.size();i_charge++){
			variance+=(pow(singlePMT_ang.at(i_charge)-barycenter_angle,2)*singlePMT_q.at(i_charge)/total_charge);
			skewness+=((pow(singlePMT_ang.at(i_charge)-barycenter_angle,3)*singlePMT_q.at(i_charge)/total_charge)/*/(pow((pow(singlePMT_ang.at(i_charge)-barycenter_angle,2)*singlePMT_q.at(i_charge)/total_charge),1.5))*/);
			kurtosis+=((pow(singlePMT_ang.at(i_charge)-barycenter_angle,4)*singlePMT_q.at(i_charge)/total_charge)/*/(pow((pow(singlePMT_ang.at(i_charge)-barycenter_angle,2)*singlePMT_q.at(i_charge)/total_charge),2))-3*/);	
		}
		kurtosis-=(3*variance*variance);

		all_npmts->Fill(nhits);
		all_pe_total->Fill(total_charge);
		all_theta_rms->Fill(rms_angle);
		all_theta_var->Fill(variance);
		all_theta_skew->Fill(skewness);
		all_theta_kurt->Fill(kurtosis);
		all_theta_bary->Fill(barycenter_angle);
		all_ratio_ring->Fill(ratio_ring);
		all_ratio_downstream->Fill(ratio_downstream);
		all_ratio_ring_noweight->Fill(ratio_ring_noweight);

		if (nhits > 10) {		//require at least 10 PMTs to be hit for PID
		
			if (efound) singleEvent_pid.push_back(11);
			else if (mufound) singleEvent_pid.push_back(13);
			if (efound) singleEvent_energy.push_back(electronenergy);
			else if (mufound) singleEvent_energy.push_back(muonenergy);
			singleEvent_distance.push_back(singlePMT_dist);
			singleEvent_angle.push_back(singlePMT_ang);
			singleEvent_distVert.push_back(distWallVert);
			singleEvent_distHor.push_back(distWallHor);
			singleEvent_distInnerStrVert.push_back(distInnerStrVert);
			singleEvent_distInnerStrHor.push_back(distInnerStrHor);
			singleEvent_charge.push_back(singlePMT_q);
			singleEvent_time.push_back(singlePMT_t);
			singleEvent_id.push_back(singlePMT_id);
			singleEvent_vtxX.push_back(vtx_x);
			singleEvent_vtxY.push_back(vtx_y);
			singleEvent_vtxZ.push_back(vtx_z);
			singleEvent_dirX.push_back(dir_x);
			singleEvent_dirY.push_back(dir_y);
			singleEvent_dirZ.push_back(dir_z);
			singleEvent_total_charge.push_back(total_charge);
			singleEvent_average_time.push_back(average_time);
			singleEvent_average_distance.push_back(average_dist);
			singleEvent_nhits.push_back(nhits);
			//std::cout << "barycenter (x,y,z) = ("<<barycenter_x<<","<<barycenter_y<<","<<barycenter_z<<"), distance: "<<barycenter_distance<<", angle: "<<barycenter_angle<<std::endl;
			singleEvent_barycenter.push_back(barycenter_angle);
			singleEvent_rms.push_back(rms_angle);
			singleEvent_variance.push_back(variance);
			singleEvent_skewness.push_back(skewness);
			singleEvent_kurtosis.push_back(kurtosis);
			singleEvent_evnum.push_back(EventNumber);
			singleEvent_ratioChargeRing.push_back(ratio_ring);
			singleEvent_rms_bary.push_back(rms_bary);
			singleEvent_ratioChargeBary.push_back(ratio_charge_bary);
			singleEvent_highestQ.push_back(ratio_highestQ);
		}	
		for (int i_pmt=1;i_pmt<=n_tank_pmts;i_pmt++){			//accomodate unhit PMTs in PDF building
			if (PMT_Hit.at(i_pmt-1)==0){
				double x_pmt = pmts_x.at(i_pmt);
                		double y_pmt = pmts_y.at(i_pmt);
                		double z_pmt = pmts_z.at(i_pmt);
                		double x_dir_pmt = x_pmt - vtx_x;
                		double y_dir_pmt = y_pmt - vtx_y;
                		double z_dir_pmt = z_pmt - vtx_z;
                		double distance_pmt = sqrt(pow(x_pmt-vtx_x,2)+pow(y_pmt-vtx_y,2)+pow(z_pmt-vtx_z,2));			
				double angle_pmt = acos((x_dir_pmt*dir_x+y_dir_pmt*dir_y+z_dir_pmt*dir_z)/distance_pmt);

				if (mufound){
					n_pe_muons_total[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
					//don't fill the unhit pmts in the same histogram (pdf bin)
					//charge_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(0);
				}
				if (efound){
					n_pe_electrons_total[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;			
					//charge_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]->Fill(0);
				}	
			}	
		}			 
	}else {
		cout <<"No MCHits"<<endl;			//in the case of 0 PMT hits, we do not count the times or charges anyhow
	}

	//-----------------------------------------------------------------------
	//--------------------- Read out MC - LAPPD hits ------------------------
	//-----------------------------------------------------------------------

	if (MCHits && MCLAPPDHits && TDCData && t_reco > -100. && t_reco < 100. && distInnerStrVert > 0.2 && distInnerStrHor > 0.2){

	        std::vector<double> lappd_ang, lappd_ang_center, lappd_q, lappd_t, lappd_x, lappd_y, lappd_z;
		Log("CreatePIDPDF Tool: Num LAPPD Digits = "+to_string(MCLAPPDHits->size()),v_message,verbosity);
		for(std::pair<unsigned long,std::vector<LAPPDHit>>&& apair : *MCLAPPDHits){

			unsigned long chankey = apair.first;
			Detector* det = geom->ChannelToDetector(chankey);
			if(det==nullptr){
			Log("CreatePIDPDF Tool: LAPPD Detector not found! ",v_message,verbosity);
				continue;
			}
			int detkey = det->GetDetectorID();
			int LAPPDId = detectorkey_to_lappdid.at(detkey);
			//std::cout << "LAPPD chankey: "<<chankey<<", det key: "<<detkey<<", Tube DetectorElement: "<<det->GetDetectorElement()<<std::endl;

			std::vector<LAPPDHit>& hits = apair.second;
			for(LAPPDHit& ahit : hits){
				double lappd_charge = ahit.GetCharge();
				lappd_charge = 1.0; 				//for now just use 1 hit for the LAPPD, since the digitization is not implemented yet
				std::vector<double> temp_pos = ahit.GetPosition();
				double x_lappd = temp_pos.at(0)-tank_center_x;
				double y_lappd = temp_pos.at(1)-tank_center_y;
				double z_lappd = temp_pos.at(2)-tank_center_z;
				double x_dir_lappd = x_lappd - vtx_x;
				double y_dir_lappd = y_lappd - vtx_y;
				double z_dir_lappd = z_lappd - vtx_z;
				double distance_lappd = sqrt(pow(x_lappd-vtx_x,2)+pow(y_lappd-vtx_y,2)+pow(z_lappd-vtx_z,2));
				double angle_lappd = acos((x_dir_lappd*dir_x+y_dir_lappd*dir_y+z_dir_lappd*dir_z)/distance_lappd);
				double t_lappd = ahit.GetTime();
				t_lappd = frand.Gaus(t_lappd, 0.1);		//Gaussian smearing for time information in LAPPD hits

				

				if (mufound){
					if (muonenergy>=2000.) continue;
					pe_muons_lappd[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]+=lappd_charge;
					n_pe_muons_lappd[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]++;
					times_muons_lappd[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]+=t_lappd;
				}

				if (efound){
					if (electronenergy>=2000.) continue;
					pe_electrons_lappd[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]+=lappd_charge;
					n_pe_electrons_lappd[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]++;
					times_electrons_lappd[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]+=t_lappd;
				}
				lappd_x.push_back(x_lappd);
				lappd_y.push_back(y_lappd);
				lappd_z.push_back(z_lappd);
				lappd_q.push_back(lappd_charge);
				lappd_t.push_back(t_lappd);
				lappd_ang.push_back(angle_lappd);

				all_time_lappd->Fill(t_lappd);
				all_theta_lappd->Fill(angle_lappd);						
			}

		}
                average_time_lappd=0.;
                nhits_lappd = 0;
                double barycenter_x_lappd=0.;
                double barycenter_y_lappd=0.;
                double barycenter_z_lappd=0.;
                double rms_angle_lappd=0.;
		double charge_ring_lappd = 0.;
		double charge_downstream_lappd = 0.;
		double ratio_ring_lappd = 0.;
		double ratio_downstream_lappd = 0.;
		double rms_bary_lappd=0.;
		double charge_bary_lappd=0.;
		double ratio_charge_bary_lappd=0.;


                for (int i_charge = 0; i_charge < lappd_q.size(); i_charge++){
                	average_time_lappd+=lappd_t.at(i_charge);
                        nhits_lappd++;
                        rms_angle_lappd+=(lappd_ang.at(i_charge)*lappd_ang.at(i_charge));
                        barycenter_x_lappd += (lappd_q.at(i_charge)*lappd_x.at(i_charge));
                        barycenter_y_lappd += (lappd_q.at(i_charge)*lappd_y.at(i_charge));
                        barycenter_z_lappd += (lappd_q.at(i_charge)*lappd_z.at(i_charge));
			if (lappd_ang.at(i_charge) < TMath::Pi()/2.) charge_downstream_lappd+=lappd_q.at(i_charge);
			if (lappd_ang.at(i_charge) < cherenkov_angle) charge_ring_lappd+=lappd_q.at(i_charge);
                }
                if (lappd_t.size() != 0.){
                         average_time_lappd /= lappd_t.size();
                         rms_angle_lappd /= lappd_t.size();
                         rms_angle_lappd = sqrt(rms_angle_lappd);
                         barycenter_x_lappd /= nhits_lappd;
                         barycenter_y_lappd /= nhits_lappd;
                         barycenter_z_lappd /= nhits_lappd;
                	 ratio_ring_lappd = charge_ring_lappd/double(nhits_lappd);
			 ratio_downstream_lappd = charge_downstream_lappd/double(nhits_lappd);
		}
                double x_dir_bary_lappd = barycenter_x_lappd - vtx_x;
                double y_dir_bary_lappd = barycenter_y_lappd - vtx_y;
                double z_dir_bary_lappd = barycenter_z_lappd - vtx_z;
                double barycenter_distance_lappd = sqrt(pow(barycenter_x_lappd-vtx_x,2)+pow(barycenter_y_lappd-vtx_y,2)+pow(barycenter_z_lappd-vtx_z,2));
                double barycenter_angle_lappd = acos((x_dir_bary_lappd*dir_x+y_dir_bary_lappd*dir_y+z_dir_bary_lappd*dir_z)/barycenter_distance_lappd);
                double barycenter_distance_lappd_center = sqrt(pow(barycenter_x_lappd,2)+pow(barycenter_y_lappd,2)+pow(barycenter_z_lappd,2));
		for (int i_charge =0.;i_charge < lappd_q.size(); i_charge++){

                        double x_lappd = lappd_x.at(i_charge);
                        double y_lappd = lappd_y.at(i_charge);
                        double z_lappd = lappd_z.at(i_charge);
                        double distance_lappd = sqrt(x_lappd*x_lappd+y_lappd*y_lappd+z_lappd*z_lappd);
                        double angle_bary_lappd = acos((x_lappd*barycenter_x_lappd+y_lappd*barycenter_y_lappd+z_lappd*barycenter_z_lappd)/(distance_lappd*barycenter_distance_lappd_center));
                        //std::cout <<"barycenter_distance_lappd_center: "<<barycenter_distance_lappd_center<<std::endl;
			//std::cout <<"distance_lappd: "<<distance_lappd<<std::endl;
			//std::cout <<"angle bary LAPPD: "<<angle_bary_lappd/TMath::Pi()<<" Pi"<<std::endl;
                        lappd_ang_center.push_back(angle_bary_lappd);
                        if (angle_bary_lappd < TMath::Pi()/2.) charge_bary_lappd += lappd_q.at(i_charge);
                        rms_bary_lappd += (angle_bary_lappd*angle_bary_lappd);
			//std::cout <<"rms_bary_lappd: "<<rms_bary_lappd<<std::endl;
		}

		//std::cout <<"lappd_q.size: "<<lappd_q.size()<<std::endl;
		//std::cout <<"lappd_q[0]: "<<lappd_q.at(0)<<", lappd x: "<<lappd_x.at(0)<<std::endl;
		//std::cout <<"rms_angle: "<<rms_bary_lappd<<std::endl;
		if (lappd_q.size()!=0) rms_bary_lappd/=lappd_q.size();
		rms_bary_lappd = sqrt(rms_bary_lappd);
		//std::cout <<"rms bary lappd: "<<rms_bary_lappd/TMath::Pi()<<" Pi"<<std::endl;
		ratio_charge_bary_lappd = charge_bary_lappd / double(nhits_lappd);
		//std::cout <<"charge bary lappd: "<<charge_bary_lappd<<", ratio charge bary lappd: "<<ratio_charge_bary_lappd<<std::endl;

		double variance_lappd = 0.;
                double skewness_lappd = 0.;
                double kurtosis_lappd = 0.;
                for (int i_charge = 0; i_charge < lappd_q.size();i_charge++){
                        variance_lappd+=(pow(lappd_ang.at(i_charge)-barycenter_angle_lappd,2)*lappd_q.at(i_charge)/nhits_lappd);
                        skewness_lappd+=((pow(lappd_ang.at(i_charge)-barycenter_angle_lappd,3)*lappd_q.at(i_charge)/nhits_lappd));
                        kurtosis_lappd+=((pow(lappd_ang.at(i_charge)-barycenter_angle_lappd,4)*lappd_q.at(i_charge)/nhits_lappd));
                }
                kurtosis_lappd-=(3*variance_lappd*variance_lappd);
		
		all_pe_total_lappd->Fill(nhits_lappd);
		all_theta_rms_lappd->Fill(rms_angle_lappd);
		all_theta_var_lappd->Fill(variance_lappd);
		all_theta_skew_lappd->Fill(skewness_lappd);
		all_theta_kurt_lappd->Fill(kurtosis_lappd);
		all_theta_bary_lappd->Fill(barycenter_angle_lappd);
		all_ratio_ring_lappd->Fill(ratio_ring_lappd);
		all_ratio_downstream_lappd->Fill(ratio_downstream_lappd);
		
		//fill LAPPD hit information into vectors
		if (nhits > 10) { //need to apply the same selection cuts as for PMT vector containers
                        singleEvent_nhits_lappd.push_back(nhits_lappd);
                        singleEvent_average_time_lappd.push_back(average_time_lappd);
                        singleEvent_barycenter_lappd.push_back(barycenter_angle_lappd);
                        singleEvent_rms_lappd.push_back(rms_angle_lappd);
                        singleEvent_variance_lappd.push_back(variance_lappd);
                       	singleEvent_skewness_lappd.push_back(skewness_lappd);
                        singleEvent_kurtosis_lappd.push_back(kurtosis_lappd);
			singleEvent_ratioChargeRing_lappd.push_back(ratio_ring_lappd);
			singleEvent_rms_bary_lappd.push_back(rms_bary_lappd);
			singleEvent_ratioChargeBary_lappd.push_back(ratio_charge_bary_lappd);
		}
	} else std::cout <<"No MCLAPPDHits"<<endl;

        //-----------------------------------------------------------------------
        //----------------------- Read out MC - MRD hits ------------------------
        //-----------------------------------------------------------------------	

  	if(!TDCData){
    		std::cout<<"No TDC data to plot in Event Display!"<<std::endl;
  	} else if (MCHits && MCLAPPDHits && t_reco > -100. && t_reco < 100. && distInnerStrVert > 0.2 && distInnerStrHor > 0.2){
    		if(TDCData->size()==0){
      			std::cout<<"No TDC hits to plot in Event Display!"<<std::endl;
			all_pmts_mrd->Fill(0);
			if (nhits > 10) singleEvent_num_mrd_hits.push_back(0);
    		} else {
			int num_mrd_pmts=0;
      			for(auto&& anmrdpmt : (*TDCData)){
        			unsigned long chankey = anmrdpmt.first;
        			Detector *thedetector = geom->ChannelToDetector(chankey);
				if(thedetector->GetDetectorElement()!="MRD") {
            				continue;                 // this is a veto hit, not an MRD hit.
        			}
      				num_mrd_pmts++;
			}
    		all_pmts_mrd->Fill(num_mrd_pmts);
		if (nhits > 10) singleEvent_num_mrd_hits.push_back(num_mrd_pmts);
    		}
  	}


	std::cout <<"CreatePIDPDF tool: Execution finished"<<std::endl;

	return true;

}


bool CreatePIDPDF::Finalise(){

  std::cout <<"CreatePIDPDF tool: Finalisation started"<<std::endl;
 
	//-----------------------------------------------------------------------
	//--------------------- Fill & Write PDFs -------------------------------
	//-----------------------------------------------------------------------

  for (int i_energy=0;i_energy<num_energy;i_energy++){
	for (int i_distance=0;i_distance<num_distance;i_distance++){
		for (int i_angle=0;i_angle<num_angle;i_angle++){
			n_electrons->SetBinContent(i_distance+1,i_angle+1,i_energy+1,n_pe_electrons[i_energy][i_distance][i_angle]);
			n_muons->SetBinContent(i_distance+1,i_angle+1,i_energy+1,n_pe_muons[i_energy][i_distance][i_angle]);
			n_electrons_total->SetBinContent(i_distance+1,i_angle+1,i_energy+1,n_pe_electrons_total[i_energy][i_distance][i_angle]);
			n_muons_total->SetBinContent(i_distance+1,i_angle+1,i_energy+1,n_pe_muons_total[i_energy][i_distance][i_angle]);
			n_electrons_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,n_pe_electrons_lappd[i_energy][i_distance][i_angle]);
			n_muons_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,n_pe_muons_lappd[i_energy][i_distance][i_angle]);
			if (n_pe_electrons_total[i_energy][i_distance][i_angle]==0) n_pe_electrons_total[i_energy][i_distance][i_angle]=1;
			if (n_pe_muons_total[i_energy][i_distance][i_angle]==0) n_pe_muons_total[i_energy][i_distance][i_angle]=1;
			electrons_prob_unhit->SetBinContent(i_distance+1,i_angle+1,i_energy+1,(n_pe_electrons_total[i_energy][i_distance][i_angle]-n_pe_electrons[i_energy][i_distance][i_angle])/n_pe_electrons_total[i_energy][i_distance][i_angle]);
			muons_prob_unhit->SetBinContent(i_distance+1,i_angle+1,i_energy+1,(n_pe_muons_total[i_energy][i_distance][i_angle]-n_pe_muons[i_energy][i_distance][i_angle])/n_pe_muons_total[i_energy][i_distance][i_angle]);
			if (n_pe_electrons[i_energy][i_distance][i_angle]==0) n_pe_electrons[i_energy][i_distance][i_angle]=1;
			if (n_pe_muons[i_energy][i_distance][i_angle]==0) n_pe_muons[i_energy][i_distance][i_angle]=1;
			if (n_times_electrons[i_energy][i_distance][i_angle]==0) n_times_electrons[i_energy][i_distance][i_angle]=1;
			if (n_times_muons[i_energy][i_distance][i_angle]==0) n_times_muons[i_energy][i_distance][i_angle]=1;
			electrons_pe->SetBinContent(i_distance+1,i_angle+1,i_energy+1,pe_electrons[i_energy][i_distance][i_angle]/n_pe_electrons[i_energy][i_distance][i_angle]);
			electrons_time->SetBinContent(i_distance+1,i_angle+1,i_energy+1,times_electrons[i_energy][i_distance][i_angle]/n_times_electrons[i_energy][i_distance][i_angle]);
			muons_pe->SetBinContent(i_distance+1,i_angle+1,i_energy+1,pe_muons[i_energy][i_distance][i_angle]/n_pe_muons[i_energy][i_distance][i_angle]);
			muons_time->SetBinContent(i_distance+1,i_angle+1,i_energy+1,times_muons[i_energy][i_distance][i_angle]/n_times_muons[i_energy][i_distance][i_angle]);
			if (n_pe_electrons_lappd[i_energy][i_distance][i_angle]==0) n_pe_electrons_lappd[i_energy][i_distance][i_angle]=1;
			if (n_pe_muons_lappd[i_energy][i_distance][i_angle]==0) n_pe_muons_lappd[i_energy][i_distance][i_angle]=1;
			electrons_pe_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,pe_electrons_lappd[i_energy][i_distance][i_angle]/n_pe_electrons_lappd[i_energy][i_distance][i_angle]);
			electrons_time_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,times_electrons_lappd[i_energy][i_distance][i_angle]/n_pe_electrons_lappd[i_energy][i_distance][i_angle]);
			muons_pe_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,pe_muons_lappd[i_energy][i_distance][i_angle]/n_pe_muons_lappd[i_energy][i_distance][i_angle]);
			muons_time_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,times_muons_lappd[i_energy][i_distance][i_angle]/n_pe_muons_lappd[i_energy][i_distance][i_angle]);
		}
	}
  }
  
  std::cout <<"Writing ROOT file (overview plots)..."<<std::endl;
  TFile *file_overview = new TFile(filename_root_overview.c_str(),"RECREATE");
  file_overview->cd();
  all_npmts->Write();
  all_time->Write();
  all_pe->Write();
  all_pe_total->Write();
  all_theta->Write();
  all_theta_bary->Write();
  all_theta_rms->Write();
  all_theta_var->Write();
  all_theta_skew->Write();
  all_theta_kurt->Write();
  all_ratio_ring->Write();
  all_ratio_downstream->Write();
  all_ratio_ring_noweight->Write();
  all_time_lappd->Write();
  all_pe_total_lappd->Write();
  all_theta_lappd->Write();
  all_theta_bary_lappd->Write();
  all_theta_rms_lappd->Write();
  all_theta_var_lappd->Write();
  all_theta_skew_lappd->Write();
  all_theta_kurt_lappd->Write();
  all_ratio_ring_lappd->Write();
  all_ratio_downstream_lappd->Write();
  all_pmts_mrd->Write();
  file_overview->Close();
	
  std::cout <<"writing hists to TFile..."<<std::endl;

  out_file_pdf->cd();
  electrons_time->Write();
  electrons_pe->Write();
  muons_time->Write();
  muons_pe->Write();
  n_electrons->Write();
  n_muons->Write();
  n_electrons_total->Write();
  n_muons_total->Write();
  electrons_prob_unhit->Write();
  muons_prob_unhit->Write();
  electrons_time_lappd->Write();
  electrons_pe_lappd->Write();
  muons_time_lappd->Write();
  muons_pe_lappd->Write();
  n_electrons_lappd->Write();
  n_muons_lappd->Write();
  log_tubenr->Write();

  if (filename_root != "None"){
  out_file_root = new TFile(filename_root.c_str(),"RECREATE");
  out_file_root->cd(); 
  for (int i_distance = 0; i_distance < num_distance; i_distance++){
  	for (int i_angle = 0; i_angle < num_angle; i_angle++){
		for (int i_energy = 0; i_energy < num_energy; i_energy++){
			charge_electrons[i_energy][i_distance][i_angle]->Write();
			time_electrons[i_energy][i_distance][i_angle]->Write();
			charge_muons[i_energy][i_distance][i_angle]->Write();
			time_muons[i_energy][i_distance][i_angle]->Write();
		}
	}
  }
  out_file_root->Close();
  }
  out_file_pdf->Close(); 

  //
  //fill csv file
  //

  std::cout <<"writing CSV file..."<<std::endl;

  if (filename_csv != "None"){
  ofstream csv_file(filename_csv.c_str());
  for (int i_vector=0; i_vector < singleEvent_pid.size(); i_vector++){
	for (int i_pmt=0; i_pmt < singleEvent_charge.at(i_vector).size(); i_pmt++){
		csv_file << singleEvent_distance.at(i_vector).at(i_pmt) <<","<< singleEvent_angle.at(i_vector).at(i_pmt) <<","<< singleEvent_charge.at(i_vector).at(i_pmt) <<",";
	}
	for (int i_pmt =0; i_pmt < n_tank_pmts; i_pmt++){
	//	std::cout <<"checking if PMT "<<i_pmt+1<<"is included in the vector already..."<<std::endl;
		if (std::find(singleEvent_id.at(i_vector).begin(), singleEvent_id.at(i_vector).end(), i_pmt+1) != singleEvent_id.at(i_vector).end()) {
			//already included those events
	//		std::cout <<"yup"<<std::endl;
		}else {
			//include zero hits
	//		std::cout <<"nope"<<std::endl;
			double distance0 = sqrt(pow(pmts_x.at(i_pmt+1)-singleEvent_vtxX.at(i_vector),2)+pow(pmts_y.at(i_pmt+1)-singleEvent_vtxY.at(i_vector),2)+pow(pmts_z.at(i_pmt+1)-singleEvent_vtxZ.at(i_vector),2));
			double angle0 = acos(((pmts_x.at(i_pmt+1)-singleEvent_vtxX.at(i_vector))*singleEvent_dirX.at(i_vector)+(pmts_y.at(i_pmt+1)-singleEvent_vtxY.at(i_vector))*singleEvent_dirY.at(i_vector)+(pmts_z.at(i_pmt+1)-singleEvent_vtxZ.at(i_vector))*singleEvent_dirZ.at(i_vector))/distance0);
			double charge0 = 0.;
			csv_file << distance0 <<","<<angle0 <<","<< charge0<<",";
		}
	}
	csv_file << singleEvent_distHor.at(i_vector) <<","<< singleEvent_distVert.at(i_vector)<<","<<singleEvent_energy.at(i_vector) << std::endl;
  	//std::cout <<"written to csv file..."<<std::endl;
	}
  //std::cout <<"closing file..."<<std::endl;
  csv_file.close();
  //std::cout <<"closed..."<<std::endl;
  }

  //
  //fill csv file (average)
  //

  std::cout <<"writing CSV file (average)..."<<std::endl;

  if (filename_csv_average != "None"){
  std::cout <<"Sizes of vectors: "<<std::endl;
  std::cout <<"singleEvent_pid: "<<singleEvent_pid.size()<<std::endl;
  std::cout <<"singleEvent_nhits: "<<singleEvent_nhits.size()<<std::endl;
  std::cout <<"singleEvent_total_charge: "<<singleEvent_total_charge.size()<<std::endl;
  std::cout <<"singleEvent_average_time: "<<singleEvent_average_time.size()<<std::endl;
  std::cout <<"singleEvent_barycenter: "<<singleEvent_barycenter.size()<<std::endl;
  std::cout <<"singleEvent_rms: "<<singleEvent_rms.size()<<std::endl;
  std::cout <<"singleEvent_variance: "<<singleEvent_variance.size()<<std::endl;
  std::cout <<"singleEvent_skewness: "<<singleEvent_skewness.size()<<std::endl;
  std::cout <<"singleEvent_kurtosis: "<<singleEvent_kurtosis.size()<<std::endl;
  std::cout <<"singleEvent_ratioChargeRing: "<<singleEvent_ratioChargeRing.size()<<std::endl;
  std::cout <<"singleEvent_distVert: "<<singleEvent_distVert.size()<<std::endl;
  std::cout <<"singleEvent_distHor: "<<singleEvent_distHor.size()<<std::endl;
  std::cout <<"singleEvent_distInnerStrVert: "<<singleEvent_distInnerStrVert.size()<<std::endl;
  std::cout <<"singleEvent_distInnerStrHor: "<<singleEvent_distInnerStrHor.size()<<std::endl;
  std::cout <<"singleEvent_average_distance: "<<singleEvent_average_distance.size()<<std::endl;
  std::cout <<"singleEvent_energy: "<<singleEvent_energy.size()<<std::endl;
  std::cout <<"singleEvent_nhits_lappd: "<<singleEvent_nhits_lappd.size()<<std::endl;
  std::cout <<"singleEvent_average_time_lappd: "<<singleEvent_average_time_lappd.size()<<std::endl;
  std::cout <<"singleEvent_barycenter_lappd: "<<singleEvent_barycenter_lappd.size()<<std::endl;
  std::cout <<"singleEvent_rms_lappd: "<<singleEvent_rms_lappd.size()<<std::endl;
  std::cout <<"singleEvent_variance_lappd: "<<singleEvent_variance_lappd.size()<<std::endl;
  std::cout <<"singleEvent_skewness_lappd: "<<singleEvent_skewness_lappd.size()<<std::endl;
  std::cout <<"singleEvent_kurtosis_lappd: "<<singleEvent_kurtosis_lappd.size()<<std::endl;
  std::cout <<"singleEvent_ratioChargeRing_lappd: "<<singleEvent_ratioChargeRing_lappd.size()<<std::endl;
  std::cout <<"singleEvent_num_mrd_hits: "<<singleEvent_num_mrd_hits.size()<<std::endl;
  std::cout <<"singleEvent_evnum: "<<singleEvent_evnum.size()<<std::endl;
  std::cout <<"singleEvent_highestQ: "<<singleEvent_highestQ.size()<<std::endl;

  ofstream csv_file_average(filename_csv_average.c_str());
  for (int i_vector = 0; i_vector < singleEvent_pid.size(); i_vector++){
	std::cout <<"i_vector (average): "<<i_vector<<std::endl;
	csv_file_average << singleEvent_nhits.at(i_vector) <<","<<singleEvent_total_charge.at(i_vector)<<","<<singleEvent_average_time.at(i_vector)<<","<<singleEvent_barycenter.at(i_vector)<<","<<singleEvent_rms.at(i_vector)<<","<<singleEvent_variance.at(i_vector)<<","<<singleEvent_skewness.at(i_vector)<<","<<singleEvent_kurtosis.at(i_vector)<<","<<singleEvent_ratioChargeRing.at(i_vector)<<","<<singleEvent_distVert.at(i_vector)<<","<<singleEvent_distHor.at(i_vector)<<","<<singleEvent_distInnerStrVert.at(i_vector)<<","<<singleEvent_distInnerStrHor.at(i_vector)<<","<<singleEvent_average_distance.at(i_vector)<<","<<singleEvent_energy.at(i_vector)<<","<<singleEvent_nhits_lappd.at(i_vector)<<","<<singleEvent_average_time_lappd.at(i_vector)<<","<<singleEvent_barycenter_lappd.at(i_vector)<<","<<singleEvent_rms_lappd.at(i_vector)<<","<<singleEvent_variance_lappd.at(i_vector)<<","<<singleEvent_skewness_lappd.at(i_vector)<<","<<singleEvent_kurtosis_lappd.at(i_vector)<<","<<singleEvent_ratioChargeRing_lappd.at(i_vector)<<","<<singleEvent_highestQ.at(i_vector)<<","<<singleEvent_num_mrd_hits.at(i_vector)<<","<<singleEvent_evnum.at(i_vector)<<std::endl;
  }
  csv_file_average.close();
  }

  std::cout <<"Cross check data sizes: PMT singleEvent vectors: "<<singleEvent_nhits.size()<<", LAPPD singleEvent vectors: "<<singleEvent_nhits_lappd.size()<<std::endl;

  std::cout <<"writing CSV file (data)..."<<std::endl;

  if (filename_csv_data != "None"){
  ofstream csv_file_data(filename_csv_data.c_str());
  for (int i_vector = 0; i_vector < singleEvent_pid.size(); i_vector++){
	std::cout <<"i_vector (data): "<<i_vector<<std::endl; 
       csv_file_data << singleEvent_nhits.at(i_vector) <<","<<singleEvent_total_charge.at(i_vector)<<","<<singleEvent_average_time.at(i_vector)<<","<<singleEvent_rms_bary.at(i_vector)<<","<<singleEvent_ratioChargeBary.at(i_vector)<<","<<singleEvent_nhits_lappd.at(i_vector)<<","<<singleEvent_average_time_lappd.at(i_vector)<<","<<singleEvent_rms_bary_lappd.at(i_vector)<<","<<singleEvent_ratioChargeBary_lappd.at(i_vector)<<","<<singleEvent_highestQ.at(i_vector)<<","<<singleEvent_num_mrd_hits.at(i_vector)<<","<<singleEvent_evnum.at(i_vector)<<std::endl;
  }
  csv_file_data.close();
  }
  
  std::cout << "CSV Data file also written"<<std::endl;


 std::cout <<"CreatePIDPDF tool: Finalisation complete"<<std::endl;

  return true;
}

