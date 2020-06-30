#include "SelectMuonSample.h"

SelectMuonSample::SelectMuonSample():Tool(){}


bool SelectMuonSample::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////

  m_variables.Get("verbosity",verbosity);
  m_variables.Get("OutputFile",filename);
  m_variables.Get("OutputFile_PDF",filename_pdf);

  m_data->Stores["ANNIEEvent"]->Header->Get("AnnieGeometry",geom);

  Position detector_center = geom.GetTankCentre();
  tank_center_x = detector_center.X();
  tank_center_y = detector_center.Y();
  tank_center_z = detector_center.Z();

  n_tank_pmts = geom.GetNumTankPMTs();
  n_mrd_pmts = geom.GetNumMrdPMTs();
  n_veto_pmts = geom.GetNumVetoPMTs();

  std::map<ChannelKey,Detector>::iterator it;
  std::map<ChannelKey,Detector> *Detectors = geom.GetDetectors();
  ChannelKey ck;
  Position position_PMT;
  int detector_id;

  for (it = Detectors->begin(); it!= Detectors->end();it++){
	
	ck = it->first;
	if (it->second.GetDetectorElement() !="Tank" || it->second.GetDetectorType()=="lappd_v1"){
		continue;
	}
    	std::cout <<it->second.GetDetectorType()<<std::endl;

	position_PMT = it->second.GetDetectorPosition();
	detector_id = it->second.GetDetectorId();	
	pmts_x.insert(std::pair<int,double>(detector_id,position_PMT.X()-tank_center_x));
	pmts_y.insert(std::pair<int,double>(detector_id,position_PMT.Y()-tank_center_y));
	pmts_z.insert(std::pair<int,double>(detector_id,position_PMT.Z()-tank_center_z));

  }  

  std::cout <<"outputfilename: "<<filename<<std::endl;

  out_file_pdf = new TFile(filename_pdf.c_str(),"RECREATE");
  out_file = new TFile(filename.c_str(),"RECREATE");
  h_muon_energy = new TH1F("h_muon_energy","Muon energy",100,1,0);
  h_muon_pe = new TH1F("h_muon_pe","Muon detected p.e.",200,0,1000);
  h_muon_lappd_pe = new TH1F("h_muon_lappd_pe","Muon detected hits (LAPPD)",100,1,0);
  h_muon_vertex = new TH3F("h_muon_vertex","Muon vertex",100,1,0,100,1,0,100,1,0);
  h_muon_dir = new TH3F("h_muon_dir","Muon direction",100,1,0,100,1,0,100,1,0);
  h_electron_energy = new TH1F("h_electron_energy","Electron energy",100,1,0);
  h_electron_pe = new TH1F("h_electron_pe","Electron detected p.e.",100,0,500);
  h_electron_lappd_pe = new TH1F("h_electron_lappd_pe","Electron detected hits (LAPPD)",100,1,0);
  h_electron_vertex = new TH3F("h_electron_vertex","Electron vertex",100,1,0,100,1,0,100,1,0);
  h_electron_dir = new TH3F("h_electron_dir","Electron direction",100,1,0,100,1,0,100,1,0);
  
  //histograms for likelihood calculation (pid)
  //x-axis: energy (100 MeV - 2000 MeV)
  //y-axis: distance PMT - vertex (0.0-5.0m)
  //z-axis: angle track-PMT: (0 ... Pi)
  electrons_pe = new TH3F("electrons_pe","Electrons photoelectrons",100,0,5.0,100,0,TMath::Pi(),21,100,2000);
  electrons_time = new TH3F("electrons_time","Electrons times",100,0,5.0,100,0,TMath::Pi(),21,100,2000);
  muons_pe = new TH3F("muons_pe","Muons photoelectrons",100,0,5.0,100,0,TMath::Pi(),21,100,2000);
  muons_time = new TH3F("muons_time","Muons times",100,0,5.0,100,0,TMath::Pi(),21,100,2000);
  electrons_pe_lappd = new TH3F("electrons_pe_lappd","Electrons photoelectrons LAPPD",100,0,5.0,100,0,TMath::Pi(),21,100,2000);
  electrons_time_lappd = new TH3F("electrons_time_lappd","Electrons times LAPPD",100,0,5.0,100,0,TMath::Pi(),21,100,2000);
  muons_pe_lappd = new TH3F("muons_pe_lappd","Muons photoelectrons LAPPD",100,0,5.0,100,0,TMath::Pi(),21,100,2000);
  muons_time_lappd = new TH3F("muons_time_lappd","Muons times LAPPD",100,0,5.0,100,0,TMath::Pi(),21,100,2000);

  tree_beam = new TTree("tree_beam","Tree of beam interaction");
  TBranch *tree_mu_energy = tree_beam->Branch("mu_energy",&mu_energy,"mu_energy/D");
  TBranch *tree_mu_pe = tree_beam->Branch("mu_pe",&mu_pe,"mu_pe/D");
  TBranch *tree_mu_pe_lappd = tree_beam->Branch("mu_pe_lappd",&mu_pe_lappd,"mu_pe_lappd/D");
  TBranch *tree_mu_pos = tree_beam->Branch("mu_pos",&mu_pos,"mu_pos[3]/D");
  TBranch *tree_mu_stop_pos = tree_beam->Branch("mu_stop_pos",&mu_stop_pos,"mu_stop_pos[3]/D");
  TBranch *tree_mu_dir = tree_beam->Branch("mu_dir",&mu_dir,"mu_dir[3]/D");
  TBranch *tree_mu_start_stop = tree_beam->Branch("mu_start_stop",&mu_start_stop,"mu_start_stop[3]/D");
  TBranch *tree_is_electron = tree_beam->Branch("is_electron",&is_electron,"is_electron/I");
  TBranch *tree_mu_pe_single = tree_beam->Branch("mu_pe_single",&mu_pe_single);
  TBranch *tree_mu_angle_single = tree_beam->Branch("mu_angle_single",&mu_angle_single);
  TBranch *tree_mu_distance_single = tree_beam->Branch("mu_distance_single",&mu_distance_single);
  TBranch *tree_mu_pmt_single = tree_beam->Branch("mu_pmt_single",&mu_pmt_single);
  TBranch *tree_mu_lappd = tree_beam->Branch("mu_lappd",&mu_lappd);
  TBranch *tree_mu_distance_lappd = tree_beam->Branch("mu_distance_lappd",&mu_distance_lappd);
  TBranch *tree_mu_angle_lappd = tree_beam->Branch("mu_angle_lappd",&mu_angle_lappd);
  TBranch *tree_time_pmt = tree_beam->Branch("time_pmt",&time_pmt);
  TBranch *tree_time_lappd = tree_beam->Branch("time_lappd",&time_lappd);


  energy_step = 1900./21;
  distance_step = 5./100;
  angle_step = TMath::Pi()/100.;

  Log("SelectMuonSample Tool: Initializing",v_message,verbosity);

  return true;
}


bool SelectMuonSample::Execute(){

   	Log("SelectMuonSample Tool: Executing",v_debug,verbosity);
        get_ok = m_data->Stores.count("ANNIEEvent");
        if(!get_ok){
                Log("SelectMuonSample Tool: No ANNIEEvent store!",v_error,verbosity);
                return false;
        };
	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCTriggernum",MCTriggernum);
        if(not get_ok){ Log("SelectMuonSample Tool: Error retrieving MCTriggernum from ANNIEEvent!",v_error,verbosity); return false; }
 	if(MCTriggernum>0){ Log("SelectMuonSample Tool: Skipping delayed trigger",v_debug,verbosity); return true;}
 	get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCFile",MCFile);
        if(not get_ok){ Log("SelectMuonSample Tool: Error retrieving MCFile from ANNIEEvent!",v_error,verbosity); return false; }
        get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCEventNum",MCEventNum);
        if(not get_ok){ Log("SelectMuonSample Tool: Error retrieving MCEventNum from ANNIEEvent!",v_error,verbosity); return false; }
        get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventNumber",EventNumber);
        if(not get_ok){ Log("SelectMuonSample Tool: Error retrieving EventNumber from ANNIEEvent!",v_error,verbosity); return false; }
        get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCParticles",MCParticles);
        if(not get_ok){ Log("SelectMuonSample Tool: Error retrieving MCParticles,true from ANNIEEvent!",v_error,verbosity); return false; }
        get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCHits",MCHits);
        if(not get_ok){ Log("SelectMuonSample Tool: Error retrieving MCHits,true from ANNIEEvent!",v_error,verbosity); return false; }
        get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCLAPPDHits",MCLAPPDHits);
        if(not get_ok){ Log("SelectMuonSample Tool: Error retrieving MCLAPPDHits,true from ANNIEEvent!",v_error,verbosity); return false; }
        get_ok = m_data->Stores.at("ANNIEEvent")->Get("TDCData",TDCData);
        if(not get_ok){ Log("SelectMuonSample Tool: Error retrieving TDCData,true from ANNIEEvent!",v_error,verbosity); return false; }
        get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventTime",EventTime);
        if(not get_ok){ Log("SelectMuonSample Tool: Error retrieving EventTime,true from ANNIEEvent!",v_error,verbosity); return false; }

        logmessage = "SelectMuonSample Tool: Processing truth tracks and digits for "+MCFile
                                +", MCEvent "+to_string(MCEventNum)+", MCTrigger "+to_string(MCTriggernum);
        Log(logmessage,v_debug,verbosity);
     
	//find primary muon

	std::cout <<"Start of Event!"<<std::endl;
	std::cout <<"Trigger Number: "<<MCTriggernum<<", Event Number: "<<EventNumber<<std::endl;

	MCParticle primarymuon, primaryelectron; 
        bool mufound=false;
	bool efound=false;

        if(MCParticles){
                Log("SelectMuonSample Tool: Num MCParticles = "+to_string(MCParticles->size()),v_message,verbosity);
                for(int particlei=0; particlei<MCParticles->size(); particlei++){
                        MCParticle aparticle = MCParticles->at(particlei);
			if(aparticle.GetParentPdg()!=0) continue; 
		        if (aparticle.GetFlag()!=0) continue;
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
			Log("SelectMuon Tool: No MCParticles in this event!",v_error,verbosity);
		}
		if (not mufound && not efound){
			Log("SelectMuon Tool: Neither electrons nor muons in this event!",v_error,verbosity);
			return true;
		}

	if (mufound){
	const Position neutrinovtx = primarymuon.GetStartVertex();
	vtx_x = neutrinovtx.X()-tank_center_x;
	vtx_y = neutrinovtx.Y()-tank_center_y;
	vtx_z = neutrinovtx.Z()-tank_center_z;
	const Position neutrinovtx_end = primaryelectron.GetStopVertex();
	vtx_stop_x = neutrinovtx_end.X()-tank_center_x;
	vtx_stop_y = neutrinovtx_end.Y()-tank_center_y;
	vtx_stop_z = neutrinovtx_end.Z()-tank_center_z;
	vtx_dir_x = vtx_stop_x-vtx_x;
	vtx_dir_y = vtx_stop_y-vtx_y;
	vtx_dir_z = vtx_stop_z-vtx_z;

	const Direction muondirection = primarymuon.GetStartDirection();
	dir_x = muondirection.X();
	dir_y = muondirection.Y();
	dir_z = muondirection.Z();
/*	const Direction stopmuon = primarymuon.GetStopDirection();
	dir_stop_x = stopmuon.X();
	dir_stop_y = stopmuon.Y();
	dir_stop_z = stopmuon.Z();
*/
	muonenergy = primarymuon.GetStartEnergy();
	logmessage = "SelectMuonSample Tool: Interaction Vertex is at ("+to_string(neutrinovtx.X())
                +", "+to_string(neutrinovtx.Y())+", "+to_string(neutrinovtx.Z())+")\n"
                +"Primary muon has energy "+to_string(muonenergy)+"GeV and direction ("
                +to_string(muondirection.X())+", "+to_string(muondirection.Y())+", "+to_string(muondirection.Z())+")\n"
		+"Pimary muon has stop vertex ("+to_string(vtx_stop_x)+", "+to_string(vtx_stop_y)+", "+to_string(vtx_stop_z)+" )";
        Log(logmessage,v_debug,verbosity);
		h_muon_energy->Fill(muonenergy);
		h_muon_vertex->Fill(neutrinovtx.X(),neutrinovtx.Y(),neutrinovtx.Z());
		h_muon_dir->Fill(muondirection.X(),muondirection.Y(),muondirection.Z());
		mu_energy = muonenergy;
		mu_pos[0]=vtx_x;//neutrinovtx.X()-tank_center_x;
		mu_pos[1]=vtx_y;//neutrinovtx.Y()-tank_center_y;
		mu_pos[2]=vtx_z;//neutrinovtx.Z()-tank_center_z;
		mu_stop_pos[0]=vtx_stop_x;
		mu_stop_pos[1]=vtx_stop_y;
		mu_stop_pos[2]=vtx_stop_z;
		mu_start_stop[0]=vtx_dir_x;
		mu_start_stop[1]=vtx_dir_y;
		mu_start_stop[2]=vtx_dir_z;
		mu_dir[0]=muondirection.X();
		mu_dir[1]=muondirection.Y();
		mu_dir[2]=muondirection.Z();
	}
	else if (efound){
	const Position neutrinovtx = primaryelectron.GetStartVertex();
	vtx_x = neutrinovtx.X()-tank_center_x;
	vtx_y = neutrinovtx.Y()-tank_center_y;
	vtx_z = neutrinovtx.Z()-tank_center_z;

	const Position neutrinovtx_end = primaryelectron.GetStopVertex();
	vtx_stop_x = neutrinovtx_end.X()-tank_center_x;
	vtx_stop_y = neutrinovtx_end.Y()-tank_center_y;
	vtx_stop_z = neutrinovtx_end.Z()-tank_center_z;
	vtx_dir_x = vtx_stop_x-vtx_x;
	vtx_dir_y = vtx_stop_y-vtx_y;
	vtx_dir_z = vtx_stop_z-vtx_z;

	const Direction electrondirection = primaryelectron.GetStartDirection();
	dir_x = electrondirection.X();
	dir_y = electrondirection.Y();
	dir_z = electrondirection.Z();
	electronenergy = primaryelectron.GetStartEnergy();
	logmessage = "SelectMuonSample Tool: Electron Interaction Vertex is at ("+to_string(neutrinovtx.X())
                +", "+to_string(neutrinovtx.Y())+", "+to_string(neutrinovtx.Z())+")\n"
                +"Primary electron has energy "+to_string(electronenergy)+"GeV and direction ("
                +to_string(electrondirection.X())+", "+to_string(electrondirection.Y())+", "+to_string(electrondirection.Z())+")\n"
		+"Pimary electron has stop vertex ("+to_string(vtx_stop_x)+", "+to_string(vtx_stop_y)+", "+to_string(vtx_stop_z)+" )";
        Log(logmessage,v_debug,verbosity);
		h_electron_energy->Fill(electronenergy);
		h_electron_vertex->Fill(neutrinovtx.X(),neutrinovtx.Y(),neutrinovtx.Z());
		h_electron_dir->Fill(electrondirection.X(),electrondirection.Y(),electrondirection.Z());
		mu_energy = electronenergy;
		mu_pos[0]=vtx_x;//neutrinovtx.X();
		mu_pos[1]=vtx_y;//neutrinovtx.Y();
		mu_pos[2]=vtx_z;//neutrinovtx.Z();
		mu_stop_pos[0]=vtx_stop_x;
		mu_stop_pos[1]=vtx_stop_y;
		mu_stop_pos[2]=vtx_stop_z;
		mu_start_stop[0]=vtx_dir_x;
		mu_start_stop[1]=vtx_dir_y;
		mu_start_stop[2]=vtx_dir_z;
		mu_dir[0]=electrondirection.X();
		mu_dir[1]=electrondirection.Y();
		mu_dir[2]=electrondirection.Z();
	}

	if (MCHits){


		mu_pmt_single.clear();
		mu_pe_single.clear();
		mu_distance_single.clear();
		mu_angle_single.clear();
		time_pmt.clear();
	        mu_pe = 0;
		Log("SelectMuonSample Tool: Num PMT Digits = "+to_string(MCHits->size()),v_message,verbosity);
		for(std::pair<ChannelKey,std::vector<Hit>>&& apair : *MCHits){

                	ChannelKey chankey = apair.first;
         //       	if(chankey.GetSubDetectorType()==subdetector::ADC){
                        	std::vector<Hit>& hits = apair.second;
                        	for(Hit& ahit : hits){
		
					//std::cout <<"mc hits found: ";
					int det_id = ahit.GetTubeId();
					mu_pmt_single.push_back(det_id);
					//std::cout <<det_id<<", ";
			
					double x_pmt = pmts_x.at(det_id);
					double y_pmt = pmts_y.at(det_id);
					double z_pmt = pmts_z.at(det_id);
					double x_dir_pmt = x_pmt - vtx_x;
					double y_dir_pmt = y_pmt - vtx_y;
					double z_dir_pmt = z_pmt - vtx_z;

					double distance_pmt = sqrt(pow(x_pmt-vtx_x,2)+pow(y_pmt-vtx_y,2)+pow(z_pmt-vtx_z,2));
					//std::cout <<distance_pmt<<", ";
					mu_distance_single.push_back(distance_pmt);
					//double angle_pmt = acos((x_dir_pmt*dir_x+y_dir_pmt*dir_y+z_dir_pmt*dir_z)/(sqrt(pow(x_dir_pmt,2)+pow(y_dir_pmt,2)+pow(z_dir_pmt,2))*sqrt(pow(dir_x,2)+pow(dir_y,2)+pow(dir_z,2))));
					double angle_pmt = acos((x_dir_pmt*dir_x+y_dir_pmt*dir_y+z_dir_pmt*dir_z)/distance_pmt);
					//std::cout <<"angle pmt: "<<angle_pmt<<std::endl;	
					mu_angle_single.push_back(angle_pmt);
					double charge_pmt = ahit.GetCharge();
					//std::cout <<charge_pmt<<std::endl;
					mu_pe+=charge_pmt;
					mu_pe_single.push_back(charge_pmt);
					//mu_pe_single.push_back(charge_pmt);
					double t_pmt = ahit.GetTime();
					time_pmt.push_back(t_pmt);
					if (mufound) h_muon_pe->Fill(charge_pmt);
					else if (efound) h_electron_pe->Fill(charge_pmt);
					//fill arrays for 3D histograms
					if (mufound){
						//std::cout <<"muonenergy: "<<muonenergy<<", bin: "<<int(floor((muonenergy-100.)/energy_step))<<std::endl;
						if (muonenergy>=2000.) {/*std::cout <<"muon energy > 2000 MeV!"<<std::endl;*/ continue;}
						//std::cout <<"distance: "<<distance_pmt<<", bin: "<<int(floor((distance_pmt/distance_step)))<<std::endl;
						//std::cout <<"angle: "<<angle_pmt<<", bin: "<<int(floor(angle_pmt/angle_step))<<std::endl;
						pe_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=charge_pmt;
						n_pe_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
						times_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=t_pmt;
						//n_times_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=1;
					}
					
					if (efound){
						if (electronenergy>=2000.) {/*std::cout <<"electron energy > 2000 MeV!"<<std::endl;*/ continue;}
						
						pe_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=charge_pmt;
						n_pe_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
						times_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=t_pmt;
		//				n_times_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
					}
					//mu_pe=charge_pmt;
				}
		//	}
		}
		if (!MCLAPPDHits) tree_beam->Fill();			
		} else {
			cout <<"No MCHits"<<endl;
		}

	if (MCLAPPDHits){

		mu_distance_lappd.clear();
		mu_angle_lappd.clear();
		mu_lappd.clear();
		time_lappd.clear();
		mu_pe_lappd=0;

		Log("SelectMuonTool: Num LAPPD Digits = "+to_string(MCLAPPDHits->size()),v_message,verbosity);

		for(std::pair<ChannelKey,std::vector<LAPPDHit>>&& apair : *MCLAPPDHits){
                	ChannelKey chankey = apair.first;

				std::vector<LAPPDHit>& hits = apair.second;
                        	for(LAPPDHit& ahit : hits){
					double lappd_charge = ahit.GetCharge();
					lappd_charge = 1.0; 				//for now just use 1 hit for the LAPPD, since the digitization is not implemented yet
					mu_pe_lappd+=lappd_charge;
					if (mufound) h_muon_lappd_pe->Fill(lappd_charge);
					else if (efound) h_electron_lappd_pe->Fill(lappd_charge);
					mu_lappd.push_back(ahit.GetTubeId());
					std::vector<double> temp_pos = ahit.GetPosition();
					double x_lappd = temp_pos.at(0)-tank_center_x;
					double y_lappd = temp_pos.at(1)-tank_center_y;
					double z_lappd = temp_pos.at(2)-tank_center_z;
					double x_dir_lappd = x_lappd - vtx_x;
					double y_dir_lappd = y_lappd - vtx_y;
					double z_dir_lappd = z_lappd - vtx_z;
					double distance_lappd = sqrt(pow(x_lappd-vtx_x,2)+pow(y_lappd-vtx_y,2)+pow(z_lappd-vtx_z,2));
					//double angle_lappd = acos((x_dir_lappd*dir_x+y_dir_lappd*dir_y+z_dir_lappd*dir_z)/(sqrt(pow(x_dir_lappd,2)+pow(y_dir_lappd,2)+pow(z_dir_lappd,2))*sqrt(pow(dir_x,2)+pow(dir_y,2)+pow(dir_z,2))));
					double angle_lappd = acos((x_dir_lappd*dir_x+y_dir_lappd*dir_y+z_dir_lappd*dir_z)/distance_lappd);
		//			std::cout <<"angle lappd: "<<angle_lappd<<std::endl;	
					double t_lappd = ahit.GetTime();
					mu_distance_lappd.push_back(distance_lappd);
					mu_angle_lappd.push_back(angle_lappd);
					time_lappd.push_back(t_lappd);			
					if (mufound){
						if (muonenergy>=2000.) continue;
						pe_muons_lappd[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]+=lappd_charge;
						n_pe_muons_lappd[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]++;
						times_muons_lappd[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]+=t_lappd;
						//n_times_muons[int(floor((muonenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]+=1;
					}
					
					if (efound){
						if (electronenergy>=2000.) continue;
						pe_electrons_lappd[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]+=lappd_charge;
						n_pe_electrons_lappd[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]++;
						times_electrons_lappd[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_lappd/distance_step))][int(floor(angle_lappd/angle_step))]+=t_lappd;
		//				n_times_electrons[int(floor((electronenergy-100.)/energy_step))][int(floor(distance_pmt/distance_step))][int(floor(angle_pmt/angle_step))]++;
					}
					//mu_pe=charge_pmt;
				
				}
		}

	tree_beam->Fill();

	}

	else cout <<"No MCLAPPDHits"<<endl;

	std::cout <<"SelectMuonSample tool: Execution finished"<<std::endl;

	return true;

}


bool SelectMuonSample::Finalise(){

  std::cout <<"SelectMuonSample tool: Finalisation started"<<std::endl;
 

  //fill 3D pdfs:
  //

  for (int i_energy=0;i_energy<21;i_energy++){
	for (int i_distance=0;i_distance<100;i_distance++){
		for (int i_angle=0;i_angle<100;i_angle++){
			if (n_pe_electrons[i_energy][i_distance][i_angle]==0) n_pe_electrons[i_energy][i_distance][i_angle]=1;
			if (n_pe_muons[i_energy][i_distance][i_angle]==0) n_pe_muons[i_energy][i_distance][i_angle]=1;
			electrons_pe->SetBinContent(i_distance+1,i_angle+1,i_energy+1,pe_electrons[i_energy][i_distance][i_angle]/n_pe_electrons[i_energy][i_distance][i_angle]);
			electrons_time->SetBinContent(i_distance+1,i_angle+1,i_energy+1,times_electrons[i_energy][i_distance][i_angle]/n_pe_electrons[i_energy][i_distance][i_angle]);
			muons_pe->SetBinContent(i_distance+1,i_angle+1,i_energy+1,pe_muons[i_energy][i_distance][i_angle]/n_pe_muons[i_energy][i_distance][i_angle]);
			muons_time->SetBinContent(i_distance+1,i_angle+1,i_energy+1,times_muons[i_energy][i_distance][i_angle]/n_pe_muons[i_energy][i_distance][i_angle]);


			if (n_pe_electrons_lappd[i_energy][i_distance][i_angle]==0) n_pe_electrons_lappd[i_energy][i_distance][i_angle]=1;
			if (n_pe_muons_lappd[i_energy][i_distance][i_angle]==0) n_pe_muons_lappd[i_energy][i_distance][i_angle]=1;
			electrons_pe_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,pe_electrons_lappd[i_energy][i_distance][i_angle]/n_pe_electrons_lappd[i_energy][i_distance][i_angle]);
			electrons_time_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,times_electrons_lappd[i_energy][i_distance][i_angle]/n_pe_electrons_lappd[i_energy][i_distance][i_angle]);
			muons_pe_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,pe_muons_lappd[i_energy][i_distance][i_angle]/n_pe_muons_lappd[i_energy][i_distance][i_angle]);
			muons_time_lappd->SetBinContent(i_distance+1,i_angle+1,i_energy+1,times_muons_lappd[i_energy][i_distance][i_angle]/n_pe_muons_lappd[i_energy][i_distance][i_angle]);

		}
	}
  }

  out_file_pdf->cd();
  electrons_time->Write();
  electrons_pe->Write();
  muons_time->Write();
  muons_pe->Write();
  electrons_time_lappd->Write();
  electrons_pe_lappd->Write();
  muons_time_lappd->Write();
  muons_pe_lappd->Write();

  out_file->cd();

  h_muon_energy->Write();
  h_muon_pe->Write();
  h_muon_lappd_pe->Write();
  h_muon_vertex->Write();
  h_muon_dir->Write();
  h_electron_energy->Write();
  h_electron_pe->Write();
  h_electron_lappd_pe->Write();
  h_electron_vertex->Write();
  h_electron_dir->Write();	
  tree_beam->Write();

  out_file->Close();
  out_file_pdf->Close(); 

  std::cout <<"SelectMuonSample tool: Finalisation complete"<<std::endl;

  return true;
}
