#include "DigitBuilder.h"


static DigitBuilder* fgDigitBuilder = 0;
DigitBuilder* DigitBuilder::Instance()
{
  if( !fgDigitBuilder ){
    fgDigitBuilder = new DigitBuilder();
  }

  return fgDigitBuilder;
}

DigitBuilder::DigitBuilder():Tool(){}
DigitBuilder::~DigitBuilder() {
}

bool DigitBuilder::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(verbosity) cout<<"Initializing Tool DigitBuilder"<<endl;
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  
  ///////////////////// Defaults for Config ///////////////
  fPhotodetectorConfiguration = "All";
  fParametricModel = 0;
  fIsMC = 1;
  fDigitChargeThr = 10;


  /// Get the Tool configuration variables
  m_variables.Get("verbosity",verbosity);
  m_variables.Get("isMC",fIsMC);
  m_variables.Get("ParametricModel", fParametricModel);
  m_variables.Get("PhotoDetectorConfiguration", fPhotodetectorConfiguration);
  m_variables.Get("xshift", xshift);
  m_variables.Get("yshift", yshift);
  m_variables.Get("zshift", zshift);
  m_variables.Get("LAPPDIDFile", fLAPPDIDFile);
  m_variables.Get("DigitChargeThr",fDigitChargeThr);


  /// Construct the other objects we'll be setting at event level,
  fDigitList = new std::vector<RecoDigit>;


  // Make the RecoDigit Store if it doesn't exist
  int recoeventexists = m_data->Stores.count("RecoEvent");
  if(recoeventexists==0) m_data->Stores["RecoEvent"] = new BoostStore(false,2);
  
  // Some hard-coded values of old WCSim LAPPDIDs are in this Tool
  // I would recommend moving away from the use of WCSim IDs if possible as they are liable to change
  // but for tools that need them, in the LoadWCSim tool I put a map of WCSim TubeId to channelkey
  m_data->CStore.Get("detectorkey_to_lappdid",detectorkey_to_lappdid);
  m_data->CStore.Get("channelkey_to_pmtid",channelkey_to_pmtid);

  //Read the LAPPDID file, if given
  if(fLAPPDIDFile!="none"){
    if(verbosity>2) std::cout << "Loading digits from LAPPD IDs in file " << fLAPPDIDFile << std::endl;
    this->ReadLAPPDIDFile();
  } else {
    if(verbosity>2) std::cout << "Loading digits from all LAPPDs" << std::endl;
  }
  return true;
}

bool DigitBuilder::Execute(){
	Log("===========================================================================================",v_debug,verbosity);
	
	/// Reset everything
	this->Reset();
	
	// see if "ANNIEEvent" exists
	auto get_annieevent = m_data->Stores.count("ANNIEEvent");
	if(!get_annieevent){
		Log("DigitBuilder Tool: No ANNIEEvent store!",v_error,verbosity); 
		return false;
	};
	
	/// see if "RecoEvent" exists.  If not, make it
 	auto get_recoevent = m_data->Stores.count("RecoEvent");
 	if(!get_recoevent){
  		Log("DigitBuilder Tool: No RecoEvent store!",v_error,verbosity); 
  		return false;
	};
	
    /// Retrieve necessary info from ANNIEEvent
	auto get_geometry= m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",fGeometry);
	if(!get_geometry){
		Log("DigitBuilder Tool: Error retrieving Geometry from ANNIEEvent!",v_error,verbosity); 
		return false; 
	}
    
    if (fIsMC){
	  auto get_mchits = m_data->Stores.at("ANNIEEvent")->Get("MCHits",fMCPMTHits);
	  if(!get_mchits){ 
	  	Log("DigitBuilder Tool: Error retrieving MCHits from ANNIEEvent!",v_error,verbosity); 
	  	return false; 
	  }
    } else {
      Log("DigitBuilder Tool: ERROR; Currently no non-MC hits boost store!",v_error,verbosity);
      return false;
    }

  /// Build RecoDigit
  if (fIsMC){
    this->BuildMCRecoDigit();
  } else {
    Log("DigitBuilder Tool: ERROR; Currently no non-MC hits loading implemented!",v_error,verbosity);
  }
	
  /// Hit info. to RecoEvent
  this->PushRecoDigits(true); 
  return true;
}

bool DigitBuilder::Finalise(){
  delete fDigitList; fDigitList = 0;
  if(verbosity>0) cout<<"DigitBuilder exitting"<<endl;
  return true;
}

bool DigitBuilder::BuildMCRecoDigit() {
	
	if(fPhotodetectorConfiguration == "PMT_only") {
		this->BuildMCPMTRecoDigit();
		return true;
	}
	if(fPhotodetectorConfiguration == "LAPPD_only") {
		this->BuildMCLAPPDRecoDigit();
		return true;
	}
	if(fPhotodetectorConfiguration == "All") {
		this->BuildMCPMTRecoDigit();
	  this->BuildMCLAPPDRecoDigit();
	  return true;
	}
	else {
	  cout<<"Wrong PhotoDetector Configuration! Allowed configurations: PMT_only, LAPPD_only, All"<<endl;
	  return false;
	}
	
}

bool DigitBuilder::BuildMCPMTRecoDigit() {
	
	Log("DigitBuilder Tool: Build PMT reconstructed digits",v_message,verbosity);
	/// now move to digit retrieval
	int region = -999;
	double calT;
	double calQ = 0.;
	int digitType = -999;
	Detector* det=nullptr;
	Position  pos_sim, pos_reco;
	/// MCHits is a std::map<unsigned long,std::vector<Hit>>
	if(fMCPMTHits){
		Log("DigitBuilder Tool: Num PMT Digits = "+to_string(fMCPMTHits->size()),v_message, verbosity);
		/// iterate over the map of sensors with a measurement
		for(std::pair<unsigned long,std::vector<MCHit>>&& apair : *fMCPMTHits){
			unsigned long chankey = apair.first;
			// the channel key is a unique identifier of this signal input channel
			det = fGeometry->ChannelToDetector(chankey);
			int PMTId = channelkey_to_pmtid.at(chankey);  //PMTID In WCSim
			if(det==nullptr){
				Log("DigitBuilder Tool: Detector not found! ",v_message,verbosity);
				continue;
			}
			
			// convert the WCSim coordinates to the ANNIEreco coordinates
			// convert the unit from m to cm
			pos_sim = det->GetDetectorPosition();
			pos_sim.UnitToCentimeter();
			pos_reco.SetX(pos_sim.X()+xshift);
			pos_reco.SetY(pos_sim.Y()+yshift);
			pos_reco.SetZ(pos_sim.Z()+zshift);
	
			if(det->GetDetectorElement()=="Tank"){
				std::vector<MCHit>& hits = apair.second;
        if(fParametricModel){
          if(verbosity>2) std::cout << "Using parametric model to build PMT hits" << std::endl;
          //We'll get all hit info and then define a time/charge for each digit
          std::vector<double> hitTimes;
          std::vector<double> hitCharges;
          for(MCHit& ahit : hits){
            if(verbosity>3){
              std::cout << "This HIT'S TIME AND CHARGE: " << ahit.GetTime() <<
                  "," << ahit.GetCharge() << std::endl;
            }
            double hitTime = ahit.GetTime()*1.0;
          	if(hitTime>-10 && hitTime<40) {
			  hitTimes.push_back(ahit.GetTime()*1.0); 
              hitCharges.push_back(ahit.GetCharge());
            }
          }
          // Do median and sum
          std::sort(hitTimes.begin(), hitTimes.end());
          size_t timesize = hitTimes.size();
          if (timesize == 0) continue;
          if (timesize % 2 == 0){
            calT = (hitTimes.at(timesize/2 - 1) + hitTimes.at(timesize/2))/2;
          } else {
            calT = hitTimes.at(timesize/2);
          }
          calQ = 0.;
          for(std::vector<double>::iterator it = hitCharges.begin(); it != hitCharges.end(); ++it){
            calQ += *it;
          }
          if (verbosity>4) { 
            std::cout << "PMT position (X<Y<Z): " << 
                    to_string(pos_reco.X()) << "," << to_string(pos_reco.Y()) <<
                    "," << to_string(pos_reco.Z()) << std::endl;
            std::cout << "PMT Charge,Time: " << to_string(calQ) << "," <<
                    to_string(calT) << std::endl;
          }
          if(calQ>fDigitChargeThr) {					//changed to 0 for cross-checks with other tools, change back later!
				    digitType = RecoDigit::PMT8inch;
				    RecoDigit recoDigit(region, pos_reco, calT, calQ, digitType, PMTId);
				    fDigitList->push_back(recoDigit); 
				  }
        } else {
			    for(MCHit& ahit : hits){
				  	//if(v_message<verbosity) ahit.Print(); // << VERY verbose
				  	// get calibrated PMT time (Use the MC time for now)
				  	calT = ahit.GetTime()*1.0; 
            calQ = ahit.GetCharge();
            if (verbosity>4) { 
              std::cout << "PMT position (X<Y<Z): " << 
                      to_string(pos_reco.X()) << "," << to_string(pos_reco.Y()) <<
                      "," << to_string(pos_reco.Z()) << std::endl;
              std::cout << "PMT Charge,Time: " << to_string(calQ) << "," <<
                      to_string(calT) << std::endl;
            }
				  	digitType = RecoDigit::PMT8inch;
				  	RecoDigit recoDigit(region, pos_reco, calT, calQ, digitType, PMTId);
				    //recoDigit.Print();
				    fDigitList->push_back(recoDigit); 
          }
			  }
      }
		} // end loop over MCHits
	} else {
		cout<<"No MCHits"<<endl;
		return false;
	}
	return true;
}

bool DigitBuilder::BuildMCLAPPDRecoDigit() {
	return true;
}

void DigitBuilder::PushRecoDigits(bool savetodisk) {
	Log("DigitBuilder Tool: Push reconstructed digits to the RecoEvent store",v_message,verbosity);
	m_data->Stores.at("RecoEvent")->Set("RecoDigit", fDigitList, savetodisk);  ///> Add digits to RecoEvent
}

void DigitBuilder::Reset() {
  // Reset 
  fDigitList->clear();
}

void DigitBuilder::ReadLAPPDIDFile() {
  std::string line;
  ifstream myfile(fLAPPDIDFile);
  if (myfile.is_open()){
    while(getline(myfile,line)){
      if(verbosity>0){
        std::cout << "DigitBuilder tool: Loading hits from LAPPD ID " << line << std::endl;
      }
      int thisID = std::atoi(line.c_str());
      fLAPPDId.push_back(thisID);
    }
  } else {
    Log("Unable to open given LAPPD ID File. Using all LAPPDs",v_error,verbosity);
  }
}
