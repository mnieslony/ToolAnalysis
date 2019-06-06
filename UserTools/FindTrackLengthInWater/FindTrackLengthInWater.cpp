#include "FindTrackLengthInWater.h"
#include <boost/filesystem.hpp>
#include "TMath.h"

FindTrackLengthInWater::FindTrackLengthInWater():Tool(){}


bool FindTrackLengthInWater::Initialise(std::string configfile, DataModel &data){

  /////////////////// Usefull header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();

  m_data= &data; //assigning transient data pointer
  /////////////////////////////////////////////////////////////////
  // get configuration variables for this tool
  m_variables.Get("verbosity",verbosity);

  // make the BoostStore to replace the csv file
  BoostStore* energystore = new BoostStore(true,0); // type is single-event binary file
  m_data->Stores.emplace("EnergyReco",energystore);

  // Retrieve variables from m_data to pass to the python scripts
  // FIXME supposedly there is a way to pass variables from the python script config
  // files to the python scripts directly, rather than using another Tool to put
  // them into a BoostStore...
  // ==========================
  // Variables to be passed to DNNFindTrackLenghInWater
  // --------------------------------------------------
std::cout<<"getting DNN variables"<<std::endl;
  std::string TrackLengthTrainingDataFile;
  std::string TrackLengthTestingDataFile;
  std::string TrackLengthCheckpointDir;
  std::string TrackLengthPredictionsFile;
  // retrieve from m_data
  get_ok = m_variables.Get("TrackLengthTrainingDataFile",TrackLengthTrainingDataFile);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No TrackLengthTrainingDataFile specified, will not be written",v_warning,verbosity);
    TrackLengthTrainingDataFile="NA";
  }
  get_ok = m_variables.Get("TrackLengthTestingDataFile",TrackLengthTestingDataFile);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No TrackLengthTestingDataFile specified, will not be written",v_warning,verbosity);
    TrackLengthTestingDataFile="NA";
  }
  get_ok = m_variables.Get("TrackLengthCheckpointDir",TrackLengthCheckpointDir);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No TrackLengthCheckpointDir specified! This must specify the location of the DNN model! Aborting...",v_error,verbosity);
    return false;
  }
  get_ok = m_variables.Get("TrackLengthPredictionsFile",TrackLengthPredictionsFile);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No TrackLengthPredictionsFile specified, will not be written",v_warning,verbosity);
    TrackLengthPredictionsFile="NA";
  }
  get_ok = m_variables.Get("MaxTotalHitsToDNN",maxhits0);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No MaxTotalHitsToDNN specified: assuming 1100, but this MUST match the value used for DNN training!",v_warning,verbosity);
    maxhits0=1100;
  }
  get_ok = m_variables.Get("DNNTrainingEntries",trainingentries);
  if(not get_ok){
    auto errorlevel = (TrackLengthTestingDataFile!="NA") ? v_warning : v_message;
    Log("FindTrackLengthInWater Tool: No DNNTrainingEntries specified: will not split entries into test+training files",errorlevel,verbosity);
    trainingentries=-1;
  }
std::cout<<"passing to EnergyReco store"<<std::endl;
  // Pass to EnergyReco booststore
  m_data->Stores.at("EnergyReco")->Set("TrackLengthTrainingDataFile",TrackLengthTrainingDataFile);
  m_data->Stores.at("EnergyReco")->Set("TrackLengthTestingDataFile",TrackLengthTestingDataFile);
  m_data->Stores.at("EnergyReco")->Set("TrackLengthPredictionsFile",TrackLengthPredictionsFile);
  m_data->Stores.at("EnergyReco")->Set("TrackLengthCheckpointDir",TrackLengthCheckpointDir);
  
  // variables to be passed to BDTMuonEnergyReco
  // -------------------------------------------
std::cout<<"getting BDTMuon variables"<<std::endl;
  std::string MuonEnergyTrainingDataFile;
  double BDT_NuE_threshold;
  std::string BDTMuonModelFile;
  std::string MuonEnergyTestingDataFile;
  std::string MuonEnergyPredictionsFile;
  // retrieve from m_data
  get_ok = m_variables.Get("MuonEnergyTrainingDataFile",MuonEnergyTrainingDataFile);  // input file
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No MuonEnergyTrainingDataFile specified, will not be generated",v_warning,verbosity);
    MuonEnergyTrainingDataFile="NA";
  }
  get_ok = m_variables.Get("MuonEnergyTestingDataFile",MuonEnergyTestingDataFile);    // input file
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No MuonEnergyTestingDataFile specified, will not be generated",v_warning,verbosity);
    MuonEnergyTestingDataFile="NA";
  }
  get_ok = m_variables.Get("MuonEnergyPredictionsFile",MuonEnergyPredictionsFile);    // output file
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No MuonEnergyPredictionsFile specified, will not be generated",v_warning,verbosity);
    MuonEnergyPredictionsFile="NA";
  }
  get_ok = m_variables.Get("BDTMuonModelFile",BDTMuonModelFile);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No BDTMuonModelFile specified! Cannot train or predict Muon Energy!",v_warning,verbosity);
    BDTMuonModelFile="NA";
  }
  get_ok = m_variables.Get("BDT_NuE_threshold",BDT_NuE_threshold);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No upper threshold on Neutrino energy for BDT training specified, assuming 2GeV",v_warning,verbosity);
    BDT_NuE_threshold=2.;
  }
  // Pass to EnergyReco booststore
std::cout<<"passing to store"<<std::endl;
  m_data->Stores.at("EnergyReco")->Set("MuonEnergyTrainingDataFile",MuonEnergyTrainingDataFile);
  m_data->Stores.at("EnergyReco")->Set("MuonEnergyTestingDataFile",MuonEnergyTestingDataFile);
  m_data->Stores.at("EnergyReco")->Set("MuonEnergyPredictionsFile",MuonEnergyPredictionsFile);
  m_data->Stores.at("EnergyReco")->Set("BDTMuonModelFile",BDTMuonModelFile);
  m_data->Stores.at("EnergyReco")->Set("BDT_NuE_threshold",BDT_NuE_threshold);

  // variables to be passed to BDTNeutrinoEnergyReco
  // -----------------------------------------------
std::cout<<"getting BDTNeutrino variables"<<std::endl;
  std::string NeutrinoEnergyTrainingDataFile;
  std::string BDTNeutrinoModelFile;
  std::string NeutrinoEnergyTestingDataFile;
  std::string NeutrinoEnergyPredictionsFile;
  // retrieve from m_data
  get_ok = m_variables.Get("NeutrinoEnergyTrainingDataFile",NeutrinoEnergyTrainingDataFile);  // input file
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No NeutrinoEnergyTrainingDataFile specified, will not be generated",v_warning,verbosity);
    NeutrinoEnergyTrainingDataFile="NA";
  }
  get_ok = m_variables.Get("NeutrinoEnergyTestingDataFile",NeutrinoEnergyTestingDataFile);    // input file
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No NeutrinoEnergyTestingDataFile specified, will not be generated",v_warning,verbosity);
    NeutrinoEnergyTestingDataFile="NA";
  }
  get_ok = m_variables.Get("NeutrinoEnergyPredictionsFile",NeutrinoEnergyPredictionsFile);    // output file
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No NeutrinoEnergyPredictionsFile specified, will not be generated",v_warning,verbosity);
    NeutrinoEnergyPredictionsFile="NA";
  }
  get_ok = m_variables.Get("BDTNeutrinoModelFile",BDTNeutrinoModelFile);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No BDTNeutrinoModelFile specified! Cannot train or predict Neutrino Energy!",v_warning,verbosity);
    BDTNeutrinoModelFile="NA";
  }
  // Pass to EnergyReco booststore
std::cout<<"passing to store"<<std::endl;
  m_data->Stores.at("EnergyReco")->Set("NeutrinoEnergyTrainingDataFile",NeutrinoEnergyTrainingDataFile);
  m_data->Stores.at("EnergyReco")->Set("NeutrinoEnergyTestingDataFile",NeutrinoEnergyTestingDataFile);
  m_data->Stores.at("EnergyReco")->Set("NeutrinoEnergyPredictionsFile",NeutrinoEnergyPredictionsFile);
  m_data->Stores.at("EnergyReco")->Set("BDTNeutrinoModelFile",BDTNeutrinoModelFile);
  
  // ==================================
  
      std::cout<<"max number of hits: "<<maxhits0<<std::endl;  
      if(maxhits0>1100){ 
        std::cerr<<" Please change the dim of double lambda_vec[1100]={0.}; double digitt[1100]={0.}; from 1100 to max number of hits"<<std::endl; 
      }
      if(TrackLengthTrainingDataFile!="NA"){
        // if we're not doing training, don't write the csv files
        if(trainingentries>1){
std::cout<<"writing DNN training and testing files"<<std::endl;
          // if we're splitting the run up into training and testing samples, we need to generate multiple output csvs
//          boost::filesystem::path P(TrackLengthTrainingDataFile.c_str());
//          boost::filesystem::path newPath1 = P.parent_path() / boost::filesystem::path(P.stem().string() + "_1.csv");
//          boost::filesystem::path newPath2 = P.parent_path() / boost::filesystem::path(P.stem().string() + "_2.csv");
          tracklengthtrainingfiles.push_back(TrackLengthTrainingDataFile);
          tracklengthtrainingfiles.push_back(TrackLengthTestingDataFile);
        } else {
          // otherwise just one csv
          tracklengthtrainingfiles.push_back(TrackLengthTrainingDataFile);
std::cout<<"writing DNN training file: "<<TrackLengthTrainingDataFile<<std::endl;
        }
        // loop over the csv's we're creating and write header row
        for(std::string apath : tracklengthtrainingfiles){
std::cout<<"writing header for file "<<apath<<std::endl;
          csvfile.open(apath,std::fstream::out);
          if(!csvfile.is_open()){
           Log("FindTrackLengthInWater Tool: Failed to open "+apath+" for writing headers",v_error,verbosity);
          }
          for (int i=0; i<maxhits0;++i){
             csvfile<<"l_"<<i<<",";
          }
          for (int i=0; i<maxhits0;++i){
             csvfile<<"T_"<<i<<",";
          }
          csvfile<<"lambda_max"<<",";  //first estimation of track length(using photons projection on track)
          csvfile<<"totalPMTs"<<",";   // number of PMT hits, not number of pmts.
          csvfile<<"totalLAPPDs"<<","; // number of LAPPD hits... 
          csvfile<<"lambda_max"<<",";  // ... again...?
          csvfile<<"TrueTrackLengthInWater"<<",";
          csvfile<<"neutrinoE"<<",";
          csvfile<<"trueKE"<<",";      // of the primary muon
          csvfile<<"diffDirAbs"<<",";
          csvfile<<"TrueTrackLengthInMrd"<<",";
          csvfile<<"recoDWallR"<<",";
          csvfile<<"recoDWallZ"<<",";
          csvfile<<"dirX"<<",";        // of the reconstructed muon
          csvfile<<"dirY"<<",";
          csvfile<<"dirZ"<<",";
          csvfile<<"vtxX"<<",";        // of the reconstructed event
          csvfile<<"vtxY"<<",";
          csvfile<<"vtxZ";
          csvfile<<"recoVtxFOM"<<",";
          csvfile<<"recoStatus"<<",";
          csvfile<<"deltaVtxR"<<",";
          csvfile<<"deltaAngle"<<",";
          csvfile<<'\n';
          // }
          csvfile.close();
        }
        csvfile.open(tracklengthtrainingfiles.front(),std::fstream::app);  // open (first) file for writing events
      }
std::cout<<"getting anniegeom"<<std::endl;
  
  get_ok = m_data->Stores.at("ANNIEEvent")->Header->Get("AnnieGeometry",anniegeom);
  if(not get_ok){
    Log("FindTrackLengthInWater Tool: No Geometry in ANNIEEvent!",v_error,verbosity);
    return false;
  }
  tank_radius = anniegeom->GetTankRadius()*100.;
  tank_halfheight = anniegeom->GetTankHalfheight()*100.;
  
std::cout<<"done initializing"<<std::endl;
  return true;
}

bool FindTrackLengthInWater::Execute(){
   
   // read the input hit and reconstruction info
   // ==========================================
   Int_t recoStatus;
   double vtxX,vtxY,vtxZ,dirX,dirY,dirZ,TrueTrackLengthInMrd,TrueTrackLengthInWater,TrueNeutrinoEnergy,trueEnergy, recoVtxFOM, deltaVtxR, deltaAngle;
   std::string TrueInteractionType;
   std::vector<double> digitX; std::vector<double> digitY;  std::vector<double> digitZ;
   std::vector<double> digitT;
   std::map<unsigned long,std::vector<MCHit>>* MCHits=nullptr;
   std::map<unsigned long,std::vector<MCLAPPDHit>>* MCLAPPDHits=nullptr;
   
   // Get hits from the ANNIEEvent
   get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCHits", MCHits);           // FIXME use 'Hits' to support data
   if(not get_ok){
      Log("FindTrackLengthInWater Tool: Failed to retrieve the MCHits!",v_error,verbosity);
      return false;
   }
   get_ok = m_data->Stores.at("ANNIEEvent")->Get("MCLAPPDHits", MCLAPPDHits); // FIXME as above
   if(not get_ok){
      Log("FindTrackLengthInWater Tool: Failed to retrieve the MCLAPPDHits!",v_error,verbosity);
      return false;
   }
   
   // Get reconstructed vertex from the RecoEvent
   RecoVertex theExtendedVertex;
   get_ok = m_data->Stores.at("RecoEvent")->Get("ExtendedVertex", theExtendedVertex);
   if(not get_ok){
   	Log("FindTrackLengthInWater Tool: Failed to retrieve the ExtendedVertex from RecoEvent Store!",v_error,verbosity);
   	return false;
   }
   // get the reconstructed vertex and direction
   recoStatus = theExtendedVertex.GetStatus();
   recoVtxFOM = theExtendedVertex.GetFOM();
   vtxX = theExtendedVertex.GetPosition().X();
   vtxY = theExtendedVertex.GetPosition().Y();
   vtxZ = theExtendedVertex.GetPosition().Z();
   dirX = theExtendedVertex.GetDirection().X();
   dirY = theExtendedVertex.GetDirection().Y();
   dirZ = theExtendedVertex.GetDirection().Z();
   
   /*
   // Get neutrino info from GenieEvent TODO FIXME
   get_ok = m_data->Stores.at("GenieEvent")->Get("TrueNeutrinoEnergy", TrueNeutrinoEnergy);
   if(not get_ok){
   	Log("FindTrackLengthInWater Tool: Failed to retrieve TrueNeutrinoEnergy!",v_error,verbosity);
   	return false;
   }
   get_ok = m_data->Stores.at("GenieEvent")->Get("TrueInteractionType", TrueInteractionType);
   if(not get_ok){
   	Log("FindTrackLengthInWater Tool: Failed to retrieve the TrueInteractionType!",v_error,verbosity);
   	return false;
   }
   */
   // XXX FIXME PLACEHOLDERS UNTIL WE HAVE GENIE INFO XXX
   TrueNeutrinoEnergy =1.;
   TrueInteractionType = "QES - Weak[CC]";
   
   
   // Get the primary muon information
   // ================================
   // find muon
   int PrimaryMuonIndex;
   get_ok = m_data->Stores.at("ANNIEEvent")->Get("PrimaryMuonIndex",PrimaryMuonIndex);
   std::vector<MCParticle>* MCParticles=nullptr;
   get_ok &= m_data->Stores.at("ANNIEEvent")->Get("MCParticles", MCParticles);
   MCParticle* primarymuon=nullptr;
   if((get_ok==0) || (PrimaryMuonIndex<0)){
     Log("FindTrackLengthInWater Tool: No PrimaryMuonIndex in ANNIEEvent",v_error,verbosity); // FIXME for data?
     return false;
   } else {
     primarymuon = &(MCParticles->at(PrimaryMuonIndex));
   }
   
   // Get info
   TrueTrackLengthInWater = primarymuon->GetTrackLengthInTank();
   TrueTrackLengthInMrd = primarymuon->GetTrackLengthInMrd();
   trueEnergy = primarymuon->GetStartEnergy();
   deltaVtxR = 100.*(theExtendedVertex.GetPosition()-primarymuon->GetStartVertex()).Mag();
  double cosphi = primarymuon->GetStartDirection().X()*theExtendedVertex.GetDirection().X()+
                primarymuon->GetStartDirection().Y()*theExtendedVertex.GetDirection().Y()+
                primarymuon->GetStartDirection().Z()*theExtendedVertex.GetDirection().Z();
  double phi = TMath::ACos(cosphi); // radians
  deltaAngle = phi*TMath::RadToDeg();
   //deltaAngle = (theExtendedVertex.GetDirection()-primarymuon->GetStartDirection()).Mag();
   
   // Get the PMT hit information
   // ===========================
   int totalPMTs =0; // number of PMT hits in the event
	Log("TotalLightMap Tool: Looping over PMTs with a hit",v_debug,verbosity);
	for(std::pair<const unsigned long,std::vector<MCHit>>& nextpmt : *MCHits ){
		// if it's not a tank PMT, ignore it
		Detector* thepmt = anniegeom->ChannelToDetector(nextpmt.first);
		if(thepmt->GetDetectorElement()!="Tank") continue;
		if(thepmt->GetTankLocation()=="OD") continue;  // don't utilize OD pmts?
		totalPMTs += nextpmt.second.size();
		Position PMT_position = thepmt->GetPositionInTank();
		// loop over hits on this PMT
		// ==========================
		for(MCHit& nexthit : nextpmt.second){
			double hit_time = nexthit.GetTime();
			digitT.push_back(hit_time);
		}
		digitX.resize(digitT.size(), PMT_position.X());
		digitY.resize(digitT.size(), PMT_position.Y());
		digitZ.resize(digitT.size(), PMT_position.Z());
	}
   
   // Get the LAPPD hit information
   // =============================
   int totalLAPPDs = 0; // number of LAPPD hits in the event
	Log("TotalLightMap Tool: Looping over LAPPDs with a hit",v_debug,verbosity);
	for(std::pair<const unsigned long,std::vector<MCLAPPDHit>>& nextlappd : *MCLAPPDHits ){
		// don't actually need to get the LAPPD itself; all info we need is in the hit
		totalLAPPDs += nextlappd.second.size();
		// loop over hits on this LAPPD
		// ============================
		digitX.clear(); digitY.clear(); digitZ.clear();
		for(MCLAPPDHit& nexthit : nextlappd.second){
			double hit_time = nexthit.GetTime();
			digitT.push_back(hit_time);
			std::vector<double> hitpos = nexthit.GetPosition(); // in global coords
			Position LAPPDhitpos = Position(hitpos.at(0), hitpos.at(1), hitpos.at(2));
			LAPPDhitpos -= anniegeom->GetTankCentre();
			digitX.push_back(LAPPDhitpos.X());
			digitY.push_back(LAPPDhitpos.Y());
			digitZ.push_back(LAPPDhitpos.Z());
		}
	}
   
   // Estimate the track length in the tank
   // =====================================
   uint32_t EventNumber;
   get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventNumber", EventNumber);
   std::cout<<"EventNumber: "<<EventNumber<<endl;
   //if(recoStatus == 0){ count1++;
   if(recoVtxFOM>0){ count1++;
     // XXX FIXME XXX only if Monte Carlo!
     // XXX What about for measuring error on reconstructed energy for other event toplogies? XXX
     //if((TrueInteractionType == "QES - Weak[CC]") && TrueTrackLengthInMrd>0.){ // XXX no genie, but for data?
      if(TrueTrackLengthInMrd>0.){
   	//std::cout<<"current entry: "<<EventNumber<<" with nhits: "<<digitT.size()<<std::endl;

        //calculate diff dir with (0,0,1)  
        double diffDirAbs0 = TMath::ACos(dirZ)*TMath::RadToDeg();
        //cout<<"diffDirAbs0: "<<diffDirAbs0<<endl;    
        float diffDirAbs2=diffDirAbs0/90.;
        double recoVtxR2 = vtxX*vtxX + vtxZ*vtxZ;//vtxY*vtxY;
        double recoDWallR = tank_radius-TMath::Sqrt(recoVtxR2);   // FIXME is this coordinate-system
        double recoDWallZ = tank_halfheight-TMath::Abs(vtxY);     // dependant? Is it subtracting tank origin


	// Estimate the track length
	// =========================
	double lambda_min = 10000000;  double lambda_max = -99999999.9; double lambda = 0;
	std::vector<double> lambda_vector; // booststore works better with vectors
	for(int k=0; k<digitT.size(); k++){
          //std::cout<<"k: "<<k<<", "<<digitT.at(k)<<std::endl;

	  // Estimate length as the distance between the reconstructed vertex last photon emission point
          lambda = find_lambda(vtxX,vtxY,vtxZ,dirX,dirY,dirZ,digitX.at(k),digitY.at(k),digitZ.at(k),42.);
          if( lambda <= lambda_min ){
	      lambda_min = lambda;
	  }
	  if( lambda >= lambda_max ){
	      lambda_max = lambda;
	  }
          lambda_vector.push_back(lambda);
         //m_data->Stores["ANNIEEvent"]->Set("WaterRecoTrackLength",lambda_max);
  	}
        //std::cout<<"the track length in the water tank (1st approx) is: "<<lambda_max<<std::endl;
       
       // Post-processing of variables to store
       // =====================================
       float recoDWallR2      = recoDWallR/tank_radius;
       float recoDWallZ2      = recoDWallZ/tank_halfheight;
//       float lambda_max_2     = TMath::Abs(lambda_max)/500;
//       float totalPMTs2=double(totalPMTs)/1000.;
//       float totalLAPPDs2=double(totalLAPPDs)/1000.;
       float TrueTrackLengthInWater2 = TrueTrackLengthInWater*100.;  // convert to [cm]
       float TrueTrackLengthInMrd2 = TrueTrackLengthInMrd*100.;      // convert to [cm]
       // we need to normalise the hit time and lambda vectors to fixed dimensions to match the DNN
       lambda_vector.resize(maxhits0);
       digitT.resize(maxhits0);
       
       // Write to .csv file - including variables for track length & energy reconstruction
       // ==================
       // pick which file to write to
       if(tracklengthtrainingfiles.size()){  // if any
         if((tracklengthtrainingfiles.size()>1) && (EventNumber=trainingentries)){    // once we've processed requested
           csvfile.close();                                                           // number of training entries
           csvfile.open(tracklengthtrainingfiles.back(), std::fstream::app);          // switch output to testing file
         }
         for(int i=0; i<maxhits0;++i){
            csvfile<<lambda_vector.at(i)<<",";
         }
         for(int i=0; i<maxhits0;++i){
            csvfile<<digitT.at(i)<<",";
         }
         csvfile<<lambda_max<<",";
         csvfile<<totalPMTs<<",";
         csvfile<<totalLAPPDs<<",";
         csvfile<<lambda_max<<",";
         csvfile<<TrueTrackLengthInWater2<<",";
         csvfile<<TrueNeutrinoEnergy<<",";
         csvfile<<trueEnergy<<",";
         csvfile<<diffDirAbs2<<",";
         csvfile<<TrueTrackLengthInMrd2<<",";
         csvfile<<recoDWallR2<<",";
         csvfile<<recoDWallZ2<<",";
         csvfile<<dirX<<",";
         csvfile<<dirY<<",";
         csvfile<<dirZ<<",";
         csvfile<<vtxX<<",";
         csvfile<<vtxY<<",";
         csvfile<<vtxZ<<",";
         csvfile<<recoVtxFOM<<",";
         csvfile<<recoStatus<<",";
         csvfile<<deltaVtxR<<",";
         csvfile<<deltaAngle;
         csvfile<<'\n';
       }
        
        // Put these variables in the ANNIEEvent
        // =====================================
        m_data->Stores.at("EnergyReco")->Set("lambda_vec",lambda_vector);
        m_data->Stores.at("EnergyReco")->Set("digit_ts_vec",digitT);
        m_data->Stores.at("EnergyReco")->Set("lambda_max",lambda_max);
        m_data->Stores.at("EnergyReco")->Set("num_pmt_hits",totalPMTs);
        m_data->Stores.at("EnergyReco")->Set("num_lappd_hits",totalLAPPDs);
        m_data->Stores.at("EnergyReco")->Set("lambda_max",lambda_max);
        m_data->Stores.at("EnergyReco")->Set("TrueTrackLengthInWater",TrueTrackLengthInWater);
        m_data->Stores.at("EnergyReco")->Set("trueNeuE",TrueNeutrinoEnergy);
        m_data->Stores.at("EnergyReco")->Set("trueE",trueEnergy);
        m_data->Stores.at("EnergyReco")->Set("diffDirAbs2",diffDirAbs2);
        m_data->Stores.at("EnergyReco")->Set("TrueTrackLengthInMrd2",TrueTrackLengthInMrd2);
        // FIXME naming, if training are these actually trueDWallR2?
        m_data->Stores.at("EnergyReco")->Set("recoDWallR2",recoDWallR2);
        m_data->Stores.at("EnergyReco")->Set("recoDWallZ2",recoDWallZ2);
        m_data->Stores.at("EnergyReco")->Set("dirVec",theExtendedVertex.GetDirection());
        m_data->Stores.at("EnergyReco")->Set("vtxVec",theExtendedVertex.GetPosition());
        m_data->Stores.at("EnergyReco")->Set("recoVtxFOM",recoVtxFOM);
        m_data->Stores.at("EnergyReco")->Set("recoStatus",recoStatus);
        m_data->Stores.at("EnergyReco")->Set("deltaVtxR",deltaVtxR);
        m_data->Stores.at("EnergyReco")->Set("deltaAngle",deltaAngle);
        
     }
   }

  return true;
}

//---------------------------------------------------------------------------------------------------------//
// Find Distance between Muon Vertex and Photon Production Point (lambda) 
//---------------------------------------------------------------------------------------------------------//
double FindTrackLengthInWater::find_lambda(double xmu_rec,double ymu_rec,double zmu_rec,double xrecDir,double yrecDir,double zrecDir,double x_pmtpos,double y_pmtpos,double z_pmtpos,double theta_cher)
{
     double lambda1 = 0.0;    double lambda2 = 0.0;    double length = 0.0 ;
     double xmupos_t1 = 0.0;  double ymupos_t1 = 0.0;  double zmupos_t1 = 0.0;
     double xmupos_t2 = 0.0;  double ymupos_t2 = 0.0;  double zmupos_t2 = 0.0;
     double xmupos_t = 0.0;   double ymupos_t = 0.0;   double zmupos_t = 0.0;

     double theta_muDir_track = 0.0;
     double theta_muDir_track1 = 0.0;  double theta_muDir_track2 = 0.0;
     double cos_thetacher = cos(theta_cher*TMath::DegToRad());
     double xmupos_tmin = 0.0; double ymupos_tmin = 0.0; double zmupos_tmin = 0.0;
     double xmupos_tmax = 0.0; double ymupos_tmax = 0.0; double zmupos_tmax = 0.0;
     double lambda_min = 10000000;  double lambda_max = -99999999.9;  double lambda = 0.0;

     double alpha = (xrecDir*xrecDir + yrecDir*yrecDir + zrecDir*zrecDir) * ( (xrecDir*xrecDir + yrecDir*yrecDir + zrecDir*zrecDir) - (cos_thetacher*cos_thetacher) );
     double beta = ( (-2)*(xrecDir*(x_pmtpos - xmu_rec) + yrecDir*(y_pmtpos - ymu_rec) + zrecDir*(z_pmtpos - zmu_rec) )*((xrecDir*xrecDir + yrecDir*yrecDir + zrecDir*zrecDir) - (cos_thetacher*cos_thetacher)) );
     double gamma = ( ( (xrecDir*(x_pmtpos - xmu_rec) + yrecDir*(y_pmtpos - ymu_rec) + zrecDir*(z_pmtpos - zmu_rec))*(xrecDir*(x_pmtpos - xmu_rec) + yrecDir*(y_pmtpos - ymu_rec) + zrecDir*(z_pmtpos - zmu_rec)) ) - (((x_pmtpos - xmu_rec)*(x_pmtpos - xmu_rec) + (y_pmtpos - ymu_rec)*(y_pmtpos - ymu_rec) + (z_pmtpos - zmu_rec)*(z_pmtpos - zmu_rec))*(cos_thetacher*cos_thetacher)) );


     double discriminant = ( (beta*beta) - (4*alpha*gamma) );

     lambda1 = ( (-beta + sqrt(discriminant))/(2*alpha));
     lambda2 = ( (-beta - sqrt(discriminant))/(2*alpha));

     xmupos_t1 = xmu_rec + xrecDir*lambda1;  xmupos_t2 = xmu_rec + xrecDir*lambda2;
     ymupos_t1 = ymu_rec + yrecDir*lambda1;  ymupos_t2 = ymu_rec + yrecDir*lambda2;
     zmupos_t1 = zmu_rec + zrecDir*lambda1;  zmupos_t2 = zmu_rec + zrecDir*lambda2;

     double tr1 = sqrt((x_pmtpos - xmupos_t1)*(x_pmtpos - xmupos_t1) + (y_pmtpos - ymupos_t1)*(y_pmtpos - ymupos_t1) + (z_pmtpos - zmupos_t1)*(z_pmtpos - zmupos_t1));
     double tr2 = sqrt((x_pmtpos - xmupos_t2)*(x_pmtpos - xmupos_t2) + (y_pmtpos - ymupos_t2)*(y_pmtpos - ymupos_t2) + (z_pmtpos - zmupos_t2)*(z_pmtpos - zmupos_t2));
     theta_muDir_track1 = (acos( (xrecDir*(x_pmtpos - xmupos_t1) + yrecDir*(y_pmtpos - ymupos_t1) + zrecDir*(z_pmtpos - zmupos_t1))/(tr1) )*TMath::RadToDeg());
     theta_muDir_track2 = (acos( (xrecDir*(x_pmtpos - xmupos_t2) + yrecDir*(y_pmtpos - ymupos_t2) + zrecDir*(z_pmtpos - zmupos_t2))/(tr2) )*TMath::RadToDeg());

     //---------------------------- choose lambda!! ---------------------------------
     if( theta_muDir_track1 < theta_muDir_track2 ){
       lambda = lambda1;
       xmupos_t = xmupos_t1;
       ymupos_t = ymupos_t1;
       zmupos_t = zmupos_t1;
       theta_muDir_track=theta_muDir_track1;
     }else if( theta_muDir_track2 < theta_muDir_track1 ){
       lambda = lambda2;
       xmupos_t = xmupos_t2;
       ymupos_t = ymupos_t2;
       zmupos_t = zmupos_t2;
       theta_muDir_track=theta_muDir_track2;
     }

     return lambda;
}

bool FindTrackLengthInWater::Finalise(){
  if(csvfile.is_open()) csvfile.close();
  return true;
}
