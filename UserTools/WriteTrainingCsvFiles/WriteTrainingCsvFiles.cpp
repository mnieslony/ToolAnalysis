#include "WriteTrainingCsvFiles.h"
#include <boost/filesystem.hpp>
#include "TMath.h"

WriteTrainingCsvFiles::WriteTrainingCsvFiles():Tool(){}


bool WriteTrainingCsvFiles::Initialise(std::string configfile, DataModel &data){
  
  /////////////////// Useful header ///////////////////////
  if(configfile!="")  m_variables.Initialise(configfile); //loading config file
  //m_variables.Print();
  
  m_data= &data; //assigning transient data pointer
  
  // get configuration variables for this tool
  m_variables.Get("verbosity",verbosity);
  
  // Config Variables for this tool
  // ==============================
std::cout<<"getting DNN variables"<<std::endl;
  std::string TrackLengthTrainingDataFile;
  std::string TrackLengthTestingDataFile;
  // first check if we're splitting the tchain loops into both training and testing
  get_ok = m_variables.Get("DNNTrainingEntries",trainingentries);
  if(not get_ok){
    auto errorlevel = (TrackLengthTestingDataFile!="NA") ? v_warning : v_message;
    Log("WriteTrainingCsvFiles Tool: No DNNTrainingEntries specified: no test file will be written",errorlevel,verbosity);
    trainingentries=-1;
  }
  // retrieve from m_data
  get_ok = m_variables.Get("TrackLengthTrainingDataFile",TrackLengthTrainingDataFile);
  if(not get_ok){
    Log("WriteTrainingCsvFiles Tool: No TrackLengthTrainingDataFile specified, will not be written",v_error,verbosity);
    return false;
  }
  get_ok = m_variables.Get("TrackLengthTestingDataFile",TrackLengthTestingDataFile);
  if(not get_ok){
    Log("WriteTrainingCsvFiles Tool: No TrackLengthTestingDataFile specified, will not be written",v_error,verbosity);
    return false;
  }
  get_ok = m_data->Stores.at("EnergyReco")->Get("MaxTotalHitsToDNN",maxhits0);
  if(not get_ok){
    Log("WriteTrainingCsvFiles Tool: No MaxTotalHitsToDNN in EnergyReco store!",v_error,verbosity);
    return false;
  }
  
  // Write the file header(s)
  // ========================
  if(trainingentries>1){
std::cout<<"writing DNN training and testing files"<<std::endl;
    // if we're splitting the run up into training and testing samples, we need to generate multiple output csvs
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
     Log("WriteTrainingCsvFiles Tool: Failed to open "+apath+" for writing headers",v_error,verbosity);
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
  
std::cout<<"done initializing"<<std::endl;
  return true;
}

bool WriteTrainingCsvFiles::Execute(){
  
  // Check if the event needs writing to file
  // ========================================
  uint32_t EventNumber;
  get_ok = m_data->Stores.at("ANNIEEvent")->Get("EventNumber", EventNumber);
  std::cout<<"EventNumber: "<<EventNumber<<endl;
  // Not every event is written to file: check if this one passed checks:
  uint32_t ThisEvtNum;
  m_data->Stores.at("EnergyReco")->Set("ThisEvtNum",ThisEvtNum);
  if(ThisEvtNum!=EventNumber){ return true; } // this event didn't pass checks; don't write this entry
  
  // Retrieve variables from BoostStore
  // ==================================
  m_data->Stores.at("EnergyReco")->Set("lambda_vec",lambda_vector);
  m_data->Stores.at("EnergyReco")->Set("digit_ts_vec",digitT);
  m_data->Stores.at("EnergyReco")->Set("lambda_max",lambda_max);
  m_data->Stores.at("EnergyReco")->Set("num_pmt_hits",totalPMTs);
  m_data->Stores.at("EnergyReco")->Set("num_lappd_hits",totalLAPPDs);
  m_data->Stores.at("EnergyReco")->Set("lambda_max",lambda_max);
  m_data->Stores.at("EnergyReco")->Set("TrueTrackLengthInWater",TrueTrackLengthInWater2);
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
  
  // Write to .csv file
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
  
  return true;
}

bool WriteTrainingCsvFiles::Finalise(){
  if(csvfile.is_open()) csvfile.close();
  return true;
}
