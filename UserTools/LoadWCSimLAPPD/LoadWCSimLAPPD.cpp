/* vim:set noexpandtab tabstop=4 wrap */

#include "LoadWCSimLAPPD.h"

LoadWCSimLAPPD::LoadWCSimLAPPD():Tool(){}

bool LoadWCSimLAPPD::Initialise(std::string configfile, DataModel &data){
	
	/////////////////// Useful header ///////////////////////
	
	if(verbose) cout<<"Initializing Tool LoadWCSimLAPPD"<<endl;
	
	if(configfile!="") m_variables.Initialise(configfile); //loading config file
	//m_variables.Print();
	
	m_data= &data; //assigning transient data pointer
	/////////////////////////////////////////////////////////////////
	
	int returnval=0; // All BoostStore 'Get' calls must have their return value checked!!!
	// 0 means the entry was not found, value is NOT set! 1 means OK. No compile time checking!!
	
	// Get the Tool configuration variables
	// ====================================
	m_variables.Get("verbose",verbose);
	//verbose=10;
	
	return true;
}


bool LoadWCSimLAPPD::Execute(){
	
	
	return true;
}


bool LoadWCSimLAPPD::Finalise(){
	
	
	return true;
}
