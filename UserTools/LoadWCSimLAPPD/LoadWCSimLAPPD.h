/* vim:set noexpandtab tabstop=4 wrap */
#ifndef LoadWCSimLAPPD_H
#define LoadWCSimLAPPD_H

#include "Tool.h"



class LoadWCSimLAPPD: public Tool {
	
	public:
	
	LoadWCSimLAPPD();
	bool Initialise(std::string configfile,DataModel &data);
	bool Execute();
	bool Finalise();
	
	private:
	// config file variables
	int verbose=1;
};


#endif
