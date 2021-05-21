void create_v7(){


	ofstream chkeyfile("Chankey_WCSimID_v7.txt");


	for (int i=0; i< 132; i++){
		chkeyfile << i+332 << "    " << i+1 << std::endl;
	}
	chkeyfile.close();
}
