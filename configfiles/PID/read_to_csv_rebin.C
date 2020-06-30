void read_to_csv_rebin(){


std::string electron_filename = "electron_pdf_022619.root";
std::string muon_filename = "muon_pdf_022619_correctMuonPDF.root";
TH3F *electrons_pe, *muons_pe, *electrons_time, *muons_time;

TFile *electron_file = new TFile(electron_filename.c_str(),"READ");
TFile *muon_file = new TFile(muon_filename.c_str(),"READ");

electron_file->cd();
electrons_pe = (TH3F*) electron_file->Get("electrons_pe");
electrons_time = (TH3F*) electron_file->Get("electrons_time");
muon_file->cd();
muons_pe = (TH3F*) muon_file->Get("muons_pe");
muons_time = (TH3F*) muon_file->Get("muons_time");

const int num_E = 20;
const int num_dist = 100;
const int num_angle = 100;

//rebin to 10 energy bins and 20 distance and angle bins
double charge[num_dist/5][num_angle/5][num_E/2] = {0};
double t[num_dist/5][num_angle/5][num_E/2] = {0};
double chargeMu[num_dist/5][num_angle/5][num_E/2] = {0};
double tMu[num_dist/5][num_angle/5][num_E/2] = {0};
double dist_rebin[num_dist/5] = {0};
double angle_rebin[num_angle/5] = {0};
double energy_rebin[num_E/2] = {0};
long n[num_dist/5][num_angle/5][num_E/2] = {0};


//fill electron csv file
std::cout <<"ELECTRONS: ----------------------------------------------------"<<std::endl;
for (int i_distance = 0; i_distance < num_dist; i_distance++){
	dist_rebin[i_distance/5] += electrons_pe->GetXaxis()->GetBinCenter(i_distance);
	angle_rebin[i_distance/5] += electrons_pe->GetYaxis()->GetBinCenter(i_distance);
	if (i_distance < num_E) energy_rebin[i_distance/2] += electrons_pe->GetZaxis()->GetBinCenter(i_distance);
	for (int i_angle = 0; i_angle < num_angle; i_angle++){
		for (int i_energy = 0; i_energy < num_E; i_energy++){
			charge[i_distance/5][i_angle/5][i_energy/2] += electrons_pe->GetBinContent(i_distance+1,i_angle+1,i_energy+1);
			t[i_distance/5][i_angle/5][i_energy/2] += electrons_time->GetBinContent(i_distance+1,i_angle+1,i_energy+1);
			chargeMu[i_distance/5][i_angle/5][i_energy/2] += muons_pe->GetBinContent(i_distance+1,i_angle+1,i_energy+1);
			tMu[i_distance/5][i_angle/5][i_energy/2] += muons_time->GetBinContent(i_distance+1,i_angle+1,i_energy+1);
			n[i_distance/5][i_angle/5][i_energy/2] ++;
		}
	}
}
ofstream file_e("electrons_uniform_bins_20x_20y_10z.csv");
for (int i_distance = 0; i_distance < num_dist/5; i_distance++){
	for (int i_angle = 0; i_angle < num_angle/5; i_angle++){
		for (int i_energy = 0; i_energy < num_E/2; i_energy++){
			std::cout <<"dist: "<<dist_rebin[i_distance]/5.<<", angle: "<<angle_rebin[i_angle]/5.<<", energy: "<<energy_rebin[i_energy]/2.<<", charge: "<<charge[i_distance][i_angle][i_energy]/double(n[i_distance][i_angle][i_energy])<<", time: "<<t[i_distance][i_angle][i_energy]/double(n[i_distance][i_angle][i_energy])<<std::endl;
			file_e<<dist_rebin[i_distance]/5.<<","<<angle_rebin[i_angle]/5.<<","<<energy_rebin[i_energy]/2.<<","<<charge[i_distance][i_angle][i_energy]/double(n[i_distance][i_angle][i_energy])<<","<<t[i_distance][i_angle][i_energy]/double(n[i_distance][i_angle][i_energy])<<std::endl;
		}
	}
}
file_e.close();
std::cout <<"--------------------------------------------------------------------"<<std::endl;


//fill muon csv file
std::cout <<"MUONS: ----------------------------------------------------"<<std::endl;
ofstream file_mu("muons_uniform_bins_20x_20y_10z.csv");
for (int i_distance = 0; i_distance < num_dist/5.; i_distance++){
	for (int i_angle = 0; i_angle < num_angle/5; i_angle++){
		for (int i_energy = 0; i_energy < num_E/2; i_energy++){
			std::cout <<"dist: "<<dist_rebin[i_distance]/5.<<", angle: "<<angle_rebin[i_angle]/5.<<", energy: "<<energy_rebin[i_energy]/2.<<", charge: "<<chargeMu[i_distance][i_angle][i_energy]/double(n[i_distance][i_angle][i_energy])<<", time: "<<tMu[i_distance][i_angle][i_energy]/double(n[i_distance][i_angle][i_energy])<<std::endl;
			file_mu<<dist_rebin[i_distance]/5.<<","<<angle_rebin[i_angle]/5.<<","<<energy_rebin[i_energy]/2.<<","<<chargeMu[i_distance][i_angle][i_energy]/double(n[i_distance][i_angle][i_energy])<<","<<tMu[i_distance][i_angle][i_energy]/double(n[i_distance][i_angle][i_energy])<<std::endl;
		}
	}
}
file_mu.close();
std::cout <<"--------------------------------------------------------------------"<<std::endl;

electron_file->Close();
muon_file->Close();


}
