void read_to_csv(){


std::string electron_filename = "pdf_electron_x50_y50_z10_Parametric_complete.root";
std::string muon_filename = "pdf_muon_x50_y50_z10_Parametric_complete.root";

TH3F *electrons_pe, *muons_pe, *electrons_time, *muons_time;


TFile *electron_file = new TFile(electron_filename.c_str(),"READ");
TFile *muon_file = new TFile(muon_filename.c_str(),"READ");

electron_file->cd();
electrons_pe = (TH3F*) electron_file->Get("electrons_pe");
electrons_time = (TH3F*) electron_file->Get("electrons_time");
muon_file->cd();
muons_pe = (TH3F*) muon_file->Get("muons_pe");
muons_time = (TH3F*) muon_file->Get("muons_time");

const int num_E = 10;
const int num_dist = 50;
const int num_angle = 50;

//fill electron csv file
std::cout <<"ELECTRONS: ----------------------------------------------------"<<std::endl;
ofstream file_e("electrons_uniform_x50_y50_z10_Parametric.csv");
for (int i_distance = 0; i_distance < num_dist; i_distance++){
	for (int i_angle = 0; i_angle < num_angle; i_angle++){
		for (int i_energy = 0; i_energy < num_E; i_energy++){
			std::cout <<"dist: "<<electrons_pe->GetXaxis()->GetBinCenter(i_distance+1)<<", angle: "<<electrons_pe->GetYaxis()->GetBinCenter(i_angle+1)<<", energy: "<<electrons_pe->GetZaxis()->GetBinCenter(i_energy+1)<<", charge: "<<electrons_pe->GetBinContent(i_distance+1,i_angle+1,i_energy+1)<<", time: "<<electrons_time->GetBinContent(i_distance+1,i_angle+1,i_energy+1)<<std::endl;
			file_e<<electrons_pe->GetXaxis()->GetBinCenter(i_distance+1)<<","<<electrons_pe->GetYaxis()->GetBinCenter(i_angle+1)<<","<<electrons_pe->GetZaxis()->GetBinCenter(i_energy+1)<<","<<electrons_pe->GetBinContent(i_distance+1,i_angle+1,i_energy+1)<<","<<electrons_time->GetBinContent(i_distance+1,i_angle+1,i_energy+1)<<std::endl;
		}
	}
}
file_e.close();
std::cout <<"--------------------------------------------------------------------"<<std::endl;


//fill muon csv file
std::cout <<"MUONS: ----------------------------------------------------"<<std::endl;
ofstream file_mu("muons_uniform_x50_y50_z10_Parametric.csv");
for (int i_distance = 0; i_distance < num_dist; i_distance++){
	for (int i_angle = 0; i_angle < num_angle; i_angle++){
		for (int i_energy = 0; i_energy < num_E; i_energy++){
			std::cout <<"dist: "<<muons_pe->GetXaxis()->GetBinCenter(i_distance+1)<<", angle: "<<muons_pe->GetYaxis()->GetBinCenter(i_angle+1)<<", energy: "<<muons_pe->GetZaxis()->GetBinCenter(i_energy+1)<<", charge: "<<muons_pe->GetBinContent(i_distance+1,i_angle+1,i_energy+1)<<", time: "<<muons_time->GetBinContent(i_distance+1,i_angle+1,i_energy+1)<<std::endl;
			file_mu<<muons_pe->GetXaxis()->GetBinCenter(i_distance+1)<<","<<muons_pe->GetYaxis()->GetBinCenter(i_angle+1)<<","<<muons_pe->GetZaxis()->GetBinCenter(i_energy+1)<<","<<muons_pe->GetBinContent(i_distance+1,i_angle+1,i_energy+1)<<","<<muons_time->GetBinContent(i_distance+1,i_angle+1,i_energy+1)<<std::endl;
			
		}
	}
}
file_mu.close();
std::cout <<"--------------------------------------------------------------------"<<std::endl;

electron_file->Close();
muon_file->Close();


}
