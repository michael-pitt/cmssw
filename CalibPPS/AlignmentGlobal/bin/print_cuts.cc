#include <string>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TF1.h"

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/CTPPSDetId/interface/TotemRPDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"

using namespace std;
using namespace std::chrono;

//----------------------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{

    // User settings
    string track_module_name = "ctppsLocalTrackLiteProducer"; 
    //string track_module_name = "ctppsLocalTrackLiteProducerAlCaRecoProducer"; 
	bool _DEBUG = false;
	float EPSILON = 1E0;
	
	
	// prepare input
	vector<string> input_files;
	for (int i = 1; i < argc; ++i){
		if(strcmp(argv[i],"--debug")==0){_DEBUG=true; continue;}
        cout << argv[i] << "\n";
	    input_files.push_back(argv[i]);
	}
	
	// prepare output
	TH2F * h2_a_vs_c_h[2], * h2_a_vs_c_v[2];
	TH2F * h2_xn_vs_xf_before[2], * h2_yn_vs_yf_before[2];
	TH2F * h2_xn_vs_xf_after[2], * h2_yn_vs_yf_after[2];
	TH1F * h1_Cq_h[2], * h1_Cq_v[2];
	TH1F * h1_FitC_h[2], * h1_FitC_v[2];
	TH1F * h1_FitA_h[2], * h1_FitA_v[2];
	int nbins = 100; float xmin=0, xmax=30, ymin=-20,ymax=20;
    for (int i=0;i<2;i++){
		int sector = 45+i*11;
		h2_xn_vs_xf_before[i] = new TH2F(Form("h2_%d_xf_vs_xn_before",sector),Form("Sector %d;x_{NEAR};x_{FAR}",sector),nbins,xmin,xmax,nbins,xmin,xmax);
		h2_yn_vs_yf_before[i] = new TH2F(Form("h2_%d_yf_vs_yn_before",sector),Form("Sector %d;y_{NEAR};y_{FAR}",sector),nbins,ymin,ymax,nbins,ymin,ymax);
		h1_FitC_h[i] = new TH1F(Form("h1_%d_FitC_h",sector),Form("Sector %d;x_{FAR}-x_{NEAR};entries",sector),500,-2.,2.);
		h1_FitC_v[i] = new TH1F(Form("h1_%d_FitC_v",sector),Form("Sector %d;y_{FAR}-y_{NEAR};entries",sector),500,-2.,2.);
		h1_FitA_h[i] = new TH1F(Form("h1_%d_FitA_h",sector),Form("Sector %d;#left(x_{FAR}+c#right)/x_{NEAR};entries",sector),500,-2.,2.);
		h1_FitA_v[i] = new TH1F(Form("h1_%d_FitA_v",sector),Form("Sector %d;#left(y_{FAR}+c#right)/y_{NEAR};entries",sector),500,-2.,2.);
		h1_Cq_h[i] = new TH1F(Form("h1_%d_Cq_h",sector),Form("Sector %d;Cq_{h} = x_{FAR} + a#upointx_{NEAR} + c;events",sector),nbins,-20.,20.);
		h1_Cq_v[i] = new TH1F(Form("h1_%d_Cq_v",sector),Form("Sector %d;Cq_{v} = y_{FAR} + a#upointy_{NEAR} + c;events",sector),nbins,-20.,20.);
		h2_xn_vs_xf_after[i] = new TH2F(Form("h2_%d_xf_vs_xn_after",sector),Form("Sector %d;x_{NEAR};x_{FAR}",sector),nbins,xmin,xmax,nbins,xmin,xmax);
		h2_yn_vs_yf_after[i] = new TH2F(Form("h2_%d_yf_vs_yn_after",sector),Form("Sector %d;y_{NEAR};y_{FAR}",sector),nbins,ymin,ymax,nbins,ymin,ymax);
		h2_a_vs_c_h[i] = new TH2F(Form("h2_%d_a_vs_c_h",sector),Form("Sector %d;a_{H};c_{H}",sector),1000,-20.,20.,1000,-20.,20.);
		h2_a_vs_c_v[i] = new TH2F(Form("h2_%d_a_vs_c_v",sector),Form("Sector %d;a_{V};c_{V}",sector),1000,-20.,20.,1000,-20.,20.);
	}
	
	if(!input_files.size()){
		cout << "ERROR: no input files provided, exit" << endl;
		return 0;
	}
	
	// Get CMSSW version
	string version = getenv("CMSSW_RELEASE_BASE");
	if(version.empty()){
		cout << "ERROR: first setup CMSSW (any version you like)"<<endl;
		return 0;
	}
    
	cout << "Load events from " << input_files.size() << " files"<< endl;
	fwlite::ChainEvent event(input_files);

	printf("* events in input chain: %llu\n", event.size());

    auto start = high_resolution_clock::now();
	
	// Set counters (for 45 and 56):
	float Sum_Xn[2]={0,0}, Sum_Xf[2]={0,0}, Sum_XfXn[2]={0,0}, Sum_XfXf[2]={0,0}, Sum_XnXn[2]={0,0}, N[2]={0,0};
	float Sum_Yn[2]={0,0}, Sum_Yf[2]={0,0}, Sum_YfYn[2]={0,0}, Sum_YfYf[2]={0,0}, Sum_YnYn[2]={0,0};

    // Initiate fit parameters 
    float cut_h_a[2] =       {-1.0, -1.0};
    float cut_h_c[2] =       { 0.0,  0.0};
    float cut_v_a[2] =       {-1.0, -1.0};
    float cut_v_c[2] =       { 0.0,  0.0};

	// Use fast fit to c_h and c_v (a=-1)
	for (event.toBegin(); ! event.atEnd(); ++event)
	{
      
		// load track data
		fwlite::Handle< vector<CTPPSLocalTrackLite> > tracks;
		tracks.getByLabel(event, track_module_name.c_str());

        float xf[2], xn[2], yf[2], yn[2]; int Ntracks_near[2]={0,0}, Ntracks_far[2]={0,0};
		// process track data
		for (const auto &tr : *tracks)
		{
			
			// Works with CMSSW_12_X
			CTPPSDetId rpId(tr.rpId());
            float _x = tr.x(), _y = tr.y();
			
			// Works with CMSSW_10_X
			//CTPPSDetId rpId(tr.getRPId());			
			//float _x = tr.getX(), _y = tr.getY();

			
			if(rpId.rp()!=3) continue;
			
			int arm = rpId.arm();
			if(rpId.station()==0){ // near station
				xn[arm] = _x;
				yn[arm] = _y;
				Ntracks_near[arm]++;
			}
			else if(rpId.station()==2){ // far station
				xf[arm] = _x;
				yf[arm] = _y;
				Ntracks_far[arm]++;
			}
		}
		
		for (int i=0;i<2;i++){
		    if(Ntracks_near[i]==1 &&  Ntracks_far[i]==1){
				h1_FitC_h[i]->Fill(-xf[i]+xn[i]);
				h1_FitC_v[i]->Fill(-yf[i]+yn[i]);
				if(_DEBUG) h2_xn_vs_xf_before[i]->Fill(xn[i],xf[i]);
				if(_DEBUG) h2_yn_vs_yf_before[i]->Fill(yn[i],yf[i]);				
			}
		}
	}
	
	// Fit c_h and c_v
	for (int i=0;i<2;i++){
		h1_FitC_h[i]->Fit("gaus");
		cut_h_c[i] = -h1_FitC_h[i]->GetFunction("gaus")->GetParameter(1);
		h1_FitC_v[i]->Fit("gaus");
		cut_v_c[i] = -h1_FitC_v[i]->GetFunction("gaus")->GetParameter(1);
	}
	
	// Use fast fit to a_h and a_v (using prefitted c_h and c_v)
	for (event.toBegin(); ! event.atEnd(); ++event)
	{
      
		// load track data
		fwlite::Handle< vector<CTPPSLocalTrackLite> > tracks;
		tracks.getByLabel(event, track_module_name.c_str());

        float xf[2], xn[2], yf[2], yn[2]; int Ntracks_near[2]={0,0}, Ntracks_far[2]={0,0};
		// process track data
		for (const auto &tr : *tracks)
		{
			
			// Works with CMSSW_12_X
			CTPPSDetId rpId(tr.rpId());
            float _x = tr.x(), _y = tr.y();
			
			// Works with CMSSW_10_X
			//CTPPSDetId rpId(tr.getRPId());			
			//float _x = tr.getX(), _y = tr.getY();

			
			if(rpId.rp()!=3) continue;
			
			int arm = rpId.arm();
			if(rpId.station()==0){ // near station
				xn[arm] = _x;
				yn[arm] = _y;
				Ntracks_near[arm]++;
			}
			else if(rpId.station()==2){ // far station
				xf[arm] = _x;
				yf[arm] = _y;
				Ntracks_far[arm]++;
			}
		}
		
		for (int i=0;i<2;i++){
		    if(Ntracks_near[i]==1 &&  Ntracks_far[i]==1){
				h1_FitA_h[i]->Fill((xf[i]+cut_h_c[i])/xn[i]);
				h1_FitA_v[i]->Fill((yf[i]+cut_v_c[i])/yn[i]);
			}
		}
	}
	
	for (int i=0;i<2;i++){
		h1_FitA_h[i]->Fit("gaus");
		cut_h_a[i] = -h1_FitA_h[i]->GetFunction("gaus")->GetParameter(1);
		h1_FitA_v[i]->Fit("gaus");
		cut_v_a[i] = -h1_FitA_v[i]->GetFunction("gaus")->GetParameter(1);
	}	
	
	// Loop over all events, rejecting those far from the prefit 
	for (event.toBegin(); ! event.atEnd(); ++event)
	{
      
		// load track data
		fwlite::Handle< vector<CTPPSLocalTrackLite> > tracks;
		tracks.getByLabel(event, track_module_name.c_str());

        float xf[2], xn[2], yf[2], yn[2]; int Ntracks_near[2]={0,0}, Ntracks_far[2]={0,0};
		// process track data
		for (const auto &tr : *tracks)
		{
			
			// Works with CMSSW_12_X
			CTPPSDetId rpId(tr.rpId());
            float _x = tr.x(), _y = tr.y();
			
			// Works with CMSSW_10_X
			//CTPPSDetId rpId(tr.getRPId());			
			//float _x = tr.getX(), _y = tr.getY();

			
			if(rpId.rp()!=3) continue;
			
			int arm = rpId.arm();
			if(rpId.station()==0){ // near station
				xn[arm] = _x;
				yn[arm] = _y;
				Ntracks_near[arm]++;
			}
			else if(rpId.station()==2){ // far station
				xf[arm] = _x;
				yf[arm] = _y;
				Ntracks_far[arm]++;
			}
		}
		
		for (int i=0;i<2;i++){
		    if(Ntracks_near[i]==1 &&  Ntracks_far[i]==1){
				if(_DEBUG) h1_Cq_h[i]->Fill(xf[i]+cut_h_a[i]*xn[i]+cut_h_c[i]);
				if(_DEBUG) h1_Cq_v[i]->Fill(yf[i]+cut_v_a[i]*yn[i]+cut_v_c[i]);
				if(fabs(xf[i]+cut_h_a[i]*xn[i]+cut_h_c[i])>EPSILON) continue;
				if(fabs(yf[i]+cut_v_a[i]*yn[i]+cut_v_c[i])>EPSILON) continue;
			    N[i]++;
				Sum_Xn[i]+=xn[i];
				Sum_Xf[i]+=xf[i];
				Sum_XfXf[i]+=xf[i]*xf[i];
				Sum_XfXn[i]+=xf[i]*xn[i];
				Sum_XnXn[i]+=xn[i]*xn[i];
				Sum_Yn[i]+=yn[i];
				Sum_Yf[i]+=yf[i];
				Sum_YfYn[i]+=yf[i]*yn[i];
				Sum_YnYn[i]+=yn[i]*yn[i];
				if(_DEBUG) h2_xn_vs_xf_after[i]->Fill(xn[i],xf[i]);
				if(_DEBUG) h2_yn_vs_yf_after[i]->Fill(yn[i],yf[i]);
			}
		}
		
	}
	
	// compute constants:
	float determinant_x[2] = { 
	                           N[0]*Sum_XnXn[0]-Sum_Xn[0]*Sum_Xn[0], 
	                           N[1]*Sum_XnXn[1]-Sum_Xn[1]*Sum_Xn[1]
	                         };
    float determinant_y[2] = { 
	                           N[0]*Sum_YnYn[0]-Sum_Yn[0]*Sum_Yn[0], 
	                           N[1]*Sum_YnYn[1]-Sum_Yn[1]*Sum_Yn[1]
	                         };
							 
	// Solution for a and c
    for (int i=0;i<2;i++){
		cut_h_a[i] = Sum_Xn[i]*Sum_Xf[i] - N[i]*Sum_XfXn[i];
		cut_h_c[i] = Sum_Xn[i]*Sum_XfXn[i] - Sum_XnXn[i]*Sum_Xf[i];
		cut_v_a[i] = Sum_Yn[i]*Sum_Yf[i] - N[i]*Sum_YfYn[i];
		cut_v_c[i] = Sum_Yn[i]*Sum_YfYn[i] - Sum_YnYn[i]*Sum_Yf[i];

		cut_h_a[i] /= determinant_x[i];
		cut_h_c[i] /= determinant_x[i];
		cut_v_a[i] /= determinant_y[i];
		cut_v_c[i] /= determinant_y[i];
	}
	
	for (int i=0;i<2 && _DEBUG;i++){ // Scan in a - c phase-space
		for(int ibin_c=0;ibin_c<1000;ibin_c++){
		for(int ibin_a=0;ibin_a<1000;ibin_a++){
			float a = h2_a_vs_c_h[i]->GetXaxis()->GetBinCenter(ibin_a+1);
			float c = h2_a_vs_c_h[i]->GetXaxis()->GetBinCenter(ibin_c+1);
			h2_a_vs_c_h[i]->Fill(a,c,Sum_XfXf[i]+2.*a*Sum_XfXn[i]+2.*c*Sum_Xf[i]+a*a*Sum_XnXn[i]+2.*a*c*Sum_Xn[i]+N[i]*c*c);
			h2_a_vs_c_v[i]->Fill(a,c,Sum_YfYf[i]+2.*a*Sum_YfYn[i]+2.*c*Sum_Yf[i]+a*a*Sum_YnYn[i]+2.*a*c*Sum_Yn[i]+N[i]*c*c);
		}}
	}
				

    auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	
    printf("* events accepted (45, 56): (%i, %i)\n", int(N[0]), int(N[1]));
	printf("-------------------------\n");
	printf("Update config_cff.py with the following constants:\n");
	printf("-------------------------\n");
	printf("ppsAlignmentConfigESSource = cms.ESSource(\"PPSAlignmentConfigurationESSource\",\n");
    printf("    sector_45 = cms.PSet(\n");
	printf("\t...\n");
	printf("\tcut_h_a = cms.double(%f),\n",cut_h_a[0]);
	printf("\tcut_h_c = cms.double(%f),\n",cut_h_c[0]);
	printf("\n");
	printf("\tcut_v_a = cms.double(%f),\n",cut_v_a[0]);
	printf("\tcut_v_c = cms.double(%f)\n",cut_v_c[0]);
	printf("    ),\n");
	printf("\n");
    printf("    sector_56 = cms.PSet(\n");
	printf("\t...\n");
	printf("\tcut_h_a = cms.double(%f),\n",cut_h_a[1]);
	printf("\tcut_h_c = cms.double(%f),\n",cut_h_c[1]);
	printf("\n");
	printf("\tcut_v_a = cms.double(%f),\n",cut_v_a[1]);
	printf("\tcut_v_c = cms.double(%f)\n",cut_v_c[1]);
	printf("    ),\n");
	
	float time = duration.count(); string time_str = "microseconds";
	if (time/1000 > 1) {time/=1000.; time_str = "miliseconds";}
	if (time/1000 > 1) {time/=1000.; time_str = "seconds";}
	if (time/60 > 1) {time/=60.; time_str = "minutes";}
	printf("elapsed time: %2.2f %s\n",time,time_str.c_str());
	
	if(_DEBUG) { // produce output file for debug
		TFile * _f = new TFile("output/print_cuts_debug_plots.root","recreate");
		_f->cd();
		for(int i=0;i<2;i++){
			h2_xn_vs_xf_before[i]->Write();
			h2_yn_vs_yf_before[i]->Write();	
			h2_xn_vs_xf_after[i]->Write();
			h2_yn_vs_yf_after[i]->Write();	
			h1_Cq_h[i]->Write();	
			h1_Cq_v[i]->Write();	
			h2_a_vs_c_h[i]->Write();	
			h2_a_vs_c_v[i]->Write();	
			h1_FitA_h[i]->Write();	
			h1_FitA_v[i]->Write();	
			h1_FitC_h[i]->Write();	
			h1_FitC_v[i]->Write();	
		}
		cout << "Writes " << _f->GetName() << endl;
		_f->Write();
		_f->Close();
	}
	
	return 1;
}
