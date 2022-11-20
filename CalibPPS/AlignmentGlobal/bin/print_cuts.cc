#include <string>
#include <vector>

#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"

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
	
	
	// prepare input
	vector<string> input_files;
	for (int i = 1; i < argc; ++i){
        cout << argv[i] << "\n";
	    input_files.push_back(argv[i]);
	}
	
	// prepare output
	TH2F * h2_xn_vs_xf[2], * h2_yn_vs_yf[2];
	int nbins = 100; float xmin=0, xmax=30, ymin=-20,ymax=20;
    for (int i=0;i<2;i++){
		int sector = 45+i*11;
		h2_xn_vs_xf[i] = new TH2F(Form("h2_%d_xf_vs_xn",sector),Form("Sector %d;x_{NEAR};x_{FAR}",sector),nbins,xmin,xmax,nbins,xmin,xmax);
		h2_yn_vs_yf[i] = new TH2F(Form("h2_%d_yf_vs_yn",sector),Form("Sector %d;y_{NEAR};y_{FAR}",sector),nbins,ymin,ymax,nbins,ymin,ymax);
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
	float Sum_Xn[2]={0,0}, Sum_Xf[2]={0,0}, Sum_XfXn[2]={0,0}, Sum_XnXn[2]={0,0}, N[2]={0,0};
	float Sum_Yn[2]={0,0}, Sum_Yf[2]={0,0}, Sum_YfYn[2]={0,0}, Sum_YnYn[2]={0,0};
	for (event.toBegin(); ! event.atEnd(); ++event)
	{
      
		// load track data
		fwlite::Handle< vector<CTPPSLocalTrackLite> > tracks;
		tracks.getByLabel(event, track_module_name.c_str());

        float xf[2], xn[2], yf[2], yn[2]; int Ntraks_near[2]={0,0}, Ntraks_far[2]={0,0};
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
				Ntraks_near[arm]++;
			}
			else if(rpId.station()==2){ // far station
				xf[arm] = _x;
				yf[arm] = _y;
				Ntraks_far[arm]++;
			}
		}
		
		for (int i=0;i<2;i++){
		    if(Ntraks_near[i]==1 &&  Ntraks_far[i]==1){
			    N[i]++;
				Sum_Xn[i]+=xn[i];
				Sum_Xf[i]+=xf[i];
				Sum_XfXn[i]+=xf[i]*xn[i];
				Sum_XnXn[i]+=xn[i]*xn[i];
				Sum_Yn[i]+=yn[i];
				Sum_Yf[i]+=yf[i];
				Sum_YfYn[i]+=yf[i]*yn[i];
				Sum_YnYn[i]+=yn[i]*yn[i];
				if(_DEBUG) h2_xn_vs_xf[i]->Fill(xn[i],xf[i]);
				if(_DEBUG) h2_yn_vs_yf[i]->Fill(yn[i],yf[i]);
			}
		}
		
		// Fill output histograms
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
    float cut_h_a[2] =       {
		                       Sum_Xn[0]*Sum_Xf[0] - N[0]*Sum_XfXn[0],
		                       Sum_Xn[1]*Sum_Xf[1] - N[1]*Sum_XfXn[1]
	                         };
    float cut_h_c[2] =       {
		                       Sum_Xn[0]*Sum_XfXn[0] - Sum_XnXn[0]*Sum_Xf[0],
							   Sum_Xn[1]*Sum_XfXn[1] - Sum_XnXn[1]*Sum_Xf[1]
	                         };
    float cut_v_a[2] =       {
		                       Sum_Yn[0]*Sum_Yf[0] - N[0]*Sum_YfYn[0],
		                       Sum_Yn[1]*Sum_Yf[1] - N[1]*Sum_YfYn[1]
	                         };
    float cut_v_c[2] =       {
		                       Sum_Yn[0]*Sum_YfYn[0] - Sum_YnYn[0]*Sum_Yf[0],
							   Sum_Yn[1]*Sum_YfYn[1] - Sum_YnYn[1]*Sum_Yf[1]
	                         };
	for(int i=0;i<2;i++){
		cut_h_a[i]/=determinant_x[i];
		cut_h_c[i]/=determinant_x[i];
		cut_v_a[i]/=determinant_y[i];
		cut_v_c[i]/=determinant_y[i];
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
			h2_xn_vs_xf[i]->Write();
			h2_yn_vs_yf[i]->Write();	
		}
		cout << "Writes " << _f->GetName() << endl;
		_f->Write();
		_f->Close();
	}
	
	return 1;
}
