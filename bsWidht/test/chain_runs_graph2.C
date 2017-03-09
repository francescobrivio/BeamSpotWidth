/*
    Macro to: 
	- read the tracFile.root,
	- select the tracks,
	- create the possible pairs of tracks per vertex (combinatorial)
	- create tgraphs of <dxy1*dxy2> vs cos(phi1-phi2) and vs cos(phi1+phi2)	
	- fit the graphs
*/

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TTree.h"
#include "TBranch.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TProfile.h"
#include "TF1.h"

#include <string>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#include <fstream>

// Global variables
double pt_cut = 1.;
double eta_cut = 1.;
int trk_hits_cut = 8;
int pix_hits_cut = 1;

int dxy_bins = 2000;
double dxy_min = -0.03;
double dxy_max = 0.03;

int cos_bins = 40;

void chain_runs_graph2 ()
{

	std::cout<<"------------ BEGIN OF MACRO ------------"<<std::endl;

	// INPUT file and tree
	/*int total_trks = 0;
	int tem_min_run = 375376;
	int tem_max_run = 0;
	int tem_irun;
	int tem_irun2;*/

	// TEXT output file
	ofstream myfile;
	myfile.open("Fill5199_output.txt");
	myfile<<"***************** Analysing Fill5199 data *****************\n";
	myfile<<"------------ BEGIN OF MACRO ------------"<<std::endl;

	// TChain
	TChain chain("demo/trackTree");   // name of the tree is the argument
	for (int n = 1; n< 107 ;n++)
	//for (int n = 102; n< 103 ;n++)
	{
	TString file_path = "";
	file_path.Form("root://xrootd-cms.infn.it///store/user/fbrivio/BeamSpot/ZeroBias/crab_Fill_5199_bsWidth/170221_144118/0000/tracksFile_%d.root",n);
	//file_path.Form("/gwteras/cms/store/user/fbrivio/BeamSpot/ZeroBias/crab_2016B_Fills/160701_084313/0000/tracksFile_%d.root",n);
	std::cout <<"Adding: "<< file_path <<"\n";
	myfile<<"Adding: "<< file_path <<"\n";
	//TFile* tem_file = TFile::Open(file_path,"READ");
	//TTree *tem_tree = (TTree*) tem_file->Get("demo/trackTree_");

	//total_trks = total_trks + tem_tree->GetEntriesFast();
	//std::cout<< "\t Total Tracks up to now: "<< total_trks <<"\n";
	
	chain.Add(file_path);
	}

	std::cout<<"TChain built \n";
  	std::cout << "\t Input file: " << chain.GetName() << std::endl;
	myfile<<"TChain built \n";
	myfile<< "\t Input file: " << chain.GetName() << std::endl;

	// OUTPUT file
	TFile* outfile = TFile::Open("Fill5199_graphs.root","RECREATE");
  	std::cout << "\t Output file: " << outfile->GetName() << std::endl;
	myfile<< "\t Output file: " << outfile->GetName() << std::endl;

	// TH2D histograms
	double pt_bins[63] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0 , 5.5, 6.0 , 6.5, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
	TH2D *errIP_pt  = new TH2D("errIP_pt" ,"",62,pt_bins,100,0.,2.);
	TH2D *errIP_eta = new TH2D("errIP_eta","",40.,-2.5,2.5,100,0.,2.);
	TH2D *errIP_phi = new TH2D("errIP_phi","",50,-3.5,3.5,100,0.,2.);

	// TH1D per IP
	h_bs 	 = new 	TH1D("h_bs", "h_d0_bs", 100, -0.05, 0.05);

	// TProfiles
	prof_plus_bs = new TProfile("prof_plus_bs","Profile of <d_{xy}^{1} d_{xy}^{2}> versus cos(#phi_{1}+#phi_{2})",cos_bins,-1.,1.,dxy_min,dxy_max);
	prof_plus_bs->GetXaxis()->SetTitle("cos(#phi_{1}+#phi_{2})");
	prof_plus_bs->GetYaxis()->SetTitle("<d_{xy}^{1} d_{xy}^{2}>");
	prof_plus_bs->SetMarkerStyle(22);
	prof_plus_bs->SetMarkerColor(kBlue);
	prof_minus_bs = new TProfile("prof_minus_bs","Profile of <d_{xy}^{1} d_{xy}^{2}> versus cos(#phi_{1}-#phi_{2})",cos_bins,-1.,1.,dxy_min,dxy_max);
	prof_minus_bs->GetXaxis()->SetTitle("cos(#phi_{1}#pm#phi_{2})");
	prof_minus_bs->GetYaxis()->SetTitle("<d_{xy}^{1} d_{xy}^{2}>");
	prof_minus_bs->SetMarkerStyle(20);
	prof_minus_bs->SetMarkerColor(kRed);

	prof_plus_tt_d0_bs = new TProfile("prof_plus_tt_d0_bs","Profile of <d_{xy}^{1} d_{xy}^{2}> versus cos(#phi_{1}+#phi_{2})",cos_bins,-1.,1.,dxy_min,dxy_max);
	prof_plus_tt_d0_bs->GetXaxis()->SetTitle("cos(#phi_{1}+#phi_{2})");
	prof_plus_tt_d0_bs->GetYaxis()->SetTitle("<d_{xy}^{1} d_{xy}^{2}>");
	prof_minus_tt_d0_bs = new TProfile("prof_minus_tt_d0_bs","Profile of <d_{xy}^{1} d_{xy}^{2}> versus cos(#phi_{1}+#phi_{2})",cos_bins,-1.,1.,dxy_min,dxy_max);
	prof_minus_tt_d0_bs->GetXaxis()->SetTitle("cos(#phi_{1}#pm#phi_{2})");
	prof_minus_tt_d0_bs->GetYaxis()->SetTitle("<d_{xy}^{1} d_{xy}^{2}>");

	// TGraph
	run_graph_x = new TGraphErrors();
	run_graph_x->SetName("Sigma X");
	run_graph_x->SetTitle("#sigma_{x}");
	run_graph_x->GetXaxis()->SetTitle("Run");
	run_graph_x->GetYaxis()->SetTitle("#mum");

	run_graph_y = new TGraphErrors();
	run_graph_y->SetName("Sigma Y");
	run_graph_y->SetTitle("#sigma_{y}");
	run_graph_y->GetXaxis()->SetTitle("Run");
	run_graph_y->GetYaxis()->SetTitle("#mum");

	// Check number of entries in ttree and the number of vertexes
	int nentries = chain.GetEntries();
	//int nentries = 1361496470;
  	std::cout << "\t Number of entries  = " << nentries << std::endl;
	myfile<< "\t Number of entries  = " << nentries << std::endl;
	//int nvtxs = chain.GetMaximum("VtxID");
	//std::cout << "\t Number of vertexes = " << nvtxs << std::endl;
	int min_run = chain.GetMinimum("Run");
	int max_run = chain.GetMaximum("Run");
	//int min_run = 275376;
	//int max_run = 99;
	std::cout << "\t Run numbers: "<<min_run<<" - "<<max_run<<std::endl;
	myfile<< "\t Run numbers: "<<min_run<<" - "<<max_run<<std::endl;

	int point = 0;
	int tmp_run = min_run;
	int Runi = 0;

	// Declaration of variables
	double ipt, ieta, iphi, id0_bs, itt_d0_bs, itt_d0_err_bs;
	int iPix_HITs, iTrack_HITs, iVtxID, iRun;

	// Access branches of ttree
	chain.SetBranchAddress("Run",&iRun);			// Branches for loops
	chain.SetBranchAddress("VtxID",&iVtxID);
	chain.SetBranchAddress("pt",&ipt);			// Branches for track selection
	chain.SetBranchAddress("eta",&ieta);
	chain.SetBranchAddress("Pix_HITs",&iPix_HITs);
	chain.SetBranchAddress("Track_HITs",&iTrack_HITs);
	chain.SetBranchAddress("phi",&iphi);			// Branches for measure
	chain.SetBranchAddress("d0_bs",&id0_bs);
	chain.SetBranchAddress("tt_d0_bs",&itt_d0_bs);
	chain.SetBranchAddress("tt_d0_err_bs",&itt_d0_err_bs);

	// Loop on tracks
	std::cout<<" - Begin of Track Loop - "<<std::endl;
	myfile<<" - Begin of Track Loop - "<<std::endl;
	int after_cuts = 0;

	int first_run = 0;
	int last_run = 0;
	for (int i = 0; i < nentries; ++i) 
	{
		chain.GetEntry(i);

		//if (i == 50000) break;

		if(i%100000 == 0)
		{
			std::cout<<"Track: "<<i<<" - run - tmp_run: "<<Runi<<" - "<<tmp_run<<"\n";
			myfile<<"Track: "<<i<<" - run - tmp_run: "<<Runi<<" - "<<tmp_run<<"\n";
		}

		errIP_pt->Fill(ipt,itt_d0_err_bs);
		errIP_eta->Fill(ieta,itt_d0_err_bs);
		errIP_phi->Fill(iphi,itt_d0_err_bs);

		Runi  = iRun;

		// Track selection (pt>1 - |eta|<1 - 8Hits - 1PixHit)
		if( (ipt > pt_cut) && (std::abs(ieta) < eta_cut) && (iTrack_HITs >= trk_hits_cut) && (iPix_HITs >= pix_hits_cut) ) 
		{

		after_cuts += 1;

		// Temporary variables
		int 	tmp_vtx = iVtxID;
		double 	tmp_phi = iphi;
		double	tmp_d0_bs = id0_bs;
		double	tmp_tt_d0_bs = itt_d0_bs;


		h_bs->Fill(id0_bs);

		// Second loop on tracks
		for (int j = i+1; j< nentries; ++j)
		{
			chain.GetEntry(j);

			// Check if the track is from same vertex
			if (iVtxID > tmp_vtx) break;

			// Track + vertex selection
			if ( (iVtxID == tmp_vtx) && (ipt > pt_cut) && (std::abs(ieta)< eta_cut) && (iTrack_HITs >= trk_hits_cut) && (iPix_HITs >= pix_hits_cut) )
			{

			// Fill TProfiles
			prof_plus_bs 	->Fill(std::cos(iphi + tmp_phi), id0_bs*tmp_d0_bs 	);
			prof_minus_bs	->Fill(std::cos(iphi - tmp_phi), id0_bs*tmp_d0_bs 	);
			prof_plus_tt_d0_bs ->Fill(std::cos(iphi + tmp_phi), itt_d0_bs*tmp_tt_d0_bs);
			prof_minus_tt_d0_bs->Fill(std::cos(iphi - tmp_phi), itt_d0_bs*tmp_tt_d0_bs);

			if (iRun == min_run) first_run += 1;
			//else if (iRun == max_run) last_run += 1;

			} //conditions for second loop
			
		} //end of second loop on tracks

		} //coditions for first loop

		if (tmp_run == Runi)
		{
			tmp_vtx      = 0;
			tmp_phi      = 0.;
			tmp_d0_bs    = 0.;
			tmp_tt_d0_bs = 0.;
			continue;
		}	
		else
		{
			std::cout<<"\n - NEW POINT IN GRAPH, run - tmp_run: "<<Runi<<" - "<<tmp_run<<std::endl;
			myfile<<"\n - NEW POINT IN GRAPH, run - tmp_run: "<<Runi<<" - "<<tmp_run<<std::endl;
		        try{
			prof_minus_bs->Fit("pol1");
   			prof_plus_bs ->Fit("pol1");

			TF1 *f_m_bs = prof_minus_bs->GetFunction("pol1");
			TF1 *f_p_bs = prof_plus_bs->GetFunction("pol1");
			f_p_bs->SetLineColor(kBlue);
			f_p_bs->SetLineColor(kRed);
			double m_bs = f_m_bs->GetParameter(1);
			double m_bs_err = f_m_bs->GetParError(1);
			double p_bs = f_p_bs->GetParameter(1);
			double p_bs_err = f_p_bs->GetParError(1);
			double sigma_x_err = 0.5*sqrt( ((m_bs_err*m_bs_err) + (p_bs_err * p_bs_err))/(m_bs - p_bs) );
			double sigma_y_err = 0.5*sqrt( ((m_bs_err*m_bs_err) + (p_bs_err * p_bs_err))/(m_bs + p_bs) );

			std::cout<<"Sigma_x = "<<sqrt(m_bs - p_bs) * 10000<<" +/- "<<sigma_x_err * 10000<<" um"<<std::endl;
			std::cout<<"Sigma_y = "<<sqrt(m_bs + p_bs) * 10000<<" +/- "<<sigma_y_err * 10000<<" um"<<std::endl;
			myfile<<"Sigma_x = "<<sqrt(m_bs - p_bs) * 10000<<" +/- "<<sigma_x_err * 10000<<" um"<<std::endl;
			myfile<<"Sigma_y = "<<sqrt(m_bs + p_bs) * 10000<<" +/- "<<sigma_y_err * 10000<<" um"<<std::endl;

			std::cout<<"\t Run:"<<tmp_run<<std::endl;
			std::cout<<"\t x:"<<sqrt(m_bs - p_bs)*10000<<std::endl;
			std::cout<<"\t y:"<<sqrt(m_bs + p_bs)*10000<<std::endl;
			myfile<<"\t Run:"<<tmp_run<<std::endl;
			myfile<<"\t x:"<<sqrt(m_bs - p_bs)*10000<<std::endl;
			myfile<<"\t y:"<<sqrt(m_bs + p_bs)*10000<<std::endl;
			run_graph_x->SetPoint(point, tmp_run, sqrt(m_bs - p_bs)*10000);
			run_graph_x->SetPointError(point,0., sigma_x_err*10000);
			run_graph_y->SetPoint(point, tmp_run, sqrt(m_bs + p_bs)*10000);
			run_graph_y->SetPointError(point,0., sigma_y_err*10000);

			prof_minus_bs->Reset();
			prof_plus_bs ->Reset();

			tmp_run = Runi;
			point++;

			// Clear temporary variables
			tmp_vtx      = 0;
			tmp_phi      = 0.;
			tmp_d0_bs    = 0.;
			tmp_tt_d0_bs = 0.;
			}
			catch (std::exception& e)
			{
				std::cerr << "Exception catched : " << e.what() << std::endl;
				continue;
			}
		} //end of else to print beamspot to tgraph

	} //end of first loop on tracks

	std::cout<<" - End of Track Loop - \n";
	std::cout<<"\t Number of events after cuts = "<<after_cuts<<"\n";
	std::cout<<"\t \t First Run events : "<<first_run<<"\n";
	std::cout<<"\t \t Last Run events  : "<<last_run<<"\n";
	myfile<<" - End of Track Loop - \n";
	myfile<<"\t Number of events after cuts = "<<after_cuts<<"\n";
	myfile<<"\t \t First Run events : "<<first_run<<"\n";
	myfile<<"\t \t Last Run events  : "<<last_run<<"\n";

	// Fit the TProfiles
	std::cout<<"Fitting...\n";
	myfile<<"Fitting...\n";
   	prof_minus_bs->Fit("pol1");
   	prof_plus_bs ->Fit("pol1");
	prof_minus_tt_d0_bs->Fit("pol1");
	prof_plus_tt_d0_bs ->Fit("pol1");

	// Retrieve fit parameters
	std::cout<<"Retrieving fit parameters...\n";
	myfile<<"Retrieving fit parameters...\n";
	TF1 *f_m_bs = prof_minus_bs->GetFunction("pol1");
	TF1 *f_p_bs = prof_plus_bs->GetFunction("pol1");
	f_p_bs->SetLineColor(kBlue);
	f_p_bs->SetLineColor(kRed);
	double m_bs = f_m_bs->GetParameter(1);
	double m_bs_err = f_m_bs->GetParError(1);
	double m_chi2 = (f_m_bs->GetChisquare())/(f_m_bs->GetNDF());
	double p_bs = f_p_bs->GetParameter(1);
	double p_bs_err = f_p_bs->GetParError(1);
	double p_chi2 = (f_p_bs->GetChisquare())/(f_p_bs->GetNDF());
	double sigma_x_err = 0.5*sqrt( ((m_bs_err*m_bs_err) + (p_bs_err * p_bs_err))/(m_bs - p_bs) );
	double sigma_y_err = 0.5*sqrt( ((m_bs_err*m_bs_err) + (p_bs_err * p_bs_err))/(m_bs + p_bs) );

	std::cout<<"-------- LAST ONE ---------\n";
	std::cout<<"\n - NEW POINT IN GRAPH, run - tmp_run: "<<Runi<<" - "<<tmp_run<<std::endl;
	std::cout<<"Sigma_x = "<<sqrt(m_bs - p_bs) * 10000<<" +/- "<<sigma_x_err * 10000<<" um"<<std::endl;
	std::cout<<"Sigma_y = "<<sqrt(m_bs + p_bs) * 10000<<" +/- "<<sigma_y_err * 10000<<" um"<<std::endl;
	myfile<<"-------- LAST ONE ---------\n";
	myfile<<"\n - NEW POINT IN GRAPH, run - tmp_run: "<<Runi<<" - "<<tmp_run<<std::endl;
	myfile<<"Sigma_x = "<<sqrt(m_bs - p_bs) * 10000<<" +/- "<<sigma_x_err * 10000<<" um"<<std::endl;
	myfile<<"Sigma_y = "<<sqrt(m_bs + p_bs) * 10000<<" +/- "<<sigma_y_err * 10000<<" um"<<std::endl;
	run_graph_x->SetPoint(point, tmp_run, sqrt(m_bs - p_bs)*10000);
	run_graph_x->SetPointError(point,0., sigma_x_err*10000);
	run_graph_y->SetPoint(point, tmp_run, sqrt(m_bs + p_bs)*10000);
	run_graph_y->SetPointError(point,0., sigma_y_err*10000);

	TF1 *f_m_tt = prof_minus_tt_d0_bs->GetFunction("pol1");
	TF1 *f_p_tt = prof_plus_tt_d0_bs ->GetFunction("pol1");
	double m_tt = f_m_tt->GetParameter(1);
	double p_tt = f_p_tt->GetParameter(1);

	std::cout<<"\n*************** - FIT RESULTS - ***************"<<std::endl;
	std::cout<<"	- bs - "<<std::endl;
	std::cout<<"m_bs = "<<m_bs<<" +/- "<<m_bs_err<<std::endl;
	std::cout<<"m_chi2 = "<<f_m_bs->GetChisquare()<<"/"<<f_m_bs->GetNDF()<<" = "<<m_chi2<<std::endl;
	std::cout<<"p_bs = "<<p_bs<<" +/- "<<p_bs_err<<std::endl;
	std::cout<<"p_chi2 = "<<f_p_bs->GetChisquare()<<"/"<<f_p_bs->GetNDF()<<" = "<<p_chi2<<std::endl;
	std::cout<<"Sigma_x^2 = "<<m_bs - p_bs<<std::endl;
	std::cout<<"Sigma_y^2 = "<<m_bs + p_bs<<std::endl;
	std::cout<<"Sigma_x = "<<sqrt(m_bs - p_bs) * 10000<<" +/- "<<sigma_x_err * 10000<<" um"<<std::endl;
	std::cout<<"Sigma_y = "<<sqrt(m_bs + p_bs) * 10000<<" +/- "<<sigma_y_err * 10000<<" um"<<std::endl;
	std::cout<<"\n";
	std::cout<<"	- tt_bs - "<<std::endl;
	std::cout<<"Sigma_x^2 = "<<m_tt - p_tt<<std::endl;
	std::cout<<"Sigma_y^2 = "<<m_tt + p_tt<<std::endl;
	std::cout<<"Sigma_x = "<<sqrt(m_tt - p_tt) * 10000<<" um"<<std::endl;
	std::cout<<"Sigma_y = "<<sqrt(m_tt + p_tt) * 10000<<" um"<<std::endl;
	std::cout<<"*************************************************"<<std::endl;

	myfile <<"\n*************** - FIT RESULTS - ***************\n";
	myfile<<"	- bs - "<<std::endl;
	myfile<<"m_bs = "<<m_bs<<" +/- "<<m_bs_err<<std::endl;
	myfile<<"m_chi2 = "<<f_m_bs->GetChisquare()<<"/"<<f_m_bs->GetNDF()<<" = "<<m_chi2<<std::endl;
	myfile<<"p_bs = "<<p_bs<<" +/- "<<p_bs_err<<std::endl;
	myfile<<"p_chi2 = "<<f_p_bs->GetChisquare()<<"/"<<f_p_bs->GetNDF()<<" = "<<p_chi2<<std::endl;
	myfile<<"Sigma_x^2 = "<<m_bs - p_bs<<std::endl;
	myfile<<"Sigma_y^2 = "<<m_bs + p_bs<<std::endl;
	myfile<<"Sigma_x = "<<sqrt(m_bs - p_bs) * 10000<<" +/- "<<sigma_x_err * 10000<<" um"<<std::endl;
	myfile<<"Sigma_y = "<<sqrt(m_bs + p_bs) * 10000<<" +/- "<<sigma_y_err * 10000<<" um"<<std::endl;
	myfile<<"*************************************************\n";

	// Projection of TH2D of errIP
	TH1D * proj_errIP_pt = errIP_pt->ProjectionX();
	TH1D * proj_errIP_phi = errIP_phi->ProjectionX();
	TH1D * proj_errIP_eta = errIP_eta->ProjectionX();

	// Close file and save the results
	outfile			-> cd();

	proj_errIP_pt		-> Write();
	proj_errIP_eta		-> Write();
	proj_errIP_phi		-> Write();

	h_bs			-> Write();
	
	prof_plus_bs		-> Write();
	prof_minus_bs		-> Write();
	prof_plus_tt_d0_bs 	-> Write();
	prof_minus_tt_d0_bs	-> Write();

	run_graph_x		-> Write();
	run_graph_y		-> Write();
	
  	outfile			-> Close();  
  
	std::cout<<"------------ END OF MACRO ------------"<<std::endl;
	myfile<<"------------ END OF MACRO ------------"<<std::endl;
	myfile.close();
	return;
}




















