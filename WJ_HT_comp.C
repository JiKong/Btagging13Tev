//root
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TF1.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TBox.h>
#include <TMath.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TList.h>
#include <TChain.h>
#include <TCanvas.h>
#include <THStack.h>  // for stack hist
#include <TGraphErrors.h>  // for TGraphErrors
#include <TFrame.h>  // for TFrameTImage 
#include <TImage.h>  // for TImage 
#include <TString.h> // for TString
#include <TClonesArray.h> //// for TClonesArray
#include <sstream> // for sstream
#include <TBranch.h>
#include <TRandom.h>
#include <TProfile.h>
#include <TLorentzVector.h>
#include <TSystemDirectory.h>
#include <TLegend.h>

//include interface in 13 Tev
#include "interface13/untuplizer.h"
#include "ROG_class.h"

//////////=======rog
string int_to_string(int a)
{
	stringstream ss;
	ss << a;
	string s;
	ss >> s;
	ss.clear();
	return s;
}

void save_hist(TH1D* h_, TString dirname_, TString filename_) // type: eff mr
{
	   
	TCanvas *c1 = new TCanvas(filename_,filename_,200,10,700,500);

	//c1->SetFillColor(42);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor(21);
	c1->GetFrame()->SetBorderSize(12);

	h_->Draw();  // overlay
	
	TImage *img1 = TImage::Create();
	img1->FromPad(c1);
	img1->WriteImage("./"+dirname_+"/"+filename_+".png");
	delete c1;
	delete img1;
/*	
	//save hist in root file 
	TFile* this_rf=new TFile("./"+dirname_+"_rf/"+filename_+".root", "recreate");
	h_->Write();
	this_rf->Close();
*/
}

//main =======================================================================================================================================================
void WJ_HT_comp()
{
		TH1D* h_inclusive_HT = new TH1D("inclusive_HT","inclusive_HT", 3010,-10,3000);
		TH1D* h_HT_HT = new TH1D("HT_HT","HT_HT", 3010,-10,3000);


		// WJ(W+jets) 
		vector<string> inclusive_sample;
		vector<string> HT_sample;
		
		// inclusive
	  	for (int i=1;i<=569;i++)
	  	{
	  		if(i==6 || i==93 || i==324 || i==514 || i==533 || i==567){continue;}
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_221450/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		
	  		// to count Nev of this sample
	  		inclusive_sample.push_back(path);		
	  	}	
	  	
	  	
		// HT 100 - 200
	  	for (int i=1;i<=235;i++)
	  	{
	  		//if(i==6 || i==93 || i==324 || i==514 || i==533 || i==567){continue;}
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235712/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		
	  		// to count Nev of this sample
	  		HT_sample.push_back(path);		
	  	}
	  	
		// HT 200 - 400
	  	for (int i=1;i<=114;i++)
	  	{
	  		if(i==7 || i==8 || i==9 || (i>=65 && i<=99) ){continue;}
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235758/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		
	  		// to count Nev of this sample
	  		HT_sample.push_back(path);		
	  	}
	  	
	  	// HT 1200 - 2500
	  	for (int i=1;i<=7;i++)
	  	{
	  		//if(i==6 || i==93 || i==324 || i==514 || i==533 || i==567){continue;}
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151026_000152/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		
	  		// to count Nev of this sample
	  		HT_sample.push_back(path);		
	  	}


		//=========================================================draw hist=============================================================================	  			
		
		TreeReader inner_data_for_inclusive(inclusive_sample);	
		Long64_t Nev_inclusive = inner_data_for_inclusive.GetEntriesFast();
	  	for (Long64_t Nth = 0; Nth < Nev_inclusive; Nth++)
	    	{	
	    		if(Nth==1){cout<<"****"<<endl;}
	    			
			inner_data_for_inclusive.GetEntry(Nth);
			
			Float_t HT = inner_data_for_inclusive.GetFloat("HT");
			h_inclusive_HT->Fill(HT);
		
		}
		
	
		TreeReader inner_data_for_HT(HT_sample);
		Long64_t Nev_HT = inner_data_for_HT.GetEntriesFast();
	  	for (Long64_t Nth = 0; Nth < Nev_HT; Nth++)
	    	{		
		
			inner_data_for_HT.GetEntry(Nth);
			Float_t HT = inner_data_for_HT.GetFloat("HT");
			h_HT_HT->Fill(HT);
		
		}



		save_hist(h_inclusive_HT, "pic_HT_comp", "inclusive_HT");
		save_hist(h_HT_HT, "pic_HT_comp", "HT_HT");






}
