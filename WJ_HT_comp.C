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

using namespace std;

// global variable================================================

Long_t Nev_inclusive =0;
Long_t Nev_HT100to200 =0;
Long_t Nev_HT200to400 =0;
Long_t Nev_HT400to600 =0;
Long_t Nev_HT600to800 =0;
Long_t Nev_HT800to1200 =0;
Long_t Nev_HT1200to2500 =0;

double available_event_inclusive_scaled =0;
double available_event_HT100to200_scaled =0;
double available_event_HT200to400_scaled =0;
double available_event_HT400to600_scaled =0;
double available_event_HT600to800_scaled =0;
double available_event_HT800to1200_scaled =0;
double available_event_HT1200to2500_scaled =0;

int Nev_index_inclusive =0;
int Nev_index_HT100to200 =0;
int Nev_index_HT200to400 =0;
int Nev_index_HT400to600 =0;
int Nev_index_HT600to800 =0;
int Nev_index_HT800to1200 =0;
int Nev_index_HT1200to2500 =0;

double Xs_inclusive =61526;
double Xs_HT100to200 =1627.45;
double Xs_HT200to400 =435.237;
double Xs_HT400to600 =59.18;
double Xs_HT600to800 =14.58;
double Xs_HT800to1200 =6.65;
double Xs_HT1200to2500 =1.6;

double sf_inclusive =0;
double sf_HT100to200 =0;
double sf_HT200to400 =0;
double sf_HT400to600 =0;
double sf_HT600to800 =0;
double sf_HT800to1200 =0;
double sf_HT1200to2500 =0;

double data_lumi_13_Tev = 831.7;

double Nev_HT100up_scaled =0 ;
double Nev_HT100down_scaled =0 ;

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

double get_scaled_value(Long_t Nth_)
{
	double scaled_value=0;
	if (Nth_>=Nev_index_HT1200to2500) {scaled_value=sf_HT1200to2500;}
	else if(Nth_>=Nev_index_HT800to1200) {scaled_value=sf_HT800to1200;}
	else if(Nth_>=Nev_index_HT600to800) {scaled_value=sf_HT600to800;}
	else if(Nth_>=Nev_index_HT400to600) {scaled_value=sf_HT400to600;}
	else if(Nth_>=Nev_index_HT200to400) {scaled_value=sf_HT200to400;}
	else if(Nth_>=Nev_index_HT100to200) {scaled_value=sf_HT100to200;}
	else if(Nth_>=Nev_index_inclusive) {scaled_value=sf_inclusive; }
	
	return scaled_value;
}


// HEP function
// //general
Float_t delta_R(PO O1_,PO O2_)
{
	Float_t a = (O1_.p4.Eta()-O2_.p4.Eta());
	Float_t b = (O1_.p4.Phi()-O2_.p4.Phi());
	return sqrt(a*a+b*b);
}

Float_t delta_R(TLorentzVector T1_,TLorentzVector T2_)
{
	Float_t a = (T1_.Eta()-T2_.Eta());
	Float_t b = (T1_.Phi()-T2_.Phi());
	return sqrt(a*a+b*b);
}

// // for THINjet
bool is_pass_THINjet_loose_id(THINjet j_)
{
	if (j_.is_loose==true){return true;}
	else {return false;}
}

bool is_pass_THINjet_tight_id(THINjet j_)
{
	if (j_.is_tight==true){return true;}
	else {return false;}
}

bool is_good_jet( THINjet j_)
{
	return (j_.p4.Pt()>30 && abs(j_.p4.Eta())<2.5  && is_pass_THINjet_loose_id(j_));
}


//  //  ele selection

bool is_pass_ele_loose_id(ele e_)
{
	if (e_.is_loose==true){return true;}
	else {return false;}
}

bool is_pass_ele_tight_id(ele e_)
{
	if (e_.is_tight==true){return true;}
	else {return false;}
}

//  //  muon selection  

bool is_pass_mu_idiso(mu m_)
{
	//if (!m_.is_tight==true){return false;}
	if (!m_.is_loose==true){return false;}
	double this_v = (m_.ChHadIso+max(0., m_.NeHadIso + m_.GamIso - 0.5*m_.PUPt))/m_.p4.Pt(); 
	if (! (this_v<0.15) ){return false;}
	return true;

}

//  //  evse rec level
bool is_pass_ev_selection_TTbar(event ev)
{

	// trigger path
	if (!ev.is_pass_mu_trigger){return false;}
	
	// same as theevent selection on ttbar sample, be cause we will see how ttbar selection works.....
	int channel=0; //1 for ele, 2 for mu
	int ev_flag=0;

	// 1: event must have 1 good muon, isolated, PT>50, |Eta|<2.1
	vector<ele> good_eles=ev.good_eles;
	bool have_available_ele=false;
	/*
	//electron channel
	for (int i=0;i<good_eles.size();i++)
	{
		ele this_e=good_eles[i];
		if (this_e.p4.Pt()>50 && abs(this_e.p4.Eta())<2.1)
		{ have_available_ele=true; break;}
	}
	if(!have_available_ele)
	{ev_flag=2; return false;}
	*/
	 
	// muon channel
	vector<mu> good_mus=ev.good_mus;
	bool have_available_muon=false;
	for (int i=0;i<good_mus.size();i++)
	{
		mu this_m=good_mus[i];
		if (this_m.p4.Pt()>50 && abs(this_m.p4.Eta())<2.1)
		{ have_available_muon=true; break;}
	}
	if(!have_available_muon)
	{ev_flag=2; return false;}
	


	// 3: to check if have at least 4 good jet
	// using THINjet
	vector<THINjet> good_jets;
	
	for (int i=0;i<ev.THINjets_deoverlap_with_lepton.size();i++)
	{
	
		THINjet this_j = ev.THINjets_deoverlap_with_lepton[i];		
		if(!is_good_jet(this_j)){continue;}						
		good_jets.push_back(this_j);

	}
	if (good_jets.size()<4){ ev_flag=3; return false;}

	// 4: to check if have 2 good leading jet
	// using THINjet
	bool have_good_leading_jets=false;
	for (int i=0;i<good_jets.size();i++)
	{
		for (int j=0;j<good_jets.size();j++)
		{
			if(i==j){continue;}
			if(good_jets[i].p4.Pt()>=70 && good_jets[j].p4.Pt()>=50
			//&&  is_pass_cisvv2_13Tev(good_jets[i]) &&  is_pass_cisvv2_13Tev(good_jets[j])
			)
			{have_good_leading_jets=true;}
			//if(good_jets[i].p4.Pt()>=200 && good_jets[j].p4.Pt()>=150) {have_good_leading_jets=true;}
		}
	}
	if(!have_good_leading_jets){ev_flag=4;return false;}
	
	
	// 5: to check if missing Et > 20
	//if (ev.missing_et.pt<20){ev_flag=5;return false;}
	if (ev.missing_et.pt<20){ev_flag=5;return false;}

//test
	return true;

}


//main =======================================================================================================================================================
void WJ_HT_comp()
{
		TH1D* h_inclusive_HT = new TH1D("inclusive_HT","inclusive_HT", 100,-10,2500);
		TH1D* h_HT_HT = new TH1D("HT_HT","HT_HT", 100,-10,2500);
		
		double N_available_event_for_inclusive = 0;
		double N_available_event_for_HT = 0;
		
		int Nev_index_counter = 0 ;

		// WJ(W+jets) 
		vector<string> input_files;	  	
	  	
	  			
		// inclusive
	  	bool is_save_Nev_index_inclusive=false;
	  	for (int i=1;i<=569;i++)
	  	{	
	  		//if (i>=10){continue;}
	  		if(i==6 || i==93 || i==324 || i==514 || i==533 || i==567){continue;}
			if (!is_save_Nev_index_inclusive){Nev_index_inclusive=Nev_index_counter; is_save_Nev_index_inclusive=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_221450/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
	  		Nev_inclusive+=this_reader.GetEntriesFast(); 
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
	  	// HT 100 - 200
	  	bool is_save_Nev_index_HT100to200=false;
	  	for (int i=1;i<=235;i++)
	  	{
	  		//if (i>=10){continue;}
			if (!is_save_Nev_index_HT100to200){Nev_index_HT100to200=Nev_index_counter; is_save_Nev_index_HT100to200=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235712/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
	  		Nev_HT100to200+=this_reader.GetEntriesFast(); 
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
	  	// HT 200 - 400
	  	bool is_save_Nev_index_HT200to400=false;
	  	for (int i=1;i<=114;i++)
	  	{
	  		//if (i>=10){continue;}
	  		if(i==7 || i==8 || i==9 || (i>=65 && i<=99) ){continue;}
			if (!is_save_Nev_index_HT200to400){Nev_index_HT200to400=Nev_index_counter; is_save_Nev_index_HT200to400=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235758/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
	  		Nev_HT200to400+=this_reader.GetEntriesFast(); 
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}	
	  	
	  	
	  	// HT 400 - 600
	  	bool is_save_Nev_index_HT400to600=false;
	  	for (int i=1;i<=44;i++)
	  	{
	  		//if (i>=10){continue;}
	  		if(i==38 ){continue;}
			if (!is_save_Nev_index_HT400to600){Nev_index_HT400to600=Nev_index_counter; is_save_Nev_index_HT400to600=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235853/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
	  		Nev_HT400to600+=this_reader.GetEntriesFast(); 
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}


	  	// HT 600 - 800
	  	bool is_save_Nev_index_HT600to800=false;
	  	for (int i=1;i<=97;i++)
	  	{
	  		//if (i>=10){continue;}
	  		
			if (!is_save_Nev_index_HT600to800){Nev_index_HT600to800=Nev_index_counter; is_save_Nev_index_HT600to800=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235938/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
	  		Nev_HT600to800+=this_reader.GetEntriesFast(); 
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}

	  	// HT 800 - 1200
	  	bool is_save_Nev_index_HT800to1200=false;
	  	for (int i=1;i<=39;i++)
	  	{
	  		//if (i>=10){continue;}
			if (!is_save_Nev_index_HT800to1200){Nev_index_HT800to1200=Nev_index_counter; is_save_Nev_index_HT800to1200=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151026_000033/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
	  		Nev_HT800to1200+=this_reader.GetEntriesFast(); 
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}

	  	// HT 1200 - 2500
	  	bool is_save_Nev_index_HT1200to2500=false;
	  	for (int i=1;i<=7;i++)
	  	{
	  		//if (i>=10){continue;}
			if (!is_save_Nev_index_HT1200to2500){Nev_index_HT1200to2500=Nev_index_counter; is_save_Nev_index_HT1200to2500=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151026_000152/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
	  		Nev_HT1200to2500+=this_reader.GetEntriesFast(); 
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
		// ===============================================calculate scale factor ======================================================================

		sf_inclusive =data_lumi_13_Tev/((double)Nev_inclusive/Xs_inclusive);
		sf_HT100to200 =data_lumi_13_Tev/((double)Nev_HT100to200/Xs_HT100to200);
		sf_HT200to400 =data_lumi_13_Tev/((double)Nev_HT200to400 /Xs_HT200to400 );
		sf_HT400to600 =data_lumi_13_Tev/((double)Nev_HT400to600/Xs_HT400to600);
		sf_HT600to800 =data_lumi_13_Tev/((double)Nev_HT600to800/Xs_HT600to800);
		sf_HT800to1200 =data_lumi_13_Tev/((double)Nev_HT800to1200/Xs_HT800to1200);
		sf_HT1200to2500 =data_lumi_13_Tev/((double)Nev_HT1200to2500/Xs_HT1200to2500);

		//=========================================================draw hist=============================================================================	  			
		
		TreeReader inner_data(input_files);	
		Long64_t Nev_all_WJ = inner_data.GetEntriesFast();
	  	//for (Long64_t Nth = 0; Nth < Nev_all_WJ; Nth++)  // draw all root file
	  	//for (Long64_t Nth = 0; Nth < Nev_index_HT100to200; Nth++)  // only draw HT bin sample
	  	for (Long64_t Nth = Nev_index_HT100to200; Nth < Nev_all_WJ; Nth++)  // only draw inclusive sample
	    	{	
			
			// event should pass event selection
			event this_ev;
			inner_data.GetEntry(Nth);
			// trigger var
			string* hlt_trigName = inner_data.GetPtrString("hlt_trigName");
			vector<bool> &hlt_trigResult = *((vector<bool>*) inner_data.GetPtr("hlt_trigResult"));	
			int string_size = inner_data.GetPtrStringSize();
		
	      		// THINjet var
	      		Int_t THINnJet=inner_data.GetInt("THINnJet");
	      		
	      		TClonesArray* THINjetP4 = (TClonesArray*) inner_data.GetPtrTObject("THINjetP4");
	      		Float_t *THINjetCSV=inner_data.GetPtrFloat("THINjetCSV");
			Float_t* THINjetCISVV2 = inner_data.GetPtrFloat("THINjetCISVV2");
			Int_t * THINjetHadronFlavor=inner_data.GetPtrInt("THINjetHadronFlavor"); 			
			
			vector<bool> &THINjetPassIDTight = *((vector<bool>*) inner_data.GetPtr("THINjetPassIDTight"));			
			vector<bool> &THINjetPassIDLoose  = *((vector<bool>*) inner_data.GetPtr("THINjetPassIDLoose"));

			Float_t *THINjetCEmEF=inner_data.GetPtrFloat("THINjetCEmEF");
			Float_t *THINjetCHadEF=inner_data.GetPtrFloat("THINjetCHadEF");	
			Float_t *THINjetNEmEF=inner_data.GetPtrFloat("THINjetNEmEF");	
			Float_t *THINjetNHadEF=inner_data.GetPtrFloat("THINjetNHadEF");
			Float_t *THINjetPhoEF=inner_data.GetPtrFloat("THINjetPhoEF");

			// muon var
	      		Int_t nMu =inner_data.GetInt("nMu");
	      		TClonesArray* muP4 = (TClonesArray*) inner_data.GetPtrTObject("muP4");
	      		
	    		
	    		vector<bool> &isGlobalMuon  = *((vector<bool>*) inner_data.GetPtr("isGlobalMuon"));
	    		vector<bool> &isTrackerMuon  = *((vector<bool>*) inner_data.GetPtr("isTrackerMuon"));
	    		vector<bool> &isTightMuon  = *((vector<bool>*) inner_data.GetPtr("isTightMuon"));
	    		vector<bool> &isMediumMuon  = *((vector<bool>*) inner_data.GetPtr("isMediumMuon"));
	    		vector<bool> &isLooseMuon = *((vector<bool>*) inner_data.GetPtr("isLooseMuon"));
	    		vector<bool> &isHighPtMuon  = *((vector<bool>*) inner_data.GetPtr("isHighPtMuon"));
	    		vector<bool> &isSoftMuon  = *((vector<bool>*) inner_data.GetPtr("isSoftMuon"));
	    		
	    		Float_t* muChHadIso = inner_data.GetPtrFloat("muChHadIso");
	    		Float_t* muNeHadIso = inner_data.GetPtrFloat("muNeHadIso");
	    		Float_t* muGamIso = inner_data.GetPtrFloat("muGamIso");
	    		Float_t* muPUPt = inner_data.GetPtrFloat("muPUPt");

	      		// ele var
	      		Int_t nEle =inner_data.GetInt("nEle");
	      		TClonesArray* eleP4 = (TClonesArray*) inner_data.GetPtrTObject("eleP4");
	      		
	      		vector<bool> &eleIsPassTight  = *((vector<bool>*) inner_data.GetPtr("eleIsPassTight"));
	      		vector<bool> &eleIsPassMedium  = *((vector<bool>*) inner_data.GetPtr("eleIsPassMedium"));
	      		vector<bool> &eleIsPassLoose  = *((vector<bool>*) inner_data.GetPtr("eleIsPassLoose"));


			//met
			Float_t pfMetCorrPt=inner_data.GetFloat("pfMetCorrPt");
			Float_t pfMetCorrPhi=inner_data.GetFloat("pfMetCorrPhi");

			for (Int_t i=0;i<string_size; i++)
			{
				string this_trig_name = hlt_trigName[i];
				bool this_trig_result = hlt_trigResult[i];
				if (this_trig_name.find("HLT_Ele105")!=std::string::npos && this_trig_result==1)
				{ this_ev.is_pass_ele_trigger = true;}
				if (this_trig_name.find("HLT_Mu45")!=std::string::npos && this_trig_result==1)
				{ this_ev.is_pass_mu_trigger = true;}
				
			}

	      		for (Int_t i=0;i<nEle;i++)
			{	
				ele ele_i( *((TLorentzVector*)eleP4->At(i)),eleIsPassTight[i],eleIsPassMedium[i],eleIsPassLoose[i]);
				if (ele_i.p4.Pt()==0){cout<<"e";}
				this_ev.eles.push_back(ele_i);
			}

	      		for (Int_t i=0;i<nMu;i++)
			{
	    		
	    			mu mu_i( *((TLorentzVector*)muP4->At(i)), isGlobalMuon[i], isTrackerMuon[i], isTightMuon[i], isMediumMuon[i], isLooseMuon[i],isHighPtMuon[i], isSoftMuon[i], muChHadIso[i], muNeHadIso[i], muGamIso[i], muPUPt[i]);
				if (mu_i.p4.Pt()==0){cout<<"m";}
				this_ev.mus.push_back(mu_i);
			}
	      	      				
	      		for (Int_t i=0;i<THINnJet;i++)
			{	
				if (THINjetCISVV2[i]<0 || THINjetCISVV2[i]>1){continue;}
		
				THINjet jet_i( *((TLorentzVector*)THINjetP4->At(i)), THINjetCSV[i], THINjetCISVV2[i], THINjetHadronFlavor[i], THINjetPassIDTight[i],THINjetPassIDLoose[i], THINjetCEmEF[i], THINjetCHadEF[i], THINjetNEmEF[i], THINjetNHadEF[i], THINjetPhoEF[i]);
				if (jet_i.p4.Pt()==0){cout<<"j";}
				this_ev.THINjets.push_back(jet_i);

			}

met met_i(pfMetCorrPt,pfMetCorrPhi);
			this_ev.missing_et=met_i;

			// get good ele
			for (int i=0;i<this_ev.eles.size();i++)
			{
				ele p = this_ev.eles[i];
				if ( !is_pass_ele_loose_id(p)){continue;}
				this_ev.good_eles.push_back(p);
	
			}
			
			// get good mu	
			for (int i=0;i<this_ev.mus.size();i++)
			{
	
				mu p = this_ev.mus[i];
				if ( !is_pass_mu_idiso(p) ){continue;}
				this_ev.good_mus.push_back(p);	
			}
			
			// get THINjet after deoverlap with other sub jet and good lepton in every event
	     		for (Int_t i=0;i<this_ev.THINjets.size();i++)
			{
				THINjet this_j = this_ev.THINjets[i];


				bool overlap_with_lepton=false;
				/*
				for (int elei=0;elei<this_ev.good_eles.size();elei++)
				{	
					ele this_ele = this_ev.good_eles[elei];
					if (delta_R(this_j, this_ele)<0.4){overlap_with_lepton=true;break;}
				}
				*/
				for (int mui=0;mui<this_ev.good_mus.size();mui++)
				{
					mu this_mu = this_ev.good_mus[mui];
					if (delta_R(this_j, this_mu)<0.4){overlap_with_lepton=true;break;}
				}
				
				if (overlap_with_lepton){ continue;}
			
				this_ev.THINjets_deoverlap_with_lepton.push_back(this_j);
				
			}

			if (!is_pass_ev_selection_TTbar(this_ev)){continue;}




			// =================================================event passed event selection=========================================

	    		double scaled_value = get_scaled_value(Nth);
	    		double event_weight=0;
	    		Float_t mcWeight=inner_data.GetFloat("mcWeight");
	    		if (mcWeight>0){event_weight=1;}
	    		else {event_weight=-1;}
	    		
	    		
	    		
	    		if(Nth>=Nev_index_HT100to200)  // for HT sample
	    		{	
	    			N_available_event_for_HT+=scaled_value;		
				Float_t HT = inner_data.GetFloat("HT");
				h_HT_HT->Fill(HT, scaled_value);
					
	    		}
			else if(Nth>=Nev_index_inclusive) // for inclusive sample( used event weight)
			{	
				N_available_event_for_inclusive+=scaled_value*event_weight;		
				Float_t HT = inner_data.GetFloat("HT");
				h_inclusive_HT->Fill(HT, scaled_value*event_weight);
				
				if (HT>=100){Nev_HT100up_scaled +=scaled_value*event_weight;}
				else {Nev_HT100down_scaled+=scaled_value*event_weight; }
			}
			
			

		
		}
		

		save_hist(h_inclusive_HT, "pic_HT_comp", "inclusive_HT");
		save_hist(h_HT_HT, "pic_HT_comp", "HT_HT");
		
		cout<<"===================================result================================================="<<endl;
		cout<<"N_available_event_for_inclusive="<<N_available_event_for_inclusive<<endl;
		cout<<"N_available_event_for_HT="<<N_available_event_for_HT<<endl;
		cout<<endl;
		cout<<"Nev_HT100up_scaled="<<Nev_HT100up_scaled<<endl;
		cout<<"Nev_HT100down_scaled="<<Nev_HT100down_scaled<<endl;
		cout<<"Nev_HT100up_scaled + Nev_HT100down_scaled="<<(Nev_HT100up_scaled + Nev_HT100down_scaled)<<endl;
		cout<<"fraction="<<Nev_HT100up_scaled/(Nev_HT100up_scaled + Nev_HT100down_scaled)<<endl;






}
