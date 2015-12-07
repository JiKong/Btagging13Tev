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
  	//gereral
  	Long64_t Nev_MC=0;
  	
  	Long_t Nev_ttbar =0;
	Long_t Nev_ST_s_channel =0;
	Long_t Nev_ST_t_channel_antitop =0;
	Long_t Nev_ST_t_channel_top =0;
	Long_t Nev_ST_tW_antitop =0;
	Long_t Nev_ST_tW_top =0;	  	  	
	Long_t Nev_DY_inclusive =0;
	Long_t Nev_QCD100to200 =0;
	Long_t Nev_QCD200to300 =0;
	Long_t Nev_QCD300to500 =0;
	Long_t Nev_QCD500to700 =0;
	Long_t Nev_QCD700to1000 =0;
	Long_t Nev_QCD1000to1500 =0;
	Long_t Nev_QCD1500to2000 =0;
	Long_t Nev_QCD2000toInf =0;
	//Long_t Nev_WJ =0;
	Long_t Nev_WJ100to200 =0;
	Long_t Nev_WJ200to400 =0;
	Long_t Nev_WJ400to600 =0;
	Long_t Nev_WJ600to800 =0;
	Long_t Nev_WJ800to1200 =0;
	Long_t Nev_WJ1200to2500 =0;
	Long_t Nev_WW =0;
	Long_t Nev_WZ =0;
	Long_t Nev_ZZ =0;
		
	double Nev_MC_scaled=0;
	double N_available_event_MC_scaled=0;

	double N_available_event_ttbar_scaled =0;
	double N_available_event_ST_s_channel_scaled =0;
	double N_available_event_ST_t_channel_antitop_scaled =0;
	double N_available_event_ST_t_channel_top_scaled =0;
	double N_available_event_ST_tW_antitop_scaled =0;
	double N_available_event_ST_tW_top_scaled =0;
	double N_available_event_DY_inclusive_scaled =0;
	double N_available_event_QCD100to200_scaled =0;
	double N_available_event_QCD200to300_scaled =0;
	double N_available_event_QCD300to500_scaled =0;
	double N_available_event_QCD500to700_scaled =0;
	double N_available_event_QCD700to1000_scaled =0;
	double N_available_event_QCD1000to1500_scaled =0;
	double N_available_event_QCD1500to2000_scaled =0;
	double N_available_event_QCD2000toInf_scaled =0;
	//double N_available_event_WJ_scaled =0;
	double N_available_event_WJ100to200_scaled =0;
	double N_available_event_WJ200to400_scaled =0;
	double N_available_event_WJ400to600_scaled =0;
	double N_available_event_WJ600to800_scaled =0;
	double N_available_event_WJ800to1200_scaled =0;
	double N_available_event_WJ1200to2500_scaled =0;
	double N_available_event_WW_scaled =0;
	double N_available_event_WZ_scaled =0;
	double N_available_event_ZZ_scaled =0;
	
	
	int Nev_index_ttbar=0;
	int Nev_index_ST_s_channel =0;
	int Nev_index_ST_t_channel_antitop =0;
	int Nev_index_ST_t_channel_top =0;
	int Nev_index_ST_tW_antitop =0;
	int Nev_index_ST_tW_top =0;
	int Nev_index_DY_inclusive =0;		
	int Nev_index_QCD100to200 =0;
	int Nev_index_QCD200to300 =0;
	int Nev_index_QCD300to500 =0;
	int Nev_index_QCD500to700 =0;		
	int Nev_index_QCD700to1000 =0;
	int Nev_index_QCD1000to1500 =0;
	int Nev_index_QCD1500to2000 =0;
	int Nev_index_QCD2000toInf =0;
	//int Nev_index_WJ =0;
	int Nev_index_WJ100to200 =0;
	int Nev_index_WJ200to400 =0;
	int Nev_index_WJ400to600 =0;
	int Nev_index_WJ600to800 =0;
	int Nev_index_WJ800to1200 =0;
	int Nev_index_WJ1200to2500 =0;
	int Nev_index_WW =0;
	int Nev_index_WZ =0;
	int Nev_index_ZZ =0;

	// gen level 
	vector< vector<double> > fraction_ratio;
	int max_true_b=0;
	int max_true_nb=0;
	
	double true_eff=0;
	double true_mr=0;

	double gen_j=0;
	double gen_b=0;
	double gen_nb=0;
	double gen_j_pass_cut=0;
	double gen_b_pass_cut=0;
	double gen_nb_pass_cut=0;
	double gen_j_pass_cut_btagging=0;
	double gen_b_pass_cut_btagging=0;
	double gen_nb_pass_cut_btagging=0;

	// rec level
	vector<double> Nt_from_rec(1,0) ; // # of event that have t b-tagging
	vector<double> Nk_from_rec(1,0) ; // # of event that have k jet

	double best_eff_rec=0;
	double best_mr_rec=0;

	double total_kk_rec=0;
	double total_k_rec=0;
	double total_btag_pass_cut_rec=0;
	double total_nbtag_pass_cut_rec=0;
	

	// data_lumi_13_Tev
	//double data_lumi_13_Tev = 2297.21;
	double data_lumi_13_Tev = 831.7; //  v3,v4 in Raman's sample. 
	//double data_lumi_13_Tev = 370.2; // half root file on muon channel
	//double data_lumi_13_Tev = 379.414;  // half root file on electron channel
	
	// MC sample X-section(only central value):
	double Xs_ttbar = 831.76;
	double Xs_ST_s_channel=3.38;
	double Xs_ST_t_channel_antitop=26.49;
	double Xs_ST_t_channel_top=44.51;
	double Xs_ST_tW_antitop=35.6;
	double Xs_ST_tW_top=35.6;
	double Xs_DY_inclusive=6025.2;		
	double Xs_QCD100to200=27850000;
	double Xs_QCD200to300=1717000;
	double Xs_QCD300to500=351300;
	double Xs_QCD500to700=31630;			
	double Xs_QCD700to1000=6524;
	double Xs_QCD1000to1500=1064;
	double Xs_QCD1500to2000=121.5;
	double Xs_QCD2000toInf=25.42;
	//double Xs_WJ=61526;
	double Xs_WJ100to200 =1627.45;
	double Xs_WJ200to400 =435.237;
	double Xs_WJ400to600 =59.18;
	double Xs_WJ600to800 =14.58;
	double Xs_WJ800to1200 =6.65;
	double Xs_WJ1200to2500 =1.6;
	double Xs_WW=118.7;
	double Xs_WZ=47.13;
	double Xs_ZZ = 16.523;
	
	
	// MC sample scale factor:
	
	double sf_ttbar = 0;
	double sf_ST_s_channel=0;
	double sf_ST_t_channel_antitop=0;
	double sf_ST_t_channel_top=0;
	double sf_ST_tW_antitop=0;
	double sf_ST_tW_top=0;
	double sf_DY_inclusive=0;
	double sf_QCD100to200=0;
	double sf_QCD200to300=0;
	double sf_QCD300to500=0;
	double sf_QCD500to700=0;
	double sf_QCD700to1000=0;
	double sf_QCD1000to1500=0;
	double sf_QCD1500to2000=0;
	double sf_QCD2000toInf=0;
	//double sf_WJ=0;
	double sf_WJ100to200 =0;
	double sf_WJ200to400 =0;
	double sf_WJ400to600 =0;
	double sf_WJ600to800 =0;
	double sf_WJ800to1200 =0;
	double sf_WJ1200to2500 =0;
	double sf_WW=0;
	double sf_WZ=0;
	double sf_ZZ =0;
	

// local class
// stack_hist

class stack_hist
{
	public:
	THStack *hs ;
	TH1D* for_ttbar;
	TH1D* for_ST;
	TH1D* for_DY;
	TH1D* for_QCD;
	TH1D* for_WJ;
	TH1D* for_DB;
	stack_hist(TString name_,TString title_, int N_bin_, double min_,double max_)
	{
		hs = new THStack(name_,title_);
		
		for_ttbar = new TH1D(name_+"_ttbar",title_, N_bin_, min_,max_);
		for_ST = new TH1D(name_+"_ST",title_, N_bin_, min_,max_);
		for_DY = new TH1D(name_+"_DY",title_, N_bin_, min_,max_);
		for_QCD = new TH1D(name_+"_QCD",title_, N_bin_, min_,max_);
		for_WJ = new TH1D(name_+"_WJ",title_, N_bin_, min_,max_);
		for_DB = new TH1D(name_+"_DB",title_, N_bin_, min_,max_);		
	}

	void save_sh( TString dirname_, TString filename_) 
	{
		for_ttbar->SetFillColor(2);
		for_ST->SetFillColor(3);
		for_DY->SetFillColor(4);
		for_QCD->SetFillColor(5);
		for_WJ->SetFillColor(6);
		for_DB->SetFillColor(7);
		hs->Add(for_ttbar); hs->Add(for_ST); hs->Add(for_DY); hs->Add(for_QCD); hs->Add(for_WJ); hs->Add(for_DB);
		   
		TCanvas *c1 = new TCanvas(filename_,filename_,200,10,700,500);

		//c1->SetFillColor(42);
		c1->SetGrid();
		c1->GetFrame()->SetFillColor(21);
		c1->GetFrame()->SetBorderSize(12);

		hs->Draw();  

		TLegend* leg = new TLegend(0.8,0.8,0.99,0.99);
		leg->AddEntry(for_DB,"DiBoson","f");   // "f": f means color box
		leg->AddEntry(for_WJ,"W+jet","f");
		leg->AddEntry(for_QCD,"QCD","f");
		leg->AddEntry(for_DY,"DY","f");
		leg->AddEntry(for_ST,"Single Top","f");
		leg->AddEntry(for_ttbar,"ttbar","f");
		
		leg->Draw();
		
		TImage *img1 = TImage::Create();
		img1->FromPad(c1);
		img1->WriteImage("./"+dirname_+"/"+filename_+".png");
		delete c1;
		delete img1;
		
		//save hist in root file 
		TFile* this_rf=new TFile("./"+dirname_+"_rf/"+filename_+".root", "recreate");
		for_ttbar->Write();
		for_ST->Write();
		for_DY->Write();
		for_QCD->Write();
		for_WJ->Write();
		for_DB->Write();
		this_rf->Close();
		

	}

	
	void fill_hist(double value_, double weight_, int Nth_)
	{

		if(Nth_>=Nev_index_WW) {for_DB->Fill(value_, weight_);	 }
		//else if(Nth_>=Nev_index_WJ) {for_WJ->Fill(value_, weight_);}
		else if(Nth_>=Nev_index_WJ100to200) {for_WJ->Fill(value_, weight_);}
		else if(Nth_>=Nev_index_QCD100to200) {for_QCD->Fill(value_, weight_); }
		else if(Nth_>=Nev_index_DY_inclusive) {for_DY->Fill(value_, weight_); }
		else if(Nth_>=Nev_index_ST_s_channel) {for_ST->Fill(value_, weight_);}
		else if(Nth_>=Nev_index_ttbar) {for_ttbar->Fill(value_, weight_); }
	}

	
};



//function======================================================================================================

//////////====tt
void tt(double i)
{
	cout<<"CHK"<<i<<endl;
}

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



//////////====math
Int_t factorial(Int_t n_)
{
	if (n_<0){return 0;}
	Int_t r=1;
	Int_t _n_=n_;
	for (Int_t i=1;i<=_n_;i++)
	{
		r=r*i;
	}
	return r;
}

long C(Int_t n,Int_t m)
{
	if (m>n || m<0){return 0;}
	long r_top=1;
	long r_down=1;	
	Int_t script;	
	if((n-m)<(m)){script=n-m;}
	else {script=m;}
	for (Int_t i =0;i<script;i++)
	{ r_top=r_top*(n-i); r_down=r_down*(i+1);}
	return (r_top/r_down);
//	return (factorial(n)/(factorial(m)*factorial(n-m))); //n>=18, memory overflow
}



double log_poisson(double x, double lamda)
{
	// it woks when x=200,000,000
	double term1 = (-1) * lamda;
	double term2 = x * log(lamda);
	double term3 = 0;
	for (double i = 0; i < x; i++)
	{
		term3 = term3 + log(x - i);
	}
	term3 = (-1) * term3;

	return term1 + term2 + term3;
}

// draw histo, graph
void gerrors(sample_result r, string type) // type: eff mr
{
	
   	//Draw a graph with error bars
	// To see the output of this macro, click begin_html <a href="gif/gerrors.gif">here</a>. end_html
	//Author: Rene Brun
	   
	TCanvas *c1 = new TCanvas(r.name,r.name,200,10,700,500);

	c1->SetFillColor(42);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor(21);
	c1->GetFrame()->SetBorderSize(12);

 	Int_t n = r.N_sample;
	Float_t x[n];
	Float_t ex[n];
	
	Float_t y_rec[n];
	Float_t ey_rec[n];
	
	Float_t y_gen[n];
	Float_t ey_gen[n];
	
	for (int i=0;i<n;i++)
	{
		if (i<=5)
		{ x[i]=1000+200*i;}
		else 
		{ x[i]=2000+500*(i-5);}
	   	
	   	ex[i]=0;
	   	
	   	
	   	if (type == "eff")
	   	{
	   		y_rec[i]=r.rec_effs[i];
	   		ey_rec[i]=r.rec_effs_err[i];
	   		
	   		y_gen[i]=r.gen_effs[i];
	   		ey_gen[i]=0;
	   	}
	     	else //if (type == "mr")
	   	{
	   		y_rec[i]=r.rec_mrs[i];
	   		ey_rec[i]=r.rec_mrs_err[i];
	   		
	   		y_gen[i]=r.gen_mrs[i];
	   		ey_gen[i]=0;
	   	}

	}  
	
	TGraphErrors *grec = new TGraphErrors(n,x,y_rec,ex,ey_rec);
	grec->SetName(r.name);
	grec->SetTitle(r.name);
	grec->SetMarkerColor(4);
	grec->SetLineColor(4);
	grec->SetMarkerStyle(21);
	grec->GetXaxis()->SetTitle("Zprime Mass");
	grec->Draw("AC");  // overlay
	
	TGraphErrors *ggen = new TGraphErrors(n,x,y_gen,ex,ey_gen);
	ggen->SetName(r.name);
	ggen->SetTitle(r.name);
	ggen->SetMarkerColor(2);
	ggen->SetMarkerStyle(21);
	ggen->Draw("CP");  //overlay
	
 
	TImage *img1 = TImage::Create();
	img1->FromPad(c1);
	const char* filename = ("./Zprime_sample/"+type+".png").c_str(); //string to const char*
	img1->WriteImage(filename);
	delete c1;
	delete img1;


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
	
	//save hist in root file 
	TFile* this_rf=new TFile("./"+dirname_+"_rf/"+filename_+".root", "recreate");
	h_->Write();
	this_rf->Close();

}



// general
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


//for gen par
//for Zprime
bool is_wanted_b_Zprime(gen_par g_) //in 13 Tev, the requirement that "g_.status==3" is a problem
{
	if (abs(g_.id)==5 /*&& g_.status==3*/ && abs(g_.mid)==25)
	{ return true;}
	else
	{ return false;}
}


// for TTbar sample
bool is_wanted_b_ttbar(gen_par g_)
{
	if ((g_.id==5 && g_.mid==6) ||  (g_.id==-5 && g_.mid==-6)   )
	{ return true;}
	else
	{ return false;}
}

bool is_wanted_b_allBKG(gen_par g_)
{
	if (abs(g_.id)==5 )
	{ return true;}
	else
	{ return false;}
}

bool is_gen_j(gen_par g_)
{
	if (abs(g_.id)<=6 && abs(g_.id)>=1) //no id==7(b'), id==8(t') found.
	{ return true;}
	else
	{ return false;}
}


// jet selection
bool is_pass_csvl(jet j_)
{
	if (j_.csv>0.244 && j_.csv<=1){return true;}
	else {return false;}
}



bool is_pass_csvm(jet j_)
{
	if (j_.csv>0.679 && j_.csv<=1){return true;}
	else {return false;}
}

bool is_pass_cisvv2_13Tev(jet j_)
{	
	if (j_.cisvv2>0.605 && j_.cisvv2<=1){return true;}
	else {return false;}

}

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

/*
bool is_pass_pteta(PO O1_)
{
	if(O1_.p4.Pt()<20){return false;}
	if(abs(O1_.p4.Eta())>2.4){return false;}
	return true;

}*/


//ele selection

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

//muon selection  

bool is_pass_mu_idiso(mu m_)
{
	//if (!m_.is_tight==true){return false;}
	if (!m_.is_loose==true){return false;}
	double this_v = (m_.ChHadIso+max(0., m_.NeHadIso + m_.GamIso - 0.5*m_.PUPt))/m_.p4.Pt(); 
	if (! (this_v<0.15) ){return false;}
	return true;

}



//event selection
//evse gen level



//evse rec level
bool is_pass_ev_selection_TTbar(event ev)
{

	// trigger path
	if (!ev.is_pass_mu_trigger){return false;}
	
	// same as theevent selection on ttbar sample, be cause we will see how ttbar selection works.....
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
	
	/*
	ev_flag:
		=0: this event is good event
		=2: !(have good lepton), skip this event.
		=3: !(have at least 2 good jets), skip this event.
		=4: !(have 2 good leading jets), skip this event.
		=5: !(missing Et > 20), skip this event.
	*/
}


double get_Tmass(event ev)
{
	vector<THINjet> this_THINjets = ev.THINjets_deoverlap_with_lepton;

	bool find_good_Wqqbar=false;
	int q1_index=-1;
	int q2_index=-1;
	double rec_Wmass=0;
	for (int i=0;i<this_THINjets.size();i++)
	{
		for (int j=0;j<i;j++)
		{
			if (!is_pass_cisvv2_13Tev(this_THINjets[i]) && !is_pass_cisvv2_13Tev(this_THINjets[j]) )
			{
				TLorentzVector q1l= this_THINjets[i].p4;
				TLorentzVector q2l= this_THINjets[j].p4;
				TLorentzVector Wl= q1l+q2l;
				double this_Wmass = Wl.M();
				if (abs(this_Wmass-81) <=20){find_good_Wqqbar=true;}
				if (abs(this_Wmass-81) <= abs(rec_Wmass-81))
				{
					rec_Wmass=this_Wmass;
					q1_index=i;
					q2_index=j;					
				}  
			}
		}
	}

	if (!find_good_Wqqbar){return -1;}

	bool find_good_b=false;
	double rec_Tmass=-1;
	for (int i=0;i<this_THINjets.size();i++)
	{
		if (is_pass_cisvv2_13Tev(this_THINjets[i])) 
		{
			find_good_b=true;
			TLorentzVector bq= this_THINjets[i].p4;
			TLorentzVector q1l= this_THINjets[q1_index].p4;
			TLorentzVector q2l= this_THINjets[q2_index].p4;
			TLorentzVector Tl= bq+q1l+q2l;
			double this_Tmass = Tl.M();
			if (abs(this_Tmass-174) <= abs(rec_Tmass-174))
			{
				rec_Tmass=this_Tmass;				
			}  			 
		}
	}

	if (!find_good_b){return -1;}
	else {return rec_Tmass;}


}


double get_scaled_value(Long_t Nth_)
{
	double scaled_value=0;
	if (Nth_>=Nev_index_ZZ) {scaled_value=sf_ZZ;}
	else if(Nth_>=Nev_index_WZ) {scaled_value=sf_WZ;}
	else if(Nth_>=Nev_index_WW) {scaled_value=sf_WW; }
	//else if (Nth_>=Nev_index_WJ) {scaled_value=sf_WJ;}
	else if (Nth_>=Nev_index_WJ1200to2500) {scaled_value=sf_WJ1200to2500;}
	else if(Nth_>=Nev_index_WJ800to1200) {scaled_value=sf_WJ800to1200;}
	else if(Nth_>=Nev_index_WJ600to800) {scaled_value=sf_WJ600to800;}
	else if(Nth_>=Nev_index_WJ400to600) {scaled_value=sf_WJ400to600;}
	else if(Nth_>=Nev_index_WJ200to400) {scaled_value=sf_WJ200to400;}
	else if(Nth_>=Nev_index_WJ100to200) {scaled_value=sf_WJ100to200;}
	else if(Nth_>=Nev_index_QCD2000toInf) {scaled_value=sf_QCD2000toInf; }
	else if(Nth_>=Nev_index_QCD1500to2000) {scaled_value=sf_QCD1500to2000; }
	else if(Nth_>=Nev_index_QCD1000to1500) {scaled_value=sf_QCD1000to1500; }
	else if(Nth_>=Nev_index_QCD700to1000) {scaled_value=sf_QCD700to1000; }
	else if(Nth_>=Nev_index_QCD500to700) {scaled_value=sf_QCD500to700; }
	else if(Nth_>=Nev_index_QCD300to500) {scaled_value=sf_QCD300to500; }
	else if(Nth_>=Nev_index_QCD200to300) {scaled_value=sf_QCD200to300; }
	else if(Nth_>=Nev_index_QCD100to200) {scaled_value=sf_QCD100to200; }
	else if(Nth_>=Nev_index_DY_inclusive) {scaled_value=sf_DY_inclusive; }
	else if(Nth_>=Nev_index_ST_tW_top) {scaled_value=sf_ST_tW_top; }
	else if(Nth_>=Nev_index_ST_tW_antitop) {scaled_value=sf_ST_tW_antitop; }
	else if(Nth_>=Nev_index_ST_t_channel_top) {scaled_value=sf_ST_t_channel_top; }
	else if(Nth_>=Nev_index_ST_t_channel_antitop) {scaled_value=sf_ST_t_channel_antitop; }
	else if(Nth_>=Nev_index_ST_s_channel) {scaled_value=sf_ST_s_channel;}
	else if(Nth_>=Nev_index_ttbar) {scaled_value=sf_ttbar; }
	
	return scaled_value;
}


Float_t get_event_weight(Long_t Nth_, Float_t mcWeight_)
{
	if (Nth_>=Nev_index_ZZ) {return 1;}
	else if(Nth_>=Nev_index_WZ) {return 1; }
	else if(Nth_>=Nev_index_WW) {return 1;}
	/*
	else if(Nth_>=Nev_index_WJ)  // amcatNlo sample
	{
		if (mcWeight_>0){return 1;}
		else {return-1;}
	}*/
	else if(Nth_>=Nev_index_WJ1200to2500) {return 1; }
	else if(Nth_>=Nev_index_WJ800to1200) {return 1;}
	else if(Nth_>=Nev_index_WJ600to800) {return 1; }
	else if(Nth_>=Nev_index_WJ400to600) {return 1; }
	else if(Nth_>=Nev_index_WJ200to400) {return 1; }
	else if(Nth_>=Nev_index_WJ100to200) {return 1; }		
	else if(Nth_>=Nev_index_QCD2000toInf) {return 1; }
	else if(Nth_>=Nev_index_QCD1500to2000) {return 1;}
	else if(Nth_>=Nev_index_QCD1000to1500) {return 1; }
	else if(Nth_>=Nev_index_QCD700to1000) {return 1; }
	else if(Nth_>=Nev_index_QCD500to700) {return 1; }
	else if(Nth_>=Nev_index_QCD300to500) {return 1;}
	else if(Nth_>=Nev_index_QCD200to300) {return 1; }
	else if(Nth_>=Nev_index_QCD100to200) {return 1; }			
	else if(Nth_>=Nev_index_DY_inclusive)   // amcatNlo sample
	{
		if (mcWeight_>0){return 1;}
		else {return-1;}
	}
	else if(Nth_>=Nev_index_ST_tW_top) {return 1; }
	else if(Nth_>=Nev_index_ST_tW_antitop) {return 1; }
	else if(Nth_>=Nev_index_ST_t_channel_top) {return 1;}
	else if(Nth_>=Nev_index_ST_t_channel_antitop) {return 1;}
	else if(Nth_>=Nev_index_ST_s_channel)   // amcatNlo sample
	{
		if (mcWeight_>0){return 1;}
		else {return-1;}
	}
	else if(Nth_>=Nev_index_ttbar) {return 1; }
}



// calculate likelihood



// calculate expection value of Nt
double get_expvalue_of_Nt(  double eff, double mr, int t)
{
        double this_expvalue = 0;
        for (int i = 0; i < fraction_ratio.size() ; i++)
        {
                for (int j = 0; j < fraction_ratio[i].size() ; j++)
                {
                    	double term1 = fraction_ratio[i][j]<=0 ? 0 : N_available_event_MC_scaled * fraction_ratio[i][j];

                        for (int ii = 0; ii <= i; ii++) // in i b-jet, ii b-jet pass btagging
                        {
                                for (int jj = 0; jj <= j; jj++) // in j nb-jet, jj nb-jet pass btagging
                                {
                                    	if (ii + jj != t) { continue; }

                                    	double term2 = C(i, ii) * pow(eff, ii) * pow(1 - eff, i - ii);
                                    	double term3 = C(j, jj) * pow(mr, jj) * pow(1 - mr, j - jj);
                                    	this_expvalue += (term1 *  term2 * term3);

                                    	//debug
                                    	if (term1 * term2 * term3 < 0)
                                    	{
						cout<<"in get_expvalue_of_Nt"<<endl;
                              			cout<<"<"<<term1 * term2 * term3<<". "<<"term1="<<term1<<",term2="<<term2<<",term3="<<term3<<endl;
                                		cout<<"eff="<<eff<<","<<endl;
                                		cout<<"mr="<<mr<<","<<endl;
                                		cout<<"t="<<t<<","<<endl;
                                		cout<<"N_available_event_MC_scaled="<<N_available_event_MC_scaled<<","<<endl;
                                		cout<<"fraction_ratio[i][j]="<<fraction_ratio[i][j]<<","<<endl;
                                		cout<<"i="<<i<<", "<<"j="<<j<<endl;
                                    	}
                                }
                        }

                }
        }
        return this_expvalue;
}



//for minuit

void fcn_btagging_rec(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{

	double total_log_like=0;
        for (int t = 0; t < Nt_from_rec.size(); t++)
        {
		// par[0] is eff, par[1] is mr
		if(get_expvalue_of_Nt( par[0], par[1] ,t)==0){continue;}
		double this_log_like=log_poisson(Nt_from_rec[t], get_expvalue_of_Nt( par[0], par[1] ,t));
                total_log_like += this_log_like;
		if (this_log_like < -10000000)
		{ 
			cout<<"this_log_like < -10000000"<<"when eff="<<par[0]<<",mr="<<par[1]<<",t="<<t;
			cout<<",Nt_from_rec[t]="<<Nt_from_rec[t]<<",get_expvalue_of_N="<<get_expvalue_of_Nt( par[0], par[1] ,t)<<endl;
		}  
        }

	f=-2*total_log_like;

}



vector<par_err_pair_for_minuit> get_min_par_btagging2p_rec(vector<par_for_minuit> pars_ )
{
	int Npar = pars_.size();
  	TMinuit *gMinuit = new TMinuit(Npar);  //initialize TMinuit with a maximum of 2 params
  	gMinuit->SetFCN(fcn_btagging_rec); // fcn_btagging is the function to do minimization

  	Double_t arglist[10];
  	Int_t ierflg = 0; // status flag, it will be 0 when no programing error found

  	// -- sets par error
  	arglist[0] = 1;
  	gMinuit->mnexcm("SET ERR", arglist ,Npar,ierflg);
	
	//set pars
	vector<par_err_pair_for_minuit> par_err_pairs;
	for(int i=0;i<Npar;i++)
	{
		Double_t number = pars_[i].number;
		string name = pars_[i].name;
	 	Double_t start = pars_[i].fit_start;	
  		Double_t step = pars_[i].fit_step;
		Double_t min = pars_[i].min;
		Double_t max = pars_[i].max;
  		gMinuit->mnparm(number, name, start, step, min,max,ierflg);
		
		par_err_pair_for_minuit A(0,0); 
		par_err_pairs.push_back(A);	
	}

 	// Now ready for minimization step
 	arglist[0] = 1000; //number of iteration, beyond this number, Tminuit stop. 
  	arglist[1] = 0.01;
  	gMinuit->mnexcm("MIGRAD", arglist ,Npar,ierflg);

	for (int i=0;i<Npar;i++)
	{
  		if ( ierflg == 0 ) 
   		{
      			gMinuit->GetParameter(i,par_err_pairs[i].v , par_err_pairs[i].err );  	
		}
		else
		{ cout<<"get_minuit have the error code:"<<ierflg<<endl; }
	}
	return par_err_pairs;
}



//main =======================================================================================================================================================
void Btagging_in_13Tev_MC()
{

	TFile output("output.root","recreate");

	//hist
	//TH1D* var = new TH1D("name","title",bin,xmin,xmax);	
	TH1D* H_IM = new TH1D("H_IM","H_IM",200,0,200);
	TH1D* H_IM_csvl = new TH1D("H_IM_csvl","H_IM_csvl",200,0,200);
	TH1D* H_IM_csvm = new TH1D("H_IM_csvm","H_IM_csvm",200,0,200);

	TH1D* Zee_IM = new TH1D("Zee_IM","Zee_IM",200,0,200);




	// is_pass_ev_selection function
	typedef bool (*ev_se)(event);
	//ev_se is_pass_ev_selection = is_pass_ev_selection_Zprime_signal;
	ev_se is_pass_ev_selection = is_pass_ev_selection_TTbar;	
	// Btagging function
	typedef bool (*jet_alg)(jet);
	//jet_alg btagging = is_pass_csvl;
	//jet_alg btagging = is_pass_csvm;
	jet_alg btagging = is_pass_cisvv2_13Tev;
	// is wanted b
	typedef bool (*chk_gen_b)(gen_par);
	chk_gen_b is_wanted_b = is_wanted_b_allBKG;
	//chk_gen_b is_wanted_b = is_wanted_b_TTbar;
	//chk_gen_b is_wanted_b = is_wanted_b_Zprime;

		
		
	
	
	sample_result zprime_sample_result("zprime_sample_result");
	
	for (int sample_i=0;sample_i<1; sample_i++) //collect all small samples in a big sample
	{
		// initialization variables
		N_available_event_MC_scaled=0;

		// gen level 
		max_true_b=0;
		max_true_nb=0;
	
		true_eff=0;
		true_mr=0;

		gen_j=0;
		gen_b=0;
		gen_nb=0;
		gen_j_pass_cut=0;
		gen_b_pass_cut=0;
		gen_nb_pass_cut=0;
		gen_j_pass_cut_btagging=0;
		gen_b_pass_cut_btagging=0;
		gen_nb_pass_cut_btagging=0;

		// rec level

		best_eff_rec=0;
		best_mr_rec=0;

		total_kk_rec=0;
		total_k_rec=0;
		total_btag_pass_cut_rec=0;
		total_nbtag_pass_cut_rec=0;
	


		//hist set

		
		// stack hist set	// recjet(THINjet)
		
		// ele var (insignificant var)
	
		stack_hist* sh_ele_pt_a_evcut = new stack_hist("ele_pt_a_evcut","ele_pt_a_evcut",50,0,300);		
		stack_hist* sh_ele_eta_a_evcut = new stack_hist("ele_eta_a_evcut","ele_eta_a_evcut",50,-7,7);
		stack_hist* sh_ele_dEtaAtVtx_a_evcut = new stack_hist("ele_dEtaAtVtx_a_evcut","ele_dEtaAtVtx_a_evcut",50,-0.2,0.2);
		stack_hist* sh_ele_dPhiAtVtx_a_evcut = new stack_hist("ele_dPhiAtVtx_a_evcut","ele_dPhiAtVtx_a_evcut",50,-1,1);
		stack_hist* sh_ele_D0_a_evcut = new stack_hist("ele_D0_a_evcut","ele_D0_a_evcut",50,-2,2);
		stack_hist* sh_ele_EtaseedAtVtx_a_evcut = new stack_hist("ele_EtaseedAtVtx_a_evcut","ele_EtaseedAtVtx_a_evcut",50,-0.2,0.2);
		stack_hist* sh_ele_HoverE_a_evcut = new stack_hist("ele_HoverE_a_evcut","ele_HoverE_a_evcut",50,0,6);
		stack_hist* sh_ele_MissHits_a_evcut = new stack_hist("ele_MissHits_a_evcut","ele_MissHits_a_evcut",7,0,7);
		stack_hist* sh_ele_MiniIso_a_evcut = new stack_hist("ele_MiniIso_a_evcut","ele_MiniIso_a_evcut",50,0,10);
		stack_hist* sh_ele_SigmaIEtaIEta_a_evcut = new stack_hist("ele_SigmaIEtaIEta_a_evcut","ele_SigmaIEtaIEta_a_evcut",50,0,0.08);
		
		//mu var (insignificant var)
	
		stack_hist* sh_mu_pt_a_evcut = new stack_hist("mu_pt_a_evcut","mu_pt_a_evcut",50,0,300);		
		stack_hist* sh_mu_eta_a_evcut = new stack_hist("mu_eta_a_evcut","mu_eta_a_evcut",50,-7,7);
		stack_hist* sh_mu_dxy_a_evcut = new stack_hist("mu_dxy_a_evcut","mu_dxy_a_evcut",20,-1,1);
		stack_hist* sh_mu_dz_a_evcut = new stack_hist("mu_dz_a_evcut","mu_dz_a_evcut",40,-2,2);		
		stack_hist* sh_mu_Hits_a_evcut = new stack_hist("mu_Hits_a_evcut","mu_Hits_a_evcut",50,0,100);
		stack_hist* sh_mu_Matches_a_evcut = new stack_hist("mu_Matches_a_evcut","mu_Matches_a_evcut",10,0,10);	
		stack_hist* sh_mu_MiniIso_a_evcut = new stack_hist("mu_MiniIso_a_evcut","mu_MiniIso_a_evcut",50,0,100);
		stack_hist* sh_mu_PixelHits_a_evcut = new stack_hist("mu_PixelHits_a_evcut","mu_PixelHits_a_evcut",20,0,20);
		stack_hist* sh_mu_TrkLayers_a_evcut = new stack_hist("mu_TrkLayers_a_evcut","mu_TrkLayers_a_evcut",25,0,25);
		stack_hist* sh_mu_TrkPt_a_evcut = new stack_hist("mu_TrkPt_a_evcut","mu_TrkPt_a_evcut",50,0,500);
		
	
		//THINjet var
				
		stack_hist* sh_THINjet_CEmEF_a_evcut = new stack_hist("THINjet_CEmEF_a_evcut","THINjet_CEmEF_a_evcut",50,0,1);
		stack_hist* sh_THINjet_CHadEF_a_evcut = new stack_hist("THINjet_CHadEF_a_evcut","THINjet_CHadEF_a_evcut",50,0,1);
		stack_hist* sh_THINjet_NEmEF_a_evcut = new stack_hist("THINjet_NEmEF_a_evcut","THINjet_NEmEF_a_evcut",50,0,1);
		stack_hist* sh_THINjet_NHadEF_a_evcut = new stack_hist("THINjet_NHadEF_a_evcut","THINjet_NHadEF_a_evcut",50,0,1);
		stack_hist* sh_THINjet_PhoEF_a_evcut = new stack_hist("THINjet_PhoEF_a_evcut","THINjet_PhoEF_a_evcut",50,0,1);
							
		stack_hist* sh_THINjet_cisvv2_a_evcut = new stack_hist("THINjet_cisvv2_a_evcut","THINjet_cisvv2_a_evcut",50,0,1);
		stack_hist* sh_THINjet_pt_a_evcut = new stack_hist("THINjet_pt_a_evcut","THINjet_pt_a_evcut",50,0,500);		
		stack_hist* sh_THINjet_eta_a_evcut = new stack_hist("THINjet_eta_a_evcut","THINjet_eta_a_evcut",50,-7,7);
		
		stack_hist* sh_THINBjet_CEmEF_a_evcut = new stack_hist("THINBjet_CEmEF_a_evcut","THINBjet_CEmEF_a_evcut",50,0,1);
		stack_hist* sh_THINBjet_CHadEF_a_evcut = new stack_hist("THINBjet_CHadEF_a_evcut","THINBjet_CHadEF_a_evcut",50,0,1);
		stack_hist* sh_THINBjet_NEmEF_a_evcut = new stack_hist("THINBjet_NEmEF_a_evcut","THINBjet_NEmEF_a_evcut",50,0,1);
		stack_hist* sh_THINBjet_NHadEF_a_evcut = new stack_hist("THINBjet_NHadEF_a_evcut","THINBjet_NHadEF_a_evcut",50,0,1);
		stack_hist* sh_THINBjet_PhoEF_a_evcut = new stack_hist("THINBjet_PhoEF_a_evcut","THINBjet_PhoEF_a_evcut",50,0,1);
		
		stack_hist* sh_THINBtagjet_cisvv2_a_evcut = new stack_hist("THINBtagjet_cisvv2_a_evcut","THINBtagjet_cisvv2_a_evcut",50,0,1);
		stack_hist* sh_THINBtagjet_pt_a_evcut = new stack_hist("THINBtagjet_pt_a_evcut","THINBtagjet_pt_a_evcut",50,0,500);		
		stack_hist* sh_THINBtagjet_eta_a_evcut = new stack_hist("THINBtagjet_eta_a_evcut","THINBtagjet_eta_a_evcut",50,-7,7);

		// rec t mass
		stack_hist* sh_rec_Tmass_a_evcut = new stack_hist("rec_Tmass_a_evcut","rec_Tmass_a_evcut",50,0,300);		
	
		// root file
	  	//TChain *MC = new TChain("tree/treeMaker");
	  	//MC->Add("../flattuple.root");// Zprime sample
	  	//MC->Add("../TTbar_sample/flattuple_1.root");// TTbar sample from yu
	  	//MC->Add(samples[sample_i]);
	  	
	  	vector<string> input_files;
	  	
		int Nev_index_counter=0;		

 		//key word:%5 
		int N_applied_sample_min=0;
		int N_applied_sample_max=20;

	  	// TTbar first 1000 sample
	  	bool is_save_Nev_index_ttbar =false;
	  	for (int i=1;i<=999;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			if (i==36 || i==70 ||i==304 || i==311){continue;}
			if (!is_save_Nev_index_ttbar){Nev_index_ttbar=Nev_index_counter; is_save_Nev_index_ttbar=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220812/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 			  		
	  	}

	  	// TTbar second 1000 sample
	  	for (int i=1000;i<=1999;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			if(i==1030 || i==1039 || i==1047 || i==1102){continue;}
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220812/0001/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		  			  		
	  	}

	  	// TTbar remain sample
	  	for (int i=2000;i<=2106;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_TT_TuneCUETP8M1_13TeV-powheg-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220812/0002/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		  			  		
	  	}


	  	// ST_s-channel
	  	bool is_save_Nev_index_ST_s_channel =false;
	  	for (int i=1;i<=31;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			if (!is_save_Nev_index_ST_s_channel){Nev_index_ST_s_channel=Nev_index_counter; is_save_Nev_index_ST_s_channel=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/SingleTop/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/crab_ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1MC25ns_eleIDjet_CMSSW7412_20151006/151007_220338/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 			  			  		
	  	}
	 
	 
	  	// ST_t-channel_antitop
	  	bool is_save_Nev_index_ST_t_channel_antitop =false;
	  	for (int i=1;i<=47;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			if (!is_save_Nev_index_ST_t_channel_antitop){Nev_index_ST_t_channel_antitop=Nev_index_counter; is_save_Nev_index_ST_t_channel_antitop=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/SingleTop/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1MC25ns_eleIDjet_CMSSW7412_20151006/151007_220428/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);	
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		  			  		
	  	}
	
	
	
	  	// ST_t-channel_top
	  	bool is_save_Nev_index_ST_t_channel_top=false;
	  	for (int i=1;i<=82;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
	  		if (i==22){continue;}
			if (!is_save_Nev_index_ST_t_channel_top){Nev_index_ST_t_channel_top=Nev_index_counter; is_save_Nev_index_ST_t_channel_top=true; }	
			string path ="/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/SingleTop/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1MC25ns_eleIDjet_CMSSW7412_20151006/151007_220524/0000/NCUGlobalTuples_"+int_to_string(i)+".root";	
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);	
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		  			  		
	  	}
	
	
	  	
	  	
	  	// ST_tW-antitop
	  	bool is_save_Nev_index_ST_tW_antitop=false;
	  	for (int i=1;i<=32;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			if (!is_save_Nev_index_ST_tW_antitop){Nev_index_ST_tW_antitop=Nev_index_counter; is_save_Nev_index_ST_tW_antitop=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/SingleTop/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1MC25ns_eleIDjet_CMSSW7412_20151006/151007_220611/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 			  			  		
	  	}
	  	
	  	
	  	// ST_tW-top
	  	bool is_save_Nev_index_ST_tW_top=false;
	  	for (int i=1;i<=32;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}	
			if (!is_save_Nev_index_ST_tW_top){Nev_index_ST_tW_top=Nev_index_counter; is_save_Nev_index_ST_tW_top=true; }	
			string path="/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/SingleTop/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/crab_ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1MC25ns_eleIDjet_CMSSW7412_20151006/151007_220712/0000/NCUGlobalTuples_"+int_to_string(i)+".root";	
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 	
   	
  		}
  		
	  	// DY inclusive
	  	bool is_save_Nev_index_DY_inclusive=false;
	  	for (int i=1;i<=671;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if(i>=300){continue;}

	  		if(i==124 || i==205){continue;}

			if (!is_save_Nev_index_DY_inclusive){Nev_index_DY_inclusive=Nev_index_counter; is_save_Nev_index_DY_inclusive=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_221020/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 			
	  	}

	  	// QCD 100-200
	  	bool is_save_Nev_index_QCD100to200=false;
	  	for (int i=1;i<=812;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if (i>20){continue;}
			
	  		if(i==9 ){continue;}
	  		if (!is_save_Nev_index_QCD100to200){Nev_index_QCD100to200=Nev_index_counter; is_save_Nev_index_QCD100to200=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220208/0000/NCUGlobalTuples_"+int_to_string(i)+".root";  				
			input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
	  	// QCD 200-300
	  	bool is_save_Nev_index_QCD200to300=false;
	  	for (int i=1;i<=515;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if (i>20){continue;}

	  		if( i==12 ){continue;}
	  		if (!is_save_Nev_index_QCD200to300){Nev_index_QCD200to300=Nev_index_counter; is_save_Nev_index_QCD200to300=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220428/0000/NCUGlobalTuples_"+int_to_string(i)+".root";  				
			input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
	  	// QCD 300-500
	  	bool is_save_Nev_index_QCD300to500=false;
	  	for (int i=1;i<=564;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if (i>20){continue;}

	  		if( i==9  ){continue;}
	  		if (!is_save_Nev_index_QCD300to500){Nev_index_QCD300to500=Nev_index_counter; is_save_Nev_index_QCD300to500=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/NCUGlobalTuples_"+int_to_string(i)+".root";  				
			input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
	  	// QCD 500-700
	  	bool is_save_Nev_index_QCD500to700=false;
	  	for (int i=1;i<=495;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if (i>20){continue;}

	  		//if(i==6 || i==7 || i==8 || i==9 || (i>=50 && i<=99) || i==370 ){continue;}
	  		if (!is_save_Nev_index_QCD500to700){Nev_index_QCD500to700=Nev_index_counter; is_save_Nev_index_QCD500to700=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT500to7000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/NCUGlobalTuples_"+int_to_string(i)+".root";  				
			input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}	
	  	
	  	// QCD 700-1000
	  	bool is_save_Nev_index_QCD700to1000=false;
	  	for (int i=1;i<=393;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if (i>20){continue;}

	  		if(i==6 || i==7 || i==8 || i==9 || (i>=50 && i<=99) || i==370 ){continue;}
	  		if (!is_save_Nev_index_QCD700to1000){Nev_index_QCD700to1000=Nev_index_counter; is_save_Nev_index_QCD700to1000=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/0000/NCUGlobalTuples_"+int_to_string(i)+".root";  				
			input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile); 
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}	
	  		
	  		  	
	  	// QCD 1000-1500
	  	bool is_save_Nev_index_QCD1000to1500=false;
	  	for (int i=1;i<=154;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if (i>20){continue;}

			//if(i>=400){continue;}
	  		if(i==41){continue;}

			if (!is_save_Nev_index_QCD1000to1500){Nev_index_QCD1000to1500=Nev_index_counter; is_save_Nev_index_QCD1000to1500=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220114/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);	  		
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
	  	
	  	
	  	
	  	// QCD 1500-2000
	  	bool is_save_Nev_index_QCD1500to2000=false;
	  	for (int i=1;i<=101;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if (i>20){continue;}
			//if(i>=400){continue;}
	  		if(i==1 || i==6 || i==84 || i==94 || i==95){continue;}

			if (!is_save_Nev_index_QCD1500to2000){Nev_index_QCD1500to2000=Nev_index_counter; is_save_Nev_index_QCD1500to2000=true; }
			string path="/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220249/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
	  	
	  	// QCD 2000-Inf
	  	bool is_save_Nev_index_QCD2000toInf=false;
	  	for (int i=1;i<=70;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if (i>20){continue;}
			
			if (!is_save_Nev_index_QCD2000toInf){Nev_index_QCD2000toInf=Nev_index_counter; is_save_Nev_index_QCD2000toInf=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220339/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}		  	
	  	
	  	// WJ(W+jets) 
	  	/*
	  	// WJ inclusive
	  	bool is_save_Nev_index_WJ=false;
	  	for (int i=1;i<=569;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}
			//if(i>=400){continue;}
	  		if(i==6 || i==93 || i==324 || i==514 || i==533 || i==567){continue;}
			if (!is_save_Nev_index_WJ){Nev_index_WJ=Nev_index_counter; is_save_Nev_index_WJ=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_221450/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
		*/
		
	  	// HT 100 - 200
	  	bool is_save_Nev_index_WJ100to200=false;
	  	for (int i=1;i<=235;i++)
	  	{
	  		if (i>=10){continue;}
			if (!is_save_Nev_index_WJ100to200){Nev_index_WJ100to200=Nev_index_counter; is_save_Nev_index_WJ100to200=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235712/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
	  	// HT 200 - 400
	  	bool is_save_Nev_index_WJ200to400=false;
	  	for (int i=1;i<=114;i++)
	  	{
	  		if (i>=10){continue;}
	  		if(i==7 || i==8 || i==9 || (i>=65 && i<=99) ){continue;}
			if (!is_save_Nev_index_WJ200to400){Nev_index_WJ200to400=Nev_index_counter; is_save_Nev_index_WJ200to400=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235758/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}	
	  	
	  	
	  	// HT 400 - 600
	  	bool is_save_Nev_index_WJ400to600=false;
	  	for (int i=1;i<=44;i++)
	  	{
	  		if (i>=10){continue;}
	  		if(i==38 ){continue;}
			if (!is_save_Nev_index_WJ400to600){Nev_index_WJ400to600=Nev_index_counter; is_save_Nev_index_WJ400to600=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235853/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}


	  	// HT 600 - 800
	  	bool is_save_Nev_index_WJ600to800=false;
	  	for (int i=1;i<=97;i++)
	  	{
	  		if (i>=10){continue;}
	  		
			if (!is_save_Nev_index_WJ600to800){Nev_index_WJ600to800=Nev_index_counter; is_save_Nev_index_WJ600to800=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151025_235938/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}

	  	// HT 800 - 1200
	  	bool is_save_Nev_index_WJ800to1200=false;
	  	for (int i=1;i<=39;i++)
	  	{
	  		if (i>=10){continue;}
			if (!is_save_Nev_index_WJ800to1200){Nev_index_WJ800to1200=Nev_index_counter; is_save_Nev_index_WJ800to1200=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151026_000033/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}

	  	// HT 1200 - 2500
	  	bool is_save_Nev_index_WJ1200to2500=false;
	  	for (int i=1;i<=7;i++)
	  	{
	  		if (i>=10){continue;}
			if (!is_save_Nev_index_WJ1200to2500){Nev_index_WJ1200to2500=Nev_index_counter; is_save_Nev_index_WJ1200to2500=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_ReMiniAODSIM/WJetsHTBinSampleReMiniAOD/crab_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_MC25ns_ReMiniAOD_20151026/151026_000152/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}		
		  	
	  	// WW 
	  	bool is_save_Nev_index_WW=false;
	  	for (int i=1;i<=44;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}

			if (!is_save_Nev_index_WW){Nev_index_WW=Nev_index_counter; is_save_Nev_index_WW=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_WW_TuneCUETP8M1_13TeV-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_220812/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	
	  	// WZ 
		bool is_save_Nev_index_WZ=false;
	  	for (int i=1;i<=45;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}

			if (!is_save_Nev_index_WZ){Nev_index_WZ=Nev_index_counter; is_save_Nev_index_WZ=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_WZ_TuneCUETP8M1_13TeV-pythia8MC25ns_eleIDjet_CMSSW7412_20151006/151007_221020/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}
	  	
	  	// ZZ 
	  	bool is_save_Nev_index_ZZ=false;
	  	for (int i=1;i<=36;i++)
	  	{
	  		if (i>N_applied_sample_max || i<N_applied_sample_min){continue;}

			if (!is_save_Nev_index_ZZ){Nev_index_ZZ=Nev_index_counter; is_save_Nev_index_ZZ=true; }
			string path = "/data7/khurana/NCUGlobalTuples/SPRING15_2015_10_12/crab_ZZ_TuneCUETP8M1_13TeV-pythia8MC25ns_eleIDjet_CSW7412_20151006/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		input_files.push_back(path);  		
	  		
	  		// to count Nev of this sample
	  		vector<string> this_rootfile;
	  		this_rootfile.push_back(path);
	  		TreeReader this_reader(this_rootfile);
			
			Nev_index_counter+=this_reader.GetEntriesFast(); 		
	  	}	  	

	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_100.root");
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_200.root");
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_300.root");
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_400.root");
	  	
	  				
	  	//TreeReader inner_data_for_MC(MC);
	  	TreeReader inner_data_for_MC(input_files);
	  	Nev_MC = inner_data_for_MC.GetEntriesFast();

		// ===============================================calculate # of event each sample (considered event weight)=============================================
		for (Long64_t Nth = 0; Nth < Nev_MC; Nth++)
	  	{	  		
	  		inner_data_for_MC.GetEntry(Nth);	  		
	    		Float_t mcWeight=inner_data_for_MC.GetFloat("mcWeight");
	    		int event_weight=(int)get_event_weight(Nth,mcWeight);
  		
			if (Nth>=Nev_index_ZZ) {Nev_ZZ++;}
			else if(Nth>=Nev_index_WZ) {Nev_WZ++;}
			else if(Nth>=Nev_index_WW) {Nev_WW++; }
			/*else if (Nth>=Nev_index_WJ) 
			{
				if (mcWeight>0){Nev_WJ+=1;}
				else {Nev_WJ-=1;}
			}*/
			else if (Nth>=Nev_index_WJ1200to2500) {Nev_WJ1200to2500++;}
			else if(Nth>=Nev_index_WJ800to1200) {Nev_WJ800to1200++;}
			else if(Nth>=Nev_index_WJ600to800) {Nev_WJ600to800++;}
			else if(Nth>=Nev_index_WJ400to600) {Nev_WJ400to600++;}
			else if(Nth>=Nev_index_WJ200to400) {Nev_WJ200to400++;}
			else if(Nth>=Nev_index_WJ100to200) {Nev_WJ100to200++;}
			else if(Nth>=Nev_index_QCD2000toInf) {Nev_QCD2000toInf++; }
			else if(Nth>=Nev_index_QCD1500to2000) {Nev_QCD1500to2000++; }
			else if(Nth>=Nev_index_QCD1000to1500) {Nev_QCD1000to1500++; }
			else if(Nth>=Nev_index_QCD700to1000) {Nev_QCD700to1000++; }
			else if(Nth>=Nev_index_QCD500to700) {Nev_QCD500to700++; }
			else if(Nth>=Nev_index_QCD300to500) {Nev_QCD300to500++; }
			else if(Nth>=Nev_index_QCD200to300) {Nev_QCD200to300++; }
			else if(Nth>=Nev_index_QCD100to200) {Nev_QCD100to200++; }
			else if(Nth>=Nev_index_DY_inclusive) 
			{
				if (mcWeight>0){Nev_DY_inclusive+=1;}
				else {Nev_DY_inclusive-=1;}
			}
			else if(Nth>=Nev_index_ST_tW_top) {Nev_ST_tW_top++; }
			else if(Nth>=Nev_index_ST_tW_antitop) {Nev_ST_tW_antitop++; }
			else if(Nth>=Nev_index_ST_t_channel_top) {Nev_ST_t_channel_top++; }
			else if(Nth>=Nev_index_ST_t_channel_antitop) {Nev_ST_t_channel_antitop++; }
			else if(Nth>=Nev_index_ST_s_channel)
			{
				if (mcWeight>0){Nev_ST_s_channel+=1;}
				else {Nev_ST_s_channel-=1;}
			}
			else if(Nth>=Nev_index_ttbar) {Nev_ttbar++; }
	
	
	  	}
	  	
	  	
	  	
		cout<<"============================================ number of event(each sample)  ============================================"<<endl;
		cout<<"============================================ note : before any scaling  ============================================"<<endl;
	  	cout<< "Nev_MC(all MC) : "<< Nev_MC<< endl;

	  	cout<< "Nev_ttbar : "<< Nev_ttbar<< endl;  	
	  	cout<< "Nev_ST_s_channel : "<< Nev_ST_s_channel<< endl;
	  	cout<< "Nev_ST_t_channel_antitop : "<< Nev_ST_t_channel_antitop<< endl;
	  	cout<< "Nev_ST_t_channel_top : "<< Nev_ST_t_channel_top<< endl;
	  	cout<< "Nev_ST_tW_antitop : "<< Nev_ST_tW_antitop<< endl;
	  	cout<< "Nev_ST_tW_top : "<< Nev_ST_tW_top<< endl;
	  	cout<< "Nev_DY_inclusive : "<< Nev_DY_inclusive<< endl;
	  	cout<< "Nev_QCD100to200 : "<< Nev_QCD100to200<< endl;
	  	cout<< "Nev_QCD200to300 : "<< Nev_QCD200to300<< endl;
	  	cout<< "Nev_QCD300to500 : "<< Nev_QCD300to500<< endl;
	  	cout<< "Nev_QCD500to700 : "<< Nev_QCD500to700<< endl;
	  	cout<< "Nev_QCD700to1000 : "<< Nev_QCD700to1000<< endl;
	  	cout<< "Nev_QCD1000to1500 : "<< Nev_QCD1000to1500<< endl;
	  	cout<< "Nev_QCD1500to2000 : "<< Nev_QCD1500to2000<< endl;
	  	cout<< "Nev_QCD2000toInf : "<< Nev_QCD2000toInf<< endl;
	  	//cout<< "Nev_WJ : "<< Nev_WJ << endl;
	  	cout<< "Nev_WJ100to200 : "<< Nev_WJ100to200 << endl;
	  	cout<< "Nev_WJ200to400 : "<< Nev_WJ200to400 << endl;
	  	cout<< "Nev_WJ400to600 : "<< Nev_WJ400to600 << endl;
	  	cout<< "Nev_WJ600to800 : "<< Nev_WJ600to800 << endl;
	  	cout<< "Nev_WJ800to1200 : "<< Nev_WJ800to1200 << endl;
	  	cout<< "Nev_WJ1200to2500 : "<< Nev_WJ1200to2500 << endl;
	  	cout<< "Nev_WW : "<< Nev_WW<< endl;
	  	cout<< "Nev_WZ : "<< Nev_WZ<< endl;
	  	cout<< "Nev_ZZ : "<< Nev_ZZ<< endl;


		// calculate scale factor 
		sf_ttbar = data_lumi_13_Tev/((double)Nev_ttbar/Xs_ttbar);
		sf_ST_s_channel=data_lumi_13_Tev/((double)Nev_ST_s_channel/Xs_ST_s_channel);
		sf_ST_t_channel_antitop=data_lumi_13_Tev/((double)Nev_ST_t_channel_antitop/Xs_ST_t_channel_antitop);
		sf_ST_t_channel_top=data_lumi_13_Tev/((double)Nev_ST_t_channel_top/Xs_ST_t_channel_top);
		sf_ST_tW_antitop=data_lumi_13_Tev/((double)Nev_ST_tW_antitop/Xs_ST_tW_antitop);
		sf_ST_tW_top=data_lumi_13_Tev/((double)Nev_ST_tW_top/Xs_ST_tW_top);
		sf_DY_inclusive=data_lumi_13_Tev/((double)Nev_DY_inclusive/Xs_DY_inclusive);		
		sf_QCD100to200=data_lumi_13_Tev/((double)Nev_QCD100to200/Xs_QCD100to200);
		sf_QCD200to300=data_lumi_13_Tev/((double)Nev_QCD200to300/Xs_QCD200to300);
		sf_QCD300to500=data_lumi_13_Tev/((double)Nev_QCD300to500/Xs_QCD300to500);
		sf_QCD500to700=data_lumi_13_Tev/((double)Nev_QCD500to700/Xs_QCD500to700);	
		sf_QCD700to1000=data_lumi_13_Tev/((double)Nev_QCD700to1000/Xs_QCD700to1000);
		sf_QCD1000to1500=data_lumi_13_Tev/((double)Nev_QCD1000to1500/Xs_QCD1000to1500);
		sf_QCD1500to2000=data_lumi_13_Tev/((double)Nev_QCD1500to2000/Xs_QCD1500to2000);
		sf_QCD2000toInf=data_lumi_13_Tev/((double)Nev_QCD2000toInf/Xs_QCD2000toInf);
		//sf_WJ=data_lumi_13_Tev/((double)Nev_WJ/Xs_WJ);
		sf_WJ100to200=data_lumi_13_Tev/((double)Nev_WJ100to200/Xs_WJ100to200);
		sf_WJ200to400=data_lumi_13_Tev/((double)Nev_WJ200to400/Xs_WJ200to400);
		sf_WJ400to600=data_lumi_13_Tev/((double)Nev_WJ400to600/Xs_WJ400to600);
		sf_WJ600to800=data_lumi_13_Tev/((double)Nev_WJ600to800/Xs_WJ600to800);
		sf_WJ800to1200=data_lumi_13_Tev/((double)Nev_WJ800to1200/Xs_WJ800to1200);
		sf_WJ1200to2500=data_lumi_13_Tev/((double)Nev_WJ1200to2500/Xs_WJ1200to2500);		
		sf_WW=data_lumi_13_Tev/((double)Nev_WW/Xs_WW);
		sf_WZ=data_lumi_13_Tev/((double)Nev_WZ/Xs_WZ);
		sf_ZZ =data_lumi_13_Tev/((double)Nev_ZZ/Xs_ZZ); 	
	  		
	  	
	  	cout<<"============================================ scale factor(each sample) ============================================"<<endl<<endl;

	  	cout<< "sf_ttbar : "<< sf_ttbar<< endl;  	
	  	cout<< "sf_ST_s_channel : "<< sf_ST_s_channel<< endl;
	  	cout<< "sf_ST_t_channel_antitop : "<< sf_ST_t_channel_antitop<< endl;
	  	cout<< "sf_ST_t_channel_top : "<< sf_ST_t_channel_top<< endl;
	  	cout<< "sf_ST_tW_antitop : "<< sf_ST_tW_antitop<< endl;
	  	cout<< "sf_ST_tW_top : "<< sf_ST_tW_top<< endl;
	  	cout<< "sf_DY_inclusive : "<< sf_DY_inclusive<< endl;
	  	cout<< "sf_QCD100to200 : "<< sf_QCD100to200<< endl;
	  	cout<< "sf_QCD200to300 : "<< sf_QCD200to300<< endl;
	  	cout<< "sf_QCD300to500 : "<< sf_QCD300to500<< endl;
	  	cout<< "sf_QCD500to700 : "<< sf_QCD500to700<< endl;
	  	cout<< "sf_QCD700to1000 : "<< sf_QCD700to1000<< endl;
	  	cout<< "sf_QCD1000to1500 : "<< sf_QCD1000to1500<< endl;
	  	cout<< "sf_QCD1500to2000 : "<< sf_QCD1500to2000<< endl;
	  	cout<< "sf_QCD2000toInf : "<< sf_QCD2000toInf<< endl;
	        //cout<< "sf_WJ : "<< sf_WJ<< endl;
	  	cout<< "sf_WJ100to200 : "<< sf_WJ100to200<< endl;
	  	cout<< "sf_WJ200to400 : "<< sf_WJ200to400<< endl;
	  	cout<< "sf_WJ400to600 : "<< sf_WJ400to600<< endl;
	  	cout<< "sf_WJ600to800 : "<< sf_WJ600to800<< endl;
	  	cout<< "sf_WJ800to1200 : "<< sf_WJ800to1200<< endl;
	  	cout<< "sf_WJ1200to2500 : "<< sf_WJ1200to2500<< endl;
	  	cout<< "sf_WW : "<< sf_WW<< endl;
	  	cout<< "sf_WZ : "<< sf_WZ<< endl;
	  	cout<< "sf_ZZ : "<< sf_ZZ<< endl;
	  	
	  	cout<<"============================================ number of event(each sample) after scaleing ============================================"<<endl<<endl;
	  	cout<< "Nev_ttbar_scaled : "<< Nev_ttbar*sf_ttbar << endl;  	
	  	cout<< "Nev_ST_s_channel_scaled : "<< Nev_ST_s_channel*sf_ST_s_channel << endl;
	  	cout<< "Nev_ST_t_channel_antitop_scaled : "<< Nev_ST_t_channel_antitop*sf_ST_t_channel_antitop << endl;
	  	cout<< "Nev_ST_t_channel_top_scaled : "<< Nev_ST_t_channel_top*sf_ST_t_channel_top << endl;
	  	cout<< "Nev_ST_tW_antitop_scaled : "<< Nev_ST_tW_antitop*sf_ST_tW_antitop << endl;
	  	cout<< "Nev_ST_tW_top_scaled : "<< Nev_ST_tW_top*sf_ST_tW_top << endl;
	  	cout<< "Nev_DY_inclusive_scaled : "<< Nev_DY_inclusive*sf_DY_inclusive << endl;
	  	cout<< "Nev_QCD100to200_scaled : "<< Nev_QCD100to200*sf_QCD100to200 << endl;
	  	cout<< "Nev_QCD200to300_scaled : "<< Nev_QCD200to300*sf_QCD200to300 << endl;
	  	cout<< "Nev_QCD300to500_scaled : "<< Nev_QCD300to500*sf_QCD300to500 << endl;
	  	cout<< "Nev_QCD500to700_scaled : "<< Nev_QCD500to700*sf_QCD500to700 << endl;
	  	cout<< "Nev_QCD700to1000_scaled : "<< Nev_QCD700to1000*sf_QCD700to1000 << endl;
	  	cout<< "Nev_QCD1000to1500_scaled : "<< Nev_QCD1000to1500*sf_QCD1000to1500 << endl;
	  	cout<< "Nev_QCD1500to2000_scaled : "<< Nev_QCD1500to2000*sf_QCD1500to2000 << endl;
	  	cout<< "Nev_QCD2000toInf_scaled : "<< Nev_QCD2000toInf*sf_QCD2000toInf << endl;
	  	//cout<< "Nev_WJ_scaled : "<< Nev_WJ*sf_WJ << endl;
	  	cout<< "Nev_WJ100to200_scaled : "<< Nev_WJ100to200*sf_WJ100to200 << endl;
	  	cout<< "Nev_WJ200to400_scaled : "<< Nev_WJ200to400*sf_WJ200to400 << endl;
	  	cout<< "Nev_WJ400to600_scaled : "<< Nev_WJ400to600*sf_WJ400to600 << endl;
	  	cout<< "Nev_WJ600to800_scaled : "<< Nev_WJ600to800*sf_WJ600to800 << endl;
	  	cout<< "Nev_WJ800to1200_scaled : "<< Nev_WJ800to1200*sf_WJ800to1200 << endl;
	  	cout<< "Nev_WJ1200to2500_scaled : "<< Nev_WJ1200to2500*sf_WJ1200to2500 << endl;
	  	cout<< "Nev_WW_scaled : "<< Nev_WW*sf_WW << endl;
	  	cout<< "Nev_WZ_scaled : "<< Nev_WZ*sf_WZ << endl;
	  	cout<< "Nev_ZZ_scaled : "<< Nev_ZZ*sf_ZZ << endl;
	  	cout<<endl;
	  	cout<< "Nev_MC_scaled : "<< 
	  	Nev_ttbar*sf_ttbar+
	  	Nev_ST_s_channel*sf_ST_s_channel+
	  	Nev_ST_t_channel_antitop*sf_ST_t_channel_antitop+
	  	Nev_ST_t_channel_top*sf_ST_t_channel_top+
	  	Nev_ST_tW_antitop*sf_ST_tW_antitop+
	  	Nev_ST_tW_top*sf_ST_tW_top+
	  	Nev_DY_inclusive*sf_DY_inclusive+
	  	Nev_QCD100to200*sf_QCD100to200+
	  	Nev_QCD200to300*sf_QCD200to300+
	  	Nev_QCD300to500*sf_QCD300to500+
	  	Nev_QCD500to700*sf_QCD500to700+
	  	Nev_QCD700to1000*sf_QCD700to1000+
	  	Nev_QCD1000to1500*sf_QCD1000to1500+
	  	Nev_QCD1500to2000*sf_QCD1500to2000+
	  	Nev_QCD2000toInf*sf_QCD2000toInf+
	  	//Nev_WJ*sf_WJ+
	  	Nev_WJ100to200*sf_WJ100to200+
	  	Nev_WJ200to400*sf_WJ200to400+
	  	Nev_WJ400to600*sf_WJ400to600+
	  	Nev_WJ600to800*sf_WJ600to800+
	  	Nev_WJ800to1200*sf_WJ800to1200+
	  	Nev_WJ1200to2500*sf_WJ1200to2500+
	  	Nev_WW*sf_WW+
	  	Nev_WZ*sf_WZ +
	  	Nev_ZZ*sf_ZZ<< endl;

	  		  	
  	
	  	
		//configure
		vector<event> events_MC;
	
		int n_THINj=0;
	
		//matching

		// event loop : set up MC
	  	//for (Long64_t Nth = 0; Nth < Nev_MC; Nth++)
	  	for (Long64_t Nth = 0; Nth < Nev_MC; Nth++)
	    	{
			event this_ev;

	 		inner_data_for_MC.GetEntry(Nth);

			double scaled_value=get_scaled_value(Nth);
			// amcatNlo sample should apply mcWeight
			Float_t mcWeight=inner_data_for_MC.GetFloat("mcWeight");
			
			this_ev.event_weight = get_event_weight(Nth, mcWeight);		

			// trigger var
			string* hlt_trigName = inner_data_for_MC.GetPtrString("hlt_trigName");
			vector<bool> &hlt_trigResult = *((vector<bool>*) inner_data_for_MC.GetPtr("hlt_trigResult"));	
			int string_size = inner_data_for_MC.GetPtrStringSize();

			// gen_par var
	      		Int_t nGenPar=inner_data_for_MC.GetInt("nGenPar");
	      		/*
	      		TClonesArray* genParP4_ = (TClonesArray*)data.GetPtrTObject("genParP4");
	      		TLorentzVector* genParP4=(TLorentzVector*)genParP4_;	   */
	      		TClonesArray* genParP4 = (TClonesArray*) inner_data_for_MC.GetPtrTObject("genParP4"); 		
	     			
			Int_t *genParId=inner_data_for_MC.GetPtrInt("genParId");	
			Int_t *genParIndex=inner_data_for_MC.GetPtrInt("genParIndex");
			Int_t *genMomParId=inner_data_for_MC.GetPtrInt("genMomParId");
			Int_t *genMo1=inner_data_for_MC.GetPtrInt("genMo1");  //not mom id
			Int_t *genMo2=inner_data_for_MC.GetPtrInt("genMo2");  //not mom id
			Int_t *genParSt=inner_data_for_MC.GetPtrInt("genParSt"); 



	      		// THINjet var
	      		Int_t THINnJet=inner_data_for_MC.GetInt("THINnJet");
	      		
	      		TClonesArray* THINjetP4 = (TClonesArray*) inner_data_for_MC.GetPtrTObject("THINjetP4");
	      		Float_t *THINjetCSV=inner_data_for_MC.GetPtrFloat("THINjetCSV");
	      	      //	Float_t* THINjetPrunedEn = inner_data_for_MC.GetPtrFloat("THINjetPrunedEn");
	      	      //	Int_t* THINjetPassID = inner_data_for_MC.GetPtrInt("THINjetPassID");
	      	      //	Float_t* THINjetTau1 = inner_data_for_MC.GetPtrFloat("THINjetTau1");
			//	Float_t* THINjetTau2 = inner_data_for_MC.GetPtrFloat("THINjetTau2");
			Float_t* THINjetCISVV2 = inner_data_for_MC.GetPtrFloat("THINjetCISVV2");
			Int_t * THINjetHadronFlavor=inner_data_for_MC.GetPtrInt("THINjetHadronFlavor"); 
			
			
			vector<bool> &THINjetPassIDTight = *((vector<bool>*) inner_data_for_MC.GetPtr("THINjetPassIDTight"));			
			vector<bool> &THINjetPassIDLoose  = *((vector<bool>*) inner_data_for_MC.GetPtr("THINjetPassIDLoose"));

			Float_t *THINjetCEmEF=inner_data_for_MC.GetPtrFloat("THINjetCEmEF");
			Float_t *THINjetCHadEF=inner_data_for_MC.GetPtrFloat("THINjetCHadEF");	
			Float_t *THINjetNEmEF=inner_data_for_MC.GetPtrFloat("THINjetNEmEF");	
			Float_t *THINjetNHadEF=inner_data_for_MC.GetPtrFloat("THINjetNHadEF");
			Float_t *THINjetPhoEF=inner_data_for_MC.GetPtrFloat("THINjetPhoEF");

			// muon var
	      		Int_t nMu =inner_data_for_MC.GetInt("nMu");
	      		TClonesArray* muP4 = (TClonesArray*) inner_data_for_MC.GetPtrTObject("muP4");
	      		
	    		
	    		vector<bool> &isGlobalMuon  = *((vector<bool>*) inner_data_for_MC.GetPtr("isGlobalMuon"));
	    		vector<bool> &isTrackerMuon  = *((vector<bool>*) inner_data_for_MC.GetPtr("isTrackerMuon"));
	    		vector<bool> &isTightMuon  = *((vector<bool>*) inner_data_for_MC.GetPtr("isTightMuon"));
	    		vector<bool> &isMediumMuon  = *((vector<bool>*) inner_data_for_MC.GetPtr("isMediumMuon"));
	    		vector<bool> &isLooseMuon = *((vector<bool>*) inner_data_for_MC.GetPtr("isLooseMuon"));
	    		vector<bool> &isHighPtMuon  = *((vector<bool>*) inner_data_for_MC.GetPtr("isHighPtMuon"));
	    		vector<bool> &isSoftMuon  = *((vector<bool>*) inner_data_for_MC.GetPtr("isSoftMuon"));
	    		
	    		Float_t* muChHadIso = inner_data_for_MC.GetPtrFloat("muChHadIso");
	    		Float_t* muNeHadIso = inner_data_for_MC.GetPtrFloat("muNeHadIso");
	    		Float_t* muGamIso = inner_data_for_MC.GetPtrFloat("muGamIso");
	    		Float_t* muPUPt = inner_data_for_MC.GetPtrFloat("muPUPt");
	    		
	    		
/*
	    		Int_t*   muTrkLayers  = inner_data_for_MC.GetPtrInt("muTrkLayers");
	    		Int_t*   muPixelHits  = inner_data_for_MC.GetPtrInt("muPixelHits");
	    		Int_t*   muHits       = inner_data_for_MC.GetPtrInt("muHits");
	    		Int_t*   muMatches    = inner_data_for_MC.GetPtrInt("muMatches");
	    		Float_t* mudxy        = inner_data_for_MC.GetPtrFloat("mudxy");
	    		Float_t* mudz         = inner_data_for_MC.GetPtrFloat("mudz");

	   		//Float_t *muCorrPfIso=inner_data_for_MC.GetPtrFloat("muCorrPfIso");				
			//Int_t* muPassID = inner_data_for_MC.GetPtrInt("muPassID");
			//Float_t* muCorrTrkIso = inner_data_for_MC.GetPtrFloat("muCorrTrkIso");
*/
	      		// ele var
	      		Int_t nEle =inner_data_for_MC.GetInt("nEle");
	      		TClonesArray* eleP4 = (TClonesArray*) inner_data_for_MC.GetPtrTObject("eleP4");
	      		
	      		vector<bool> &eleIsPassTight  = *((vector<bool>*) inner_data_for_MC.GetPtr("eleIsPassTight"));
	      		vector<bool> &eleIsPassMedium  = *((vector<bool>*) inner_data_for_MC.GetPtr("eleIsPassMedium"));
	      		vector<bool> &eleIsPassLoose  = *((vector<bool>*) inner_data_for_MC.GetPtr("eleIsPassLoose"));
	      		
	      		
/*
			//Float_t *eleDelEtaIn=inner_data_for_MC.GetPtrFloat("eleDelEtaIn");
			//Float_t *eleDelPhiIn=inner_data_for_MC.GetPtrFloat("eleDelPhiIn");
			//Float_t *eleSigIhIh=inner_data_for_MC.GetPtrFloat("eleSigIhIh");
			Float_t *eleHoverE=inner_data_for_MC.GetPtrFloat("eleHoverE");
			//Float_t *eleDxy=inner_data_for_MC.GetPtrFloat("eleDxy");
			Float_t *eleDz=inner_data_for_MC.GetPtrFloat("eleDz");
			Float_t *eleEoverP=inner_data_for_MC.GetPtrFloat("eleEoverP");
			//Float_t *eleCorrPfIso=inner_data_for_MC.GetPtrFloat("eleCorrPfIso");
			//Int_t *elePassConv=inner_data_for_MC.GetPtrInt("elePassConv");
			//Float_t *eleMissingHits=inner_data_for_MC.GetPtrFloat("eleMissingHits");

			//Int_t *elePassID=inner_data_for_MC.GetPtrInt("elePassID");
			//Float_t *eleUserTrkIso=inner_data_for_MC.GetPtrFloat("eleUserTrkIso");
			//Float_t *eleUserCalIso=inner_data_for_MC.GetPtrFloat("eleUserCalIso");
			Float_t eleRho=inner_data_for_MC.GetFloat("eleRho");
*/

			//met
			Float_t pfMetCorrPt=inner_data_for_MC.GetFloat("pfMetCorrPt");
			Float_t pfMetCorrPhi=inner_data_for_MC.GetFloat("pfMetCorrPhi");

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
			/*
	      			ele ele_i(elePt[i],eleEta[i],elePhi[i],eleEnergy[i],eleM[i],
					 eleDelEtaIn[i],eleDelPhiIn[i],eleSigIhIh[i],eleHoE[i],eleDxy[i],
					 eleDz[i],eleEoverP[i],eleCorrPfIso[i],elePassConv[i],eleMissingHits[i],
					 elePassID[i],eleUserTrkIso[i],eleUserCalIso[i],eleRho); */
	      			
	      		
				ele ele_i( *((TLorentzVector*)eleP4->At(i)),eleIsPassTight[i],eleIsPassMedium[i],eleIsPassLoose[i]);

				if (ele_i.p4.Pt()==0){cout<<"e";}

				this_ev.eles.push_back(ele_i);
			}

	      		for (Int_t i=0;i<nMu;i++)
			{
			/*
	      			mu mu_i(muPt[i],muEta[i],muPhi[i],0,muM[i], // no muEnergy
			 		isGlobalMuon[i],isTrackerMuon[i],muTrkLayers[i],muPixelHits[i],muHits[i],muMatches[i],mudxy[i],mudz[i],        				
					nuCorrPfIso[i],muPassID[i],muCorrTrkIso[i]);
			*/
	    		
	    			mu mu_i( *((TLorentzVector*)muP4->At(i)), isGlobalMuon[i], isTrackerMuon[i], isTightMuon[i], isMediumMuon[i], isLooseMuon[i],isHighPtMuon[i], isSoftMuon[i], muChHadIso[i], muNeHadIso[i], muGamIso[i], muPUPt[i]);
				if (mu_i.p4.Pt()==0){cout<<"m";}
				this_ev.mus.push_back(mu_i);
			}
	      	      				
	      		for (Int_t i=0;i<THINnJet;i++)
			{	
				//if (THINjetCSV[i]<0 || THINjetCSV[i]>1){continue;}
				if (THINjetCISVV2[i]<0 || THINjetCISVV2[i]>1){continue;}
		
				THINjet jet_i( *((TLorentzVector*)THINjetP4->At(i)), THINjetCSV[i], THINjetCISVV2[i], THINjetHadronFlavor[i], THINjetPassIDTight[i],THINjetPassIDLoose[i], THINjetCEmEF[i], THINjetCHadEF[i], THINjetNEmEF[i], THINjetNHadEF[i], THINjetPhoEF[i]);
				if (jet_i.p4.Pt()==0){cout<<"j";}
				this_ev.THINjets.push_back(jet_i);

			}
			for(int i=0;i<nGenPar;i++)  //gen par loop for TTbar sample
			{
				// only include final state gen particles
				if 
				( 
					( abs(genParId[i])<=6  && abs(genParId[i])>=1  )||
					( abs(genParId[i])==11 )|| 
					( abs(genParId[i])==13 )
					
				)
				{
					/*if ( (abs(genParId[i])>=1 && abs(genParId[i])<=8) || 
					(abs(genParId[i])>=11 && abs(genParId[i])<=18) || 
					(abs(genParId[i])>=21 && abs(genParId[i])<=22)){}
					else{continue;} */

					if(    (*((TLorentzVector*)genParP4->At(i))).Pt()==0   ){continue;}
					gen_par this_g = gen_par( *((TLorentzVector*)genParP4->At(i)),
						genParId[i], genParIndex[i], genMomParId[i], genMo1[i], genMo2[i], genParSt[i]);
					this_ev.gen_pars.push_back(this_g);
				}
				
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
			events_MC.push_back(this_ev);
			
		
			// draw insignificant var (these insignificant var did not saved in class event, excluding pt & eta )
			// when event pass event selection
			// // ele
			
			Float_t *eledEtaAtVtx=inner_data_for_MC.GetPtrFloat("eledEtaAtVtx");
			Float_t *eledPhiAtVtx=inner_data_for_MC.GetPtrFloat("eledPhiAtVtx");
			Float_t *eleD0=inner_data_for_MC.GetPtrFloat("eleD0");
			Float_t *eleEtaseedAtVtx=inner_data_for_MC.GetPtrFloat("eleEtaseedAtVtx");
			Float_t *eleHoverE=inner_data_for_MC.GetPtrFloat("eleHoverE");
			Int_t *eleMissHits=inner_data_for_MC.GetPtrInt("eleMissHits");
			Float_t *eleMiniIso=inner_data_for_MC.GetPtrFloat("eleMiniIso");
			Float_t *eleSigmaIEtaIEta=inner_data_for_MC.GetPtrFloat("eleSigmaIEtaIEta");
	
			if (is_pass_ev_selection(this_ev))
			{
				for (int i=0;i<nEle;i++)
				{
					sh_ele_pt_a_evcut->fill_hist( (*((TLorentzVector*)eleP4->At(i))).Pt(), scaled_value*events_MC[Nth].event_weight, Nth);
			 		sh_ele_eta_a_evcut->fill_hist( (*((TLorentzVector*)eleP4->At(i))).Eta(), scaled_value*events_MC[Nth].event_weight, Nth); 
			 		sh_ele_dEtaAtVtx_a_evcut->fill_hist(eledEtaAtVtx[i], scaled_value*events_MC[Nth].event_weight, Nth); 
			 		sh_ele_dPhiAtVtx_a_evcut->fill_hist(eledPhiAtVtx[i], scaled_value*events_MC[Nth].event_weight, Nth); 
		
					sh_ele_D0_a_evcut->fill_hist(eleD0[i], scaled_value*events_MC[Nth].event_weight, Nth);  
			 		sh_ele_EtaseedAtVtx_a_evcut ->fill_hist(eleEtaseedAtVtx[i], scaled_value*events_MC[Nth].event_weight, Nth); 
			 		sh_ele_HoverE_a_evcut->fill_hist(eleHoverE[i], scaled_value*events_MC[Nth].event_weight, Nth); 
			 		sh_ele_MissHits_a_evcut->fill_hist(eleMissHits[i], scaled_value*events_MC[Nth].event_weight, Nth); 
			 		sh_ele_MiniIso_a_evcut->fill_hist(eleMiniIso[i], scaled_value*events_MC[Nth].event_weight, Nth); 
					sh_ele_SigmaIEtaIEta_a_evcut->fill_hist(eleSigmaIEtaIEta[i], scaled_value*events_MC[Nth].event_weight, Nth); 
				}	
			}

	
			// // mu (only 1 mu decayed from ttbar)
			
			Float_t *mudxy=inner_data_for_MC.GetPtrFloat("mudxy");	
			Float_t *mudz=inner_data_for_MC.GetPtrFloat("mudz");		
			Int_t *muHits=inner_data_for_MC.GetPtrInt("muHits");		
			Int_t *muMatches=inner_data_for_MC.GetPtrInt("muMatches");		
			Float_t *muMiniIso=inner_data_for_MC.GetPtrFloat("muMiniIso");		
			Int_t *muPixelHits=inner_data_for_MC.GetPtrInt("muPixelHits");		
			Int_t *muTrkLayers=inner_data_for_MC.GetPtrInt("muTrkLayers");		
			Float_t *muTrkPt=inner_data_for_MC.GetPtrFloat("muTrkPt");

			if (is_pass_ev_selection(this_ev))
			{
				bool find_candidate_mu=false;
				for (int i=0;i<nMu;i++)
				{
					if (isTightMuon[i] && (*((TLorentzVector*)muP4->At(i))).Pt()>50){find_candidate_mu=true;}
					else{continue;}
					
					sh_mu_pt_a_evcut->fill_hist( (*((TLorentzVector*)muP4->At(i))).Pt(), scaled_value*events_MC[Nth].event_weight, Nth);	
					sh_mu_eta_a_evcut->fill_hist( (*((TLorentzVector*)muP4->At(i))).Eta(), scaled_value*events_MC[Nth].event_weight, Nth);
					sh_mu_dxy_a_evcut->fill_hist(mudxy[i], scaled_value*events_MC[Nth].event_weight, Nth); 
					sh_mu_dz_a_evcut->fill_hist(mudz[i], scaled_value*events_MC[Nth].event_weight, Nth); 	
					sh_mu_Hits_a_evcut->fill_hist(muHits[i], scaled_value*events_MC[Nth].event_weight, Nth); 
					sh_mu_Matches_a_evcut->fill_hist(muMatches[i], scaled_value*events_MC[Nth].event_weight, Nth); 
					sh_mu_MiniIso_a_evcut->fill_hist(muMiniIso[i], scaled_value*events_MC[Nth].event_weight, Nth); 
					sh_mu_PixelHits_a_evcut->fill_hist(muPixelHits[i], scaled_value*events_MC[Nth].event_weight, Nth); 
					sh_mu_TrkLayers_a_evcut->fill_hist(muTrkLayers[i], scaled_value*events_MC[Nth].event_weight, Nth); 
					sh_mu_TrkPt_a_evcut->fill_hist(muTrkPt[i], scaled_value*events_MC[Nth].event_weight, Nth); 
					
					if (find_candidate_mu){break;}
				}	
			}

			// THINjet var (only 4 good jet decayed from ttbar)
			
			if (is_pass_ev_selection(this_ev))
			{
				int N_candidate_jet=0;
				for(int i=0;i<this_ev.THINjets_deoverlap_with_lepton.size();i++)
				{
					THINjet this_j = this_ev.THINjets_deoverlap_with_lepton[i];
					if (!is_good_jet(this_j)){continue;}
					else {N_candidate_jet++;}					
					
					sh_THINjet_CEmEF_a_evcut->fill_hist(this_j.CEmEF,scaled_value*events_MC[Nth].event_weight, Nth);
					sh_THINjet_CHadEF_a_evcut->fill_hist(this_j.CHadEF,scaled_value*events_MC[Nth].event_weight, Nth);
					sh_THINjet_NEmEF_a_evcut->fill_hist(this_j.NEmEF,scaled_value*events_MC[Nth].event_weight, Nth);
					sh_THINjet_NHadEF_a_evcut->fill_hist(this_j.NHadEF,scaled_value*events_MC[Nth].event_weight, Nth);
					sh_THINjet_PhoEF_a_evcut->fill_hist(this_j.PhoEF,scaled_value*events_MC[Nth].event_weight, Nth);
								
					sh_THINjet_cisvv2_a_evcut->fill_hist(this_j.cisvv2, scaled_value*events_MC[Nth].event_weight, Nth);
					sh_THINjet_pt_a_evcut->fill_hist(this_j.p4.Pt(), scaled_value*events_MC[Nth].event_weight, Nth);
					sh_THINjet_eta_a_evcut->fill_hist(this_j.p4.Eta(), scaled_value*events_MC[Nth].event_weight, Nth);
				
					if (btagging(this_j))
					{				
						sh_THINBjet_CEmEF_a_evcut->fill_hist(this_j.CEmEF,scaled_value*events_MC[Nth].event_weight, Nth);
						sh_THINBjet_CHadEF_a_evcut->fill_hist(this_j.CHadEF,scaled_value*events_MC[Nth].event_weight, Nth);
						sh_THINBjet_NEmEF_a_evcut->fill_hist(this_j.NEmEF,scaled_value*events_MC[Nth].event_weight, Nth);
						sh_THINBjet_NHadEF_a_evcut->fill_hist(this_j.NHadEF,scaled_value*events_MC[Nth].event_weight, Nth);
						sh_THINBjet_PhoEF_a_evcut->fill_hist(this_j.PhoEF,scaled_value*events_MC[Nth].event_weight, Nth);
					
						sh_THINBtagjet_cisvv2_a_evcut->fill_hist(this_j.cisvv2, scaled_value*events_MC[Nth].event_weight, Nth);
						sh_THINBtagjet_pt_a_evcut->fill_hist(this_j.p4.Pt(), scaled_value*events_MC[Nth].event_weight, Nth);
						sh_THINBtagjet_eta_a_evcut->fill_hist(this_j.p4.Eta(), scaled_value*events_MC[Nth].event_weight, Nth);
					}
					
					if (N_candidate_jet>4){break;}

				}
			}
			

		}
	
		// event loop : get size of fraction ratio matrix
		int max_N_of_candidate_bjet=0;
		int max_N_of_candidate_nbjet=0;
		
		cout<<endl<<"=====================in get size:==============================="<<endl;
		
		
	  	for (Long64_t Nth = 0; Nth < Nev_MC; Nth++)
	    	{
			//cout<<endl<<"=====================event:"<<Nth<<"==============================="<<endl;
	
			// matching gen par and THINjet

			double scaled_value=get_scaled_value(Nth);
			
			
			
			Nev_MC_scaled+=scaled_value*events_MC[Nth].event_weight;


			if(!is_pass_ev_selection(events_MC[Nth])){continue;}

			int N_of_THINjet=0;
			int N_of_candidate_bjet=0;
			int N_of_candidate_nbjet=0;

			vector<THINjet>  this_THINjets=events_MC[Nth].THINjets_deoverlap_with_lepton;
			for (int i=0;i<this_THINjets.size();i++)
			{	
			
				THINjet this_j=this_THINjets[i];
				N_of_THINjet++;
				if ( this_j.hadron_flavor==5){ N_of_candidate_bjet++;}

				
			}
		
		
		//	 test for add event weight factor
	
			if (Nth>=Nev_index_ZZ) {N_available_event_ZZ_scaled+=scaled_value*events_MC[Nth].event_weight;}
			else if(Nth>=Nev_index_WZ) {N_available_event_WZ_scaled+=scaled_value*events_MC[Nth].event_weight;}
			else if(Nth>=Nev_index_WW) {N_available_event_WW_scaled+=scaled_value*events_MC[Nth].event_weight; }
			//else if(Nth>=Nev_index_WJ) {N_available_event_WJ_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_WJ100to200) {N_available_event_WJ100to200_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_WJ200to400) {N_available_event_WJ200to400_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_WJ400to600) {N_available_event_WJ400to600_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_WJ600to800) {N_available_event_WJ600to800_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_WJ800to1200) {N_available_event_WJ800to1200_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_WJ1200to2500) {N_available_event_WJ1200to2500_scaled+=scaled_value*events_MC[Nth].event_weight; }			
			else if(Nth>=Nev_index_QCD2000toInf) {N_available_event_QCD2000toInf_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_QCD1500to2000) {N_available_event_QCD1500to2000_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_QCD1000to1500) {N_available_event_QCD1000to1500_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_QCD700to1000) {N_available_event_QCD700to1000_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_QCD500to700) {N_available_event_QCD500to700_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_QCD300to500) {N_available_event_QCD300to500_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_QCD200to300) {N_available_event_QCD200to300_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_QCD100to200) {N_available_event_QCD100to200_scaled+=scaled_value*events_MC[Nth].event_weight; }			
			else if(Nth>=Nev_index_DY_inclusive) {N_available_event_DY_inclusive_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_ST_tW_top) {N_available_event_ST_tW_top_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_ST_tW_antitop) {N_available_event_ST_tW_antitop_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_ST_t_channel_top) {N_available_event_ST_t_channel_top_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_ST_t_channel_antitop) {N_available_event_ST_t_channel_antitop_scaled+=scaled_value*events_MC[Nth].event_weight; }
			else if(Nth>=Nev_index_ST_s_channel) {N_available_event_ST_s_channel_scaled+=scaled_value*events_MC[Nth].event_weight;}
			else if(Nth>=Nev_index_ttbar) {N_available_event_ttbar_scaled+=scaled_value*events_MC[Nth].event_weight; }


			N_of_candidate_nbjet=N_of_THINjet-N_of_candidate_bjet;
			if (max_N_of_candidate_bjet<N_of_candidate_bjet){ max_N_of_candidate_bjet=N_of_candidate_bjet; }
			if (max_N_of_candidate_nbjet<N_of_candidate_nbjet){ max_N_of_candidate_nbjet=N_of_candidate_nbjet; }
		

		}
		max_true_b=max_N_of_candidate_bjet;	
		max_true_nb=max_N_of_candidate_nbjet;
		/*	
		Nk_from_rec = vector<Int_t>(max_true_b+max_true_nb+1,0);
		Nt_from_rec = vector<Int_t>(max_true_b+max_true_nb+1,0);
		*/
		Nk_from_rec = vector<double>(50,0);
		Nt_from_rec = vector<double>(50,0);
		
	
		cout<<"max_true_b="<<max_true_b<<", max_true_nb="<<max_true_nb<<";    ";
		vector< vector<double> > mask(max_N_of_candidate_bjet+1, vector<double>(max_N_of_candidate_nbjet+1,0));
		fraction_ratio=mask;

		N_available_event_MC_scaled=
		N_available_event_ttbar_scaled+
		N_available_event_ST_s_channel_scaled+	
		N_available_event_ST_t_channel_antitop_scaled+
		N_available_event_ST_t_channel_top_scaled+
		N_available_event_ST_tW_antitop_scaled+
		N_available_event_ST_tW_top_scaled+
		N_available_event_DY_inclusive_scaled+
		N_available_event_QCD100to200_scaled+
		N_available_event_QCD200to300_scaled+
		N_available_event_QCD300to500_scaled+
		N_available_event_QCD500to700_scaled+
		N_available_event_QCD700to1000_scaled+
		N_available_event_QCD1000to1500_scaled+
		N_available_event_QCD1500to2000_scaled+
		N_available_event_QCD2000toInf_scaled+
		//N_available_event_WJ_scaled+
		N_available_event_WJ100to200_scaled+
		N_available_event_WJ200to400_scaled+
		N_available_event_WJ400to600_scaled+
		N_available_event_WJ600to800_scaled+
		N_available_event_WJ800to1200_scaled+
		N_available_event_WJ1200to2500_scaled+	
		N_available_event_WW_scaled+
		N_available_event_WZ_scaled+
		N_available_event_ZZ_scaled;

		// event loop : get fraction ratio matrix
	  	for (Long64_t Nth = 0; Nth < Nev_MC; Nth++)
	    	{
			//cout<<endl<<"=====================event:"<<Nth<<"==============================="<<endl;
	
			if(!is_pass_ev_selection(events_MC[Nth])){continue;}

			int N_of_THINjet=0;
			int N_of_candidate_bjet=0;
			int N_of_candidate_nbjet=0;


			vector<THINjet>  this_THINjets=events_MC[Nth].THINjets_deoverlap_with_lepton;
			for (int i=0;i<this_THINjets.size();i++)
			{	

				THINjet this_j=this_THINjets[i];
				N_of_THINjet++;
				if ( this_j.hadron_flavor==5){ N_of_candidate_bjet++;}			
			}

			N_of_candidate_nbjet=N_of_THINjet-N_of_candidate_bjet;
			//fraction_ratio[N_of_candidate_bjet][N_of_candidate_nbjet]+=(1/double(N_available_event_MC));


			double scaled_value=get_scaled_value(Nth);

			fraction_ratio[N_of_candidate_bjet][N_of_candidate_nbjet]+=((scaled_value*events_MC[Nth].event_weight)/N_available_event_MC_scaled);
		}	






		// event loop : counting b-tagging for MC
	  	for (Long64_t Nth = 0; Nth < Nev_MC; Nth++)
	    	{
			//cout<<endl<<"=====================event:"<<Nth<<"==============================="<<endl;
	
			// to find scaled value	
			double scaled_value=get_scaled_value(Nth);


			// draw hist before event selection
			
			int gen_b_this_ev=0;
			for(int i=0;i<events_MC[Nth].gen_pars.size();i++)
			{
				if (is_wanted_b(events_MC[Nth].gen_pars[i]))
				{ gen_b_this_ev++;}
			}
	

			if(!is_pass_ev_selection(events_MC[Nth])){continue;}

			// draw hist in event passed event selection
			// draw rec Tmass
			double this_Tmass = get_Tmass(events_MC[Nth]);
			if (this_Tmass!=-1){ sh_rec_Tmass_a_evcut->fill_hist(this_Tmass,scaled_value*events_MC[Nth].event_weight, Nth);}

			vector<THINjet> this_THINjets=events_MC[Nth].THINjets_deoverlap_with_lepton;
		

			// proccess for rec jet 
		
			// counter
			Int_t N_btagging=0;
		
			int N_THINjet_this_ev=0;
			for (int i=0;i<this_THINjets.size();i++)
			{

				THINjet this_j=this_THINjets[i];

				total_kk_rec+=scaled_value*events_MC[Nth].event_weight;
				N_THINjet_this_ev++;
				total_k_rec+=scaled_value*events_MC[Nth].event_weight;
				if ( btagging(this_j))
				{ 
					N_btagging++;  total_btag_pass_cut_rec+=scaled_value*events_MC[Nth].event_weight;
				}
				else 
				{
					total_nbtag_pass_cut_rec+=scaled_value*events_MC[Nth].event_weight;	
				}			
				
			}
			Nt_from_rec[N_btagging]+=scaled_value*events_MC[Nth].event_weight;
			Nk_from_rec[N_THINjet_this_ev]+=scaled_value*events_MC[Nth].event_weight;
	
		}
		
		 
		cout<<"Nev_MC_scaled="<<Nev_MC_scaled<<endl;    	
		cout<<"N_available_event_MC_scaled ="<<N_available_event_MC_scaled<<endl;
		cout<<"N_available_event_MC_scaled/Nev_MC_scaled ="<<double(N_available_event_MC_scaled)/double(Nev_MC_scaled)<<endl;



		cout<<"=========================  gen level  ==========================================="<<endl<<endl;
		cout<<"total_kk_gen="<<gen_j<<"(b="<<gen_b<<", nb="<<gen_nb<<"),"<<endl;
		cout<<"total_k_gen="<<gen_j_pass_cut<<"(b="<<gen_b_pass_cut<<", nb="<<gen_nb_pass_cut<<")"<<endl;
		cout<<"total_t="<<gen_j_pass_cut_btagging<<"(b="<<gen_b_pass_cut_btagging<<", nb="<<gen_nb_pass_cut_btagging<<")"<<endl;
		cout<<endl;
		true_eff=(double)gen_b_pass_cut_btagging/(double)gen_b_pass_cut;
		true_mr=(double)gen_nb_pass_cut_btagging/(double)gen_nb_pass_cut;
		
		cout<<"=========================  rec level(analysis)  ==========================================="<<endl<<endl;
		cout<<"total_kk_rec="<<total_kk_rec<<endl;
		cout<<"total_k_rec="<<total_k_rec<<endl;
		cout<<"total_btag_pass_cut_rec="<<total_btag_pass_cut_rec<<endl;
		cout<<"total_nbtag_pas_cut="<<total_nbtag_pass_cut_rec<<endl;
		cout<<endl;

		cout<<"=========================  Nk_from_rec  ==========================================="<<endl<<endl;
		double sum_of_Nk_from_rec=0;
		double N_of_jet_from_Nk_from_rec=0;
		for (int i =0;i<Nk_from_rec.size();i++)
		{ 
			cout<<"Nk_from_rec["<<i<<"]:"<<Nk_from_rec[i]<<endl;
			sum_of_Nk_from_rec+=Nk_from_rec[i];
			N_of_jet_from_Nk_from_rec+=Nk_from_rec[i]*i;
		}
		cout<<"sum_of_Nk_from_rec="<<sum_of_Nk_from_rec<<endl;
		cout<<"N_of_jet_from_Nk_from_rec="<<N_of_jet_from_Nk_from_rec<<endl;
		cout<<endl;
	
	
		cout<<"=========================  Nt_from_rec  ==========================================="<<endl<<endl;
		double sum_of_Nt_from_rec=0;
		double N_of_jet_from_Nt_from_rec=0;
		for (int i =0;i<Nt_from_rec.size();i++)
		{ 
			cout<<"Nt_from_rec["<<i<<"]:"<<Nt_from_rec[i]<<endl;
			sum_of_Nt_from_rec+=Nt_from_rec[i];
			N_of_jet_from_Nt_from_rec+=Nt_from_rec[i]*i;
		}
		cout<<"sum_of_Nt_from_rec="<<sum_of_Nt_from_rec<<endl;
		cout<<"N_of_jet_from_Nt_from_rec="<<N_of_jet_from_Nt_from_rec<<endl;
		cout<<endl;

/*
		cout<<"===================  exp value of Nt=get_expvalue_of_Nt( true_eff, true_mr, i)  =============="<<endl<<endl;
		double sum_of_exp_value_of_Nt=0;
		double N_of_jet_exp_value_of_Nt=0;	
		for ( int i=0;i< max_true_b+max_true_nb+1;i++)
		{
			double this_exp_value_of_Nt=get_expvalue_of_Nt( true_eff, true_mr, i);
			cout<<"exp value of number of event that have "<<i<<" jet pass cut & btagging:"<<this_exp_value_of_Nt<<endl;
			sum_of_exp_value_of_Nt+=this_exp_value_of_Nt;
			N_of_jet_exp_value_of_Nt+=this_exp_value_of_Nt*i;
		}
		cout<<"sum_of_exp_value_of_Nt="<<sum_of_exp_value_of_Nt<<endl;
		cout<<"N_of_jet_exp_value_of_Nt="<<N_of_jet_exp_value_of_Nt<<endl;
		cout<<endl;
*/





		//get btagging 2p for rec
		// par_for_minuit(par_number,par_name, par_starting_value, step, min, max)
		double this_step =0.001;	
		par_for_minuit eff_rec(0,"eff_rec", 0.6, this_step ,0+this_step,1-this_step);
		par_for_minuit mr_rec(1,"mr_rec", 0.2, this_step ,0+this_step,1-this_step);	
		vector<par_for_minuit> pars_btagging2p_rec;
		pars_btagging2p_rec.push_back(eff_rec);
		pars_btagging2p_rec.push_back(mr_rec);
		vector<par_err_pair_for_minuit> results_btagging2p_rec=get_min_par_btagging2p_rec(pars_btagging2p_rec);
		best_eff_rec=results_btagging2p_rec[0].v;
		best_mr_rec=results_btagging2p_rec[1].v;

		cout<<"============================  fraction ratio  ========================================="<<endl<<endl;
		double sum_of_fraction_ratio=0;
		double N_of_jet_from_fraction_ratio=0;
		for (int i=0;i<fraction_ratio.size();i++)
		{
		
			for (int j=0;j<fraction_ratio[i].size();j++)
			{ 
				if (j<=14) // j up to 14, fraction_ratio is too small to be neglected
				{ cout<<"fraction_ratio["<<i<<"]["<<j<<"]="<<fraction_ratio[i][j]<<"; ";}
				sum_of_fraction_ratio+=fraction_ratio[i][j];	
				N_of_jet_from_fraction_ratio+=(i+j)*fraction_ratio[i][j]*N_available_event_MC_scaled;
			}
			cout<<endl;
		
		}
		cout<<"sum_of_fraction_ratio="<<sum_of_fraction_ratio<<endl;
		cout<<"N_of_jet_from_fraction_ratio="<<N_of_jet_from_fraction_ratio<<endl;
		cout<<endl;
		
/*
		cout<<"=========================  2p in gen level  ==========================================="<<endl<<endl;
		cout<<"eff="<<true_eff<<endl;
		cout<<"mr="<<true_mr<<endl;
		cout<<endl;
*/
	
		cout<<"=========================  rec result  ==========================================="<<endl<<endl;
		cout<<",best_eff_rec="<<best_eff_rec<<"+-"<<results_btagging2p_rec[0].err<<endl;
		cout<<",best_mr_rec="<<best_mr_rec<<"+-"<<results_btagging2p_rec[1].err<<endl;
		cout<<endl;
/*
		cout<<"=========================  ratio rec/gen  ==========================================="<<endl<<endl;
		cout<<",best_eff_rec/true_eff="<<best_eff_rec/true_eff<<endl;
		cout<<",best_mr_rec/true_mr="<<best_mr_rec/true_mr<<endl;
		cout<<endl;
*/		
		zprime_sample_result.add_result(Nev_MC_scaled, N_available_event_MC_scaled, true_eff, true_mr, best_eff_rec, results_btagging2p_rec[0].err, best_mr_rec,results_btagging2p_rec[1].err );
		
		
		
/*	
		cout<<"=========================  data result  ==========================================="<<endl<<endl;
		cout<<",best_eff_data="<<best_eff_data<<"+-"<<results_btagging2p_data[0].err<<endl;
		cout<<",best_mr_data="<<best_mr_data<<"+-"<<results_btagging2p_data[1].err<<endl;
		cout<<endl;
	
		cout<<"=========================  ratio data/MC  ==========================================="<<endl<<endl;
		cout<<",best_eff_data/best_eff_rec="<<best_eff_data/best_eff_rec<<endl;
		cout<<",best_mr_data/best_mr_rec="<<best_mr_data/best_mr_rec<<endl;
		cout<<endl;

		cout<<"=================  counting # of THINjets before/after deoverlap with good lepton  ======================"<<endl<<endl;
		cout<<"N_THINjet_before_deoverlap:"<<N_THINjet_before_deoverlap<<endl;
		cout<<"N_THINjet_after_deoverlap:"<<N_THINjet_after_deoverlap<<endl;
		cout<<"fraction:"<<((double)N_THINjet_before_deoverlap-(double)N_THINjet_after_deoverlap)/(double)N_THINjet_before_deoverlap<<endl;
*/

		
				
		
		
	//*/

	/*
	  	//event loop : analyzer
	  	for (Long64_t Nth = 0; Nth < Nev_MC_pass_gen_selection; Nth++)
	    	{
			cout<<"=====================event:"<<Nth<<"==============================="<<endl;
		


	      		inner_data_for_MC.GetEntry(Nth);

			// genpar var
			Int_t *genParId=inner_data_for_MC.GetPtrInt("genParId");
			Int_t *genMomParId=inner_data_for_MC.GetPtrInt("genMomParId");
			Int_t nGenPar=inner_data_for_MC.GetInt("nGenPar");
			Int_t *genParSt=inner_data_for_MC.GetPtrInt("genParSt");
			Int_t nGenJet=inner_data_for_MC.GetInt("nGenJet");
			Int_t *genMo1=inner_data_for_MC.GetPtrInt("genMo1");

		

	      		// THINjet var
	      		Int_t THINnJet=inner_data_for_MC.GetInt("THINnJet");
	      		Float_t *THINjetPt=inner_data_for_MC.GetPtrFloat("THINjetPt");  
	      		Float_t *THINjetEta=inner_data_for_MC.GetPtrFloat("THINjetEta"); 
	      		Float_t *THINjetPhi=inner_data_for_MC.GetPtrFloat("THINjetPhi"); 
	     		Float_t *THINjetEn=inner_data_for_MC.GetPtrFloat("THINjetEn");
	     		Float_t *THINjetMass=inner_data_for_MC.GetPtrFloat("THINjetMass");
	      		Float_t *THINjetCSV=inner_data_for_MC.GetPtrFloat("THINjetCSV");

	  		Float_t *THINjetNHadEF=inner_data_for_MC.GetPtrFloat("THINjetNHadEF");
	  		Float_t *THINjetNEmEF=inner_data_for_MC.GetPtrFloat("THINjetNEmEF");
	  		Float_t *THINjetCHadEF=inner_data_for_MC.GetPtrFloat("THINjetCHadEF");
	  		Float_t *THINjetCEmEF=inner_data_for_MC.GetPtrFloat("THINjetCEmEF");
	  		Float_t *THINjetCMulti=inner_data_for_MC.GetPtrFloat("THINjetCMulti");;                    
	      		// ELE var
	      		Int_t nEle =inner_data_for_MC.GetInt("nEle");
	      		Float_t *elePt=inner_data_for_MC.GetPtrFloat("elePt");
	   		Float_t *eleEta=inner_data_for_MC.GetPtrFloat("eleEta");
	   		Float_t *elePhi=inner_data_for_MC.GetPtrFloat("elePhi");
	   		Float_t *eleEnergy=inner_data_for_MC.GetPtrFloat("eleEnergy");
	   		Float_t *eleM=inner_data_for_MC.GetPtrFloat("eleM");

			Float_t *eleDelEtaIn=inner_data_for_MC.GetPtrFloat("eleDelEtaIn");
			Float_t *eleDelPhiIn=inner_data_for_MC.GetPtrFloat("eleDelPhiIn");
			Float_t *eleSigIhIh=inner_data_for_MC.GetPtrFloat("eleSigIhIh");
			Float_t *eleHoE=inner_data_for_MC.GetPtrFloat("eleHoE");
			Float_t *eleDxy=inner_data_for_MC.GetPtrFloat("eleDxy");
			Float_t *eleDz=inner_data_for_MC.GetPtrFloat("eleDz");
			Float_t *eleEoverP=inner_data_for_MC.GetPtrFloat("eleEoverP");
			Float_t *eleCorrPfIso=inner_data_for_MC.GetPtrFloat("eleCorrPfIso");
			Int_t *elePassConv=inner_data_for_MC.GetPtrInt("elePassConv");
			Float_t *eleMissingHits=inner_data_for_MC.GetPtrFloat("eleMissingHits");

			// var

			Int_t N_pass_btagging=0;

			//cout<<"N_of_jet:"<<AK5nJet<<endl;


			// genpar OPRATION
			for (Int_t i=0;i<nGenPar;i++)
			{
				if (abs(genParId[i])<=8 && abs(genParId[i])>0 && genParSt[i]==3)
				{
					cout<<"*"<<i<<",";
				}
			
				if (abs(genParId[i])==5)
				{
					cout<<"Mo:"<<genMomParId[i]<<",";
				}
			
			}

			//cout<<"nGenJet:"<<nGenJet<<","<<"AK5nJet:"<<AK5nJet<<","<<"THINnJet:"<<THINnJet;
		
			for (Int_t i=0;i<nGenJet;i++)
			{
				if (abs(genParId[i])==5)
				{
				}
			
			}







	      		// jet OPERATION
	      		for (Int_t i=0;i<AK5nJet;i++)
			{
	      			AK5jet jet_i(AK5jetPt[i],AK5jetEta[i],AK5jetPhi[i],AK5jetEn[i],AK5jetMass[i],AK5jetCSV[i],
					AK5jetNHadEF[i],AK5jetNEmEF[i],AK5jetCHadEF[i],AK5jetCEmEF[i],AK5jetCMulti[i]);
			
		    		for (Int_t j =0;j<i;j++)
		      		{
	      				AK5jet jet_j(AK5jetPt[j],AK5jetEta[j],AK5jetPhi[j],AK5jetEn[j],AK5jetMass[j],AK5jetCSV[j],
					AK5jetNHadEF[j],AK5jetNEmEF[j],AK5jetCHadEF[j],AK5jetCEmEF[j],AK5jetCMulti[j]);
					// jet loose ID
					if(is_pass_jet_loose_id(jet_i) && is_pass_jet_loose_id(jet_j)){}
					else {continue;}
					TLorentzVector jet_i_lv;
					jet_i_lv.SetPtEtaPhiE(jet_i.pt, jet_i.eta, jet_i.phi, jet_i.en);
					TLorentzVector jet_j_lv;
					jet_j_lv.SetPtEtaPhiE(jet_j.pt, jet_j.eta, jet_j.phi, jet_j.en);

					TLorentzVector H_lv;
					H_lv = jet_i_lv+jet_j_lv;
					H_IM->Fill(H_lv.M());

					// jet CSV
					if (is_pass_csvl(jet_i) && is_pass_csvl(jet_j) )
					{ H_IM_csvl->Fill(H_lv.M());}
					if (is_pass_csvm(jet_i) && is_pass_csvm(jet_j) )
					{ H_IM_csvm->Fill(H_lv.M());}


					
		        
		      		}


			}




	     		// ELE OPERATION
	      		
			for (Int_t i=0;i<nEle;i++)
			{
	      			ele ele_i(elePt[i],eleEta[i],elePhi[i],eleEnergy[i],eleM[i],
					eleDelEtaIn[i],eleDelPhiIn[i],eleSigIhIh[i],eleHoE[i],
					eleDxy[i],eleDz[i],eleEoverP[i],eleCorrPfIso[i],elePassConv[i],eleMissingHits[i]);

				if(!is_pass_ele_loose_id(ele_i)){continue;}
				TLorentzVector ele_i_lv;
				ele_i_lv.SetPtEtaPhiE(ele_i.pt, ele_i.eta, ele_i.phi, ele_i.en);

				for(Int_t j=0;j<1;j++)
				{
	      				ele ele_j(elePt[j],eleEta[j],elePhi[j],eleEnergy[j],eleM[j],
						eleDelEtaIn[j],eleDelPhiIn[j],eleSigIhIh[j],eleHoE[j],
						eleDxy[j],eleDz[j],eleEoverP[j],eleCorrPfIso[j],elePassConv[j],eleMissingHits[j]);

					if(!is_pass_ele_loose_id(ele_j)){continue;}
					TLorentzVector ele_j_lv;
					ele_j_lv.SetPtEtaPhiE(ele_j.pt, ele_j.eta, ele_j.phi, ele_j.en);
					TLorentzVector Zee_lv;
					Zee_lv = ele_i_lv+ele_j_lv;
					cout<<"Zee_IM:"<<Zee_lv.M()<<endl;
					Zee_IM->Fill(Zee_lv.M());


				}
	
			}

		
			cout<<endl;

		
	    	}
	*/
	

	// insignificant var
	// ele var (insignificant var)
		
	sh_ele_pt_a_evcut->save_sh( "pic_13Tev_MC", "ele_pt_a_evcut");		
	sh_ele_eta_a_evcut->save_sh( "pic_13Tev_MC", "ele_eta_a_evcut"); 
	sh_ele_dEtaAtVtx_a_evcut->save_sh( "pic_13Tev_MC", "ele_dEtaAtVtx_a_evcut"); 
	sh_ele_dPhiAtVtx_a_evcut->save_sh( "pic_13Tev_MC", "ele_dPhiAtVtx_a_evcut"); 
	sh_ele_D0_a_evcut->save_sh( "pic_13Tev_MC", "ele_D0_a_evcut"); 
	sh_ele_EtaseedAtVtx_a_evcut->save_sh( "pic_13Tev_MC", "ele_EtaseedAtVtx_a_evcut"); 
	sh_ele_HoverE_a_evcut->save_sh( "pic_13Tev_MC", "ele_HoverE_a_evcut");
	sh_ele_MissHits_a_evcut->save_sh( "pic_13Tev_MC", "ele_MissHits_a_evcut");
	sh_ele_MiniIso_a_evcut->save_sh( "pic_13Tev_MC", "ele_MiniIso_a_evcut");
	sh_ele_SigmaIEtaIEta_a_evcut->save_sh( "pic_13Tev_MC", "ele_SigmaIEtaIEta_a_evcut");
		
	//mu var (insignificant var)
	
	sh_mu_pt_a_evcut->save_sh( "pic_13Tev_MC", "mu_pt_a_evcut");	
	sh_mu_eta_a_evcut->save_sh( "pic_13Tev_MC", "mu_eta_a_evcut"); 
	sh_mu_dxy_a_evcut->save_sh( "pic_13Tev_MC", "mu_dxy_a_evcut"); 
	sh_mu_dz_a_evcut->save_sh( "pic_13Tev_MC", "mu_dz_a_evcut"); 
	sh_mu_Hits_a_evcut->save_sh( "pic_13Tev_MC", "mu_Hits_a_evcut"); 
	sh_mu_Matches_a_evcut->save_sh( "pic_13Tev_MC", "mu_Matches_a_evcut");
	sh_mu_MiniIso_a_evcut->save_sh( "pic_13Tev_MC", "mu_MiniIso_a_evcut"); 
	sh_mu_PixelHits_a_evcut->save_sh( "pic_13Tev_MC", "mu_PixelHits_a_evcut"); 
	sh_mu_TrkLayers_a_evcut->save_sh( "pic_13Tev_MC", "mu_TrkLayers_a_evcut");
	sh_mu_TrkPt_a_evcut->save_sh( "pic_13Tev_MC", "mu_TrkPt_a_evcut"); 
		
		
		
				

	// jet var , Bjet var, Btag jet var

	sh_THINjet_CEmEF_a_evcut->save_sh("pic_13Tev_MC","THINjet_CEmEF_a_evcut");
	sh_THINjet_CHadEF_a_evcut->save_sh("pic_13Tev_MC","THINjet_CHadEF_a_evcut");
	sh_THINjet_NEmEF_a_evcut->save_sh("pic_13Tev_MC","THINjet_NEmEF_a_evcut");
	sh_THINjet_NHadEF_a_evcut->save_sh("pic_13Tev_MC","THINjet_NHadEF_a_evcut");
	sh_THINjet_PhoEF_a_evcut->save_sh("pic_13Tev_MC","THINjet_PhoEF_a_evcut");				
	sh_THINjet_cisvv2_a_evcut->save_sh(  "pic_13Tev_MC", "THINjet_cisvv2_a_evcut");
	sh_THINjet_pt_a_evcut->save_sh(  "pic_13Tev_MC", "THINjet_pt_a_evcut");
	sh_THINjet_eta_a_evcut->save_sh(  "pic_13Tev_MC", "THINjet_eta_a_evcut");

	sh_THINBjet_CEmEF_a_evcut->save_sh("pic_13Tev_MC","THINBjet_CEmEF_a_evcut");
	sh_THINBjet_CHadEF_a_evcut->save_sh("pic_13Tev_MC","THINBjet_CHadEF_a_evcut");
	sh_THINBjet_NEmEF_a_evcut->save_sh("pic_13Tev_MC","THINBjet_NEmEF_a_evcut");
	sh_THINBjet_NHadEF_a_evcut->save_sh("pic_13Tev_MC","THINBjet_NHadEF_a_evcut");
	sh_THINBjet_PhoEF_a_evcut->save_sh("pic_13Tev_MC","THINBjet_PhoEF_a_evcut");			
	sh_THINBtagjet_cisvv2_a_evcut->save_sh(  "pic_13Tev_MC", "THINBtagjet_cisvv2_a_evcut");
	sh_THINBtagjet_pt_a_evcut->save_sh(  "pic_13Tev_MC", "THINBtagjet_pt_a_evcut");
	sh_THINBtagjet_eta_a_evcut->save_sh(  "pic_13Tev_MC", "THINBtagjet_eta_a_evcut");

	sh_rec_Tmass_a_evcut->save_sh(  "pic_13Tev_MC", "rec_Tmass_a_evcut");
			


	
	}
	//gerrors(zprime_sample_result, "eff");
	//gerrors(zprime_sample_result, "mr");
 	    
}
