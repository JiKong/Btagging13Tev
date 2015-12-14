//root
#include <map>
#include <vector>
#include <string>
#include <iostream>
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
	Long64_t Nev_data=0;
	int N_available_event_data=0;

	// gen level 
	vector< vector<double> > fraction_ratio;
	int max_true_b=0;
	int max_true_nb=0;
	
	double true_eff=0;
	double true_mr=0;

	int gen_j=0;
	int gen_b=0;
	int gen_nb=0;
	int gen_j_pass_cut=0;
	int gen_b_pass_cut=0;
	int gen_nb_pass_cut=0;
	int gen_j_pass_cut_btagging=0;
	int gen_b_pass_cut_btagging=0;
	int gen_nb_pass_cut_btagging=0;

	// rec level
	vector<Int_t> Nt_from_rec(1,0) ; // # of event that have t b-tagging
	vector<Int_t> Nk_from_rec(1,0) ; // # of event that have k jet

	double best_eff_rec=0;
	double best_mr_rec=0;

	int total_kk_rec=0;
	int total_k_rec=0;
	int total_btag_pass_cut_rec=0;
	int total_nbtag_pass_cut_rec=0;
	

	// counting
	int N_THINjet_before_deoverlap=0;
	int N_THINjet_after_deoverlap=0;



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

	h_->SetMarkerColor(1);
	h_->SetMarkerStyle(2);
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

}
*/

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
				if (abs(this_Wmass-80) <=20){find_good_Wqqbar=true;}
				if (abs(this_Wmass-80) <= abs(rec_Wmass-81))
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


//event selection


//evse rec level
bool is_pass_ev_selection_TTbar(event ev)
{

	// trigger path
	//if (!ev.is_pass_mu_trigger){return false;}
	if (!ev.is_pass_ele_trigger){return false;}

	// same as theevent selection on ttbar sample, be cause we will see how ttbar selection works.....
	int channel=0; //1 for ele, 2 for mu
	int ev_flag=0;

	// 1: event must have 1 good muon, isolated, PT>50, |Eta|<2.1
	vector<ele> good_eles=ev.good_eles;
	bool have_available_ele=false;
	//electron channel
	
	for (int i=0;i<good_eles.size();i++)
	{
		ele this_e=good_eles[i];
		if (this_e.p4.Pt()>0 && abs(this_e.p4.Eta())<2.4)
		{ have_available_ele=true; break;}
	}
	if(!have_available_ele)
	{ev_flag=2; return false;}
	
	 
	// muon channel
	vector<mu> good_mus=ev.good_mus;
/*
	bool have_available_muon=false;
	for (int i=0;i<good_mus.size();i++)
	{
		mu this_m=good_mus[i];
		if (this_m.p4.Pt()>50 && abs(this_m.p4.Eta())<2.1)
		{ have_available_muon=true; break;}
	}
	if(!have_available_muon)
	{ev_flag=2; return false;}
*/	


	

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

	if (get_Tmass(ev)==-1){return false;}
	

	return true;
	

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
                    	double term1 = fraction_ratio[i][j]<=0 ? 0 : N_available_event_data * fraction_ratio[i][j];

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
                                		cout<<"N_available_event_data="<<N_available_event_data<<","<<endl;
                                		cout<<"fraction_ratio[i][j]="<<fraction_ratio[i][j]<<","<<endl;
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
		double this_log_like=log_poisson(double(Nt_from_rec[t]), get_expvalue_of_Nt( par[0], par[1] ,t));
                total_log_like += this_log_like;
		if (this_log_like < -10000000)
		{ 
			cout<<"this_log_like < -10000000"<<"when eff="<<par[0]<<",mr="<<par[1]<<",t="<<t;
			cout<<",Nt_from_rec[t]="<<double(Nt_from_rec[t])<<",get_expvalue_of_N="<<get_expvalue_of_Nt( par[0], par[1] ,t)<<endl;
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
void Btagging_in_13Tev_data()
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
	chk_gen_b is_wanted_b = is_wanted_b_ttbar;
	//chk_gen_b is_wanted_b = is_wanted_b_Zprime;

		
		
	
	
	sample_result zprime_sample_result("zprime_sample_result");
	
	for (int sample_i=0;sample_i<1; sample_i++) //collect all small samples in a big sample
	{
		// initialization variables
		N_available_event_data=0;

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
	


		
		// hist set
		// ele var (insignificant var)
	
		TH1D* h_ele_pt_a_evcut = new TH1D("ele_pt_a_evcut","ele_pt_a_evcut",50,0,300);		
		TH1D* h_ele_eta_a_evcut = new TH1D("ele_eta_a_evcut","ele_eta_a_evcut",50,-7,7);
		TH1D* h_ele_dEtaAtVtx_a_evcut = new TH1D("ele_dEtaAtVtx_a_evcut","ele_dEtaAtVtx_a_evcut",50,-0.2,0.2);
		TH1D* h_ele_dPhiAtVtx_a_evcut = new TH1D("ele_dPhiAtVtx_a_evcut","ele_dPhiAtVtx_a_evcut",50,-1,1);
		TH1D* h_ele_D0_a_evcut = new TH1D("ele_D0_a_evcut","ele_D0_a_evcut",50,-2,2);
		TH1D* h_ele_EtaseedAtVtx_a_evcut = new TH1D("ele_EtaseedAtVtx_a_evcut","ele_EtaseedAtVtx_a_evcut",50,-0.2,0.2);
		TH1D* h_ele_HoverE_a_evcut = new TH1D("ele_HoverE_a_evcut","ele_HoverE_a_evcut",50,0,6);
		TH1D* h_ele_MissHits_a_evcut = new TH1D("ele_MissHits_a_evcut","ele_MissHits_a_evcut",7,0,7);
		//TH1D* h_ele_MiniIso_a_evcut = new TH1D("ele_MiniIso_a_evcut","ele_MiniIso_a_evcut",50,0,10);
		TH1D* h_ele_SigmaIEtaIEta_a_evcut = new TH1D("ele_SigmaIEtaIEta_a_evcut","ele_SigmaIEtaIEta_a_evcut",50,0,0.08);
		
		//mu var (insignificant var)
	
		TH1D* h_mu_pt_a_evcut = new TH1D("mu_pt_a_evcut","mu_pt_a_evcut",50,0,300);		
		TH1D* h_mu_eta_a_evcut = new TH1D("mu_eta_a_evcut","mu_eta_a_evcut",50,-7,7);
		TH1D* h_mu_dxy_a_evcut = new TH1D("mu_dxy_a_evcut","mu_dxy_a_evcut",20,-1,1);
		TH1D* h_mu_dz_a_evcut = new TH1D("mu_dz_a_evcut","mu_dz_a_evcut",40,-2,2);		
		TH1D* h_mu_Hits_a_evcut = new TH1D("mu_Hits_a_evcut","mu_Hits_a_evcut",50,0,100);
		TH1D* h_mu_Matches_a_evcut = new TH1D("mu_Matches_a_evcut","mu_Matches_a_evcut",10,0,10);	
		//TH1D* h_mu_MiniIso_a_evcut = new TH1D("mu_MiniIso_a_evcut","mu_MiniIso_a_evcut",50,0,100);
		TH1D* h_mu_PixelHits_a_evcut = new TH1D("mu_PixelHits_a_evcut","mu_PixelHits_a_evcut",20,0,20);
		TH1D* h_mu_TrkLayers_a_evcut = new TH1D("mu_TrkLayers_a_evcut","mu_TrkLayers_a_evcut",25,0,25);
		TH1D* h_mu_TrkPt_a_evcut = new TH1D("mu_TrkPt_a_evcut","mu_TrkPt_a_evcut",50,0,500);
		
	
		//THINjet var
				
		TH1D* h_THINjet_CEmEF_a_evcut = new TH1D("THINjet_CEmEF_a_evcut","THINjet_CEmEF_a_evcut",50,0,1);
		TH1D* h_THINjet_CHadEF_a_evcut = new TH1D("THINjet_CHadEF_a_evcut","THINjet_CHadEF_a_evcut",50,0,1);
		TH1D* h_THINjet_NEmEF_a_evcut = new TH1D("THINjet_NEmEF_a_evcut","THINjet_NEmEF_a_evcut",50,0,1);
		TH1D* h_THINjet_NHadEF_a_evcut = new TH1D("THINjet_NHadEF_a_evcut","THINjet_NHadEF_a_evcut",50,0,1);
		TH1D* h_THINjet_PhoEF_a_evcut = new TH1D("THINjet_PhoEF_a_evcut","THINjet_PhoEF_a_evcut",50,0,1);
							
		TH1D* h_THINjet_cisvv2_a_evcut = new TH1D("THINjet_cisvv2_a_evcut","THINjet_cisvv2_a_evcut",50,0,1);
		TH1D* h_THINjet_pt_a_evcut = new TH1D("THINjet_pt_a_evcut","THINjet_pt_a_evcut",50,0,500);		
		TH1D* h_THINjet_eta_a_evcut = new TH1D("THINjet_eta_a_evcut","THINjet_eta_a_evcut",50,-7,7);
		
		TH1D* h_THINBjet_CEmEF_a_evcut = new TH1D("THINBjet_CEmEF_a_evcut","THINBjet_CEmEF_a_evcut",50,0,1);
		TH1D* h_THINBjet_CHadEF_a_evcut = new TH1D("THINBjet_CHadEF_a_evcut","THINBjet_CHadEF_a_evcut",50,0,1);
		TH1D* h_THINBjet_NEmEF_a_evcut = new TH1D("THINBjet_NEmEF_a_evcut","THINBjet_NEmEF_a_evcut",50,0,1);
		TH1D* h_THINBjet_NHadEF_a_evcut = new TH1D("THINBjet_NHadEF_a_evcut","THINBjet_NHadEF_a_evcut",50,0,1);
		TH1D* h_THINBjet_PhoEF_a_evcut = new TH1D("THINBjet_PhoEF_a_evcut","THINBjet_PhoEF_a_evcut",50,0,1);
		
		TH1D* h_THINBtagjet_cisvv2_a_evcut = new TH1D("THINBtagjet_cisvv2_a_evcut","THINBtagjet_cisvv2_a_evcut",50,0,1);
		TH1D* h_THINBtagjet_pt_a_evcut = new TH1D("THINBtagjet_pt_a_evcut","THINBtagjet_pt_a_evcut",50,0,500);		
		TH1D* h_THINBtagjet_eta_a_evcut = new TH1D("THINBtagjet_eta_a_evcut","THINBtagjet_eta_a_evcut",50,-7,7);

		// rec t mass
		TH1D* h_rec_Tmass_a_evcut = new TH1D("rec_Tmass_a_evcut","rec_Tmass_a_evcut",50,0,300);
	
	
		// root file
	  	//TChain *data = new TChain("tree/treeMaker");
	  	//data->Add("../flattuple.root");// Zprime sample
	  	//data->Add("../TTbar_sample/flattuple_1.root");// TTbar sample from yu
	  	//data->Add(samples[sample_i]);
	  	
	  	//int N_applied_sample=50;
	  	vector<string> input_files;
	  	

	  	
	  	// new version  
		// single ele channel
		// we have  102999907 events from all muon channel sample
		// we have  14580254 events from 1/7 muon channel sample
		// we have 20640253 events from 1/5 muon channel sample

	  	for (int i=1;i<=411;i++)
	  	{
	  		//if(i>205){continue;}
	  		// %5
			if (i%7!=0){continue;}
	  		TString filepath = "/data7/khurana/NCUGlobalTuples/Run2015DFullDataset/crab_SingleElectron-Run2015D-05Oct2015-v1_20151117_2p2fb_SingleEleTextFile/151117_152256/0000/NCUGlobalTuples_"+int_to_string(i)+".root";
	  		TFile f(filepath);
			if (f.IsZombie()) {continue; f.Close();} f.Close();
		
	  		input_files.push_back("/data7/khurana/NCUGlobalTuples/Run2015DFullDataset/crab_SingleElectron-Run2015D-05Oct2015-v1_20151117_2p2fb_SingleEleTextFile/151117_152256/0000/NCUGlobalTuples_"+int_to_string(i)+".root");
	  		
	  	}	  	
	  	for (int i=1;i<=855;i++)
	  	{
	  		//if(i>94){continue;}
			if (i%7!=0){continue;}
	  		TString filepath ="/data7/khurana/NCUGlobalTuples/Run2015DFullDataset/crab_SingleElectron-Run2015D-PromptReco-V420151117_2p2fb_SingleEleTextFile/151117_152534/0000/NCUGlobalTuples_"+int_to_string(i)+".root";	  		
	  		TFile f(filepath);
			if (f.IsZombie()) {continue; f.Close();} f.Close();
			
	  		input_files.push_back("/data7/khurana/NCUGlobalTuples/Run2015DFullDataset/crab_SingleElectron-Run2015D-PromptReco-V420151117_2p2fb_SingleEleTextFile/151117_152534/0000/NCUGlobalTuples_"+int_to_string(i)+".root");
	  	}


		// single muon channel
		// we have 71641966 events from all muon channel sample
		// we have 11879377 events from 1/6 muon channel sample
		// we have 14408988 events from 1/5 muon channel sample

/*
	  	for (int i=1;i<=411;i++)
	  	{
	  		//if(i>205){continue;}
	  		if (i%5!=0){continue;}
	  		TString filepath ="/data7/syu/NCUGlobalTuples/Run2015D/9b33d00/SingleMuon/crab_SingleMuon-Run2015D-05Oct2015-v1_20151119_2p2fb_SingleMuTextFile/151120_125737/0000/NCUGlobalTuples_"+int_to_string(i)+".root";	  		
	  		TFile f(filepath);
			if (f.IsZombie()) {continue; f.Close();} f.Close();
			
	  		input_files.push_back("/data7/syu/NCUGlobalTuples/Run2015D/9b33d00/SingleMuon/crab_SingleMuon-Run2015D-05Oct2015-v1_20151119_2p2fb_SingleMuTextFile/151120_125737/0000/NCUGlobalTuples_"+int_to_string(i)+".root");
	  	}
	  	for (int i=1;i<=855;i++)
	  	{
	  		//if(i>94){continue;}
	  		if (i%5!=0){continue;}
	  		TString filepath ="/data7/syu/NCUGlobalTuples/Run2015D/9b33d00/SingleMuon/crab_SingleMuon-Run2015D-PromptReco-V420151119_2p2fb_SingleMuTextFile/151119_213943/0000/NCUGlobalTuples_"+int_to_string(i)+".root";	  		
	  		TFile f(filepath);
			if (f.IsZombie()) {continue; f.Close();} f.Close();
			
	  		input_files.push_back("/data7/syu/NCUGlobalTuples/Run2015D/9b33d00/SingleMuon/crab_SingleMuon-Run2015D-PromptReco-V420151119_2p2fb_SingleMuTextFile/151119_213943/0000/NCUGlobalTuples_"+int_to_string(i)+".root");
	  	}
*/
	  	
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_100.root");
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_200.root");
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_300.root");
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_400.root");
	  	

	  	TreeReader inner_data_for_data(input_files);	  	
	  	Nev_data = inner_data_for_data.GetEntriesFast();


	  	cout<< "Nev_data: "<< Nev_data << endl;
		  	
		//configure
		vector<event> events_data;	
		int n_THINj=0;
	
		//matching

		// event loop : set up data
	  	//for (Long64_t Nth = 0; Nth < Nev_data_before_gen_selection; Nth++)
	  	for (Long64_t Nth = 0; Nth < Nev_data; Nth++)
	    	{
	 		inner_data_for_data.GetEntry(Nth);

			// trigger var
			string* hlt_trigName = inner_data_for_data.GetPtrString("hlt_trigName");
			vector<bool> &hlt_trigResult = *((vector<bool>*) inner_data_for_data.GetPtr("hlt_trigResult"));	
			int string_size = inner_data_for_data.GetPtrStringSize();
			
			// gen_par var
	      		Int_t nGenPar=inner_data_for_data.GetInt("nGenPar");
	      		/*
	      		TClonesArray* genParP4_ = (TClonesArray*)data.GetPtrTObject("genParP4");
	      		TLorentzVector* genParP4=(TLorentzVector*)genParP4_;	   */
	      		TClonesArray* genParP4 = (TClonesArray*) inner_data_for_data.GetPtrTObject("genParP4"); 		
	     			
			Int_t *genParId=inner_data_for_data.GetPtrInt("genParId");	
			Int_t *genParIndex=inner_data_for_data.GetPtrInt("genParIndex");
			Int_t *genMomParId=inner_data_for_data.GetPtrInt("genMomParId");
			Int_t *genMo1=inner_data_for_data.GetPtrInt("genMo1");  //not mom id
			Int_t *genMo2=inner_data_for_data.GetPtrInt("genMo2");  //not mom id
			Int_t *genParSt=inner_data_for_data.GetPtrInt("genParSt"); 



	      		// THINjet var
	      		Int_t THINnJet=inner_data_for_data.GetInt("THINnJet");
	      		
	      		TClonesArray* THINjetP4 = (TClonesArray*) inner_data_for_data.GetPtrTObject("THINjetP4");
	      		Float_t *THINjetCSV=inner_data_for_data.GetPtrFloat("THINjetCSV");
	      	      //	Float_t* THINjetPrunedEn = inner_data_for_data.GetPtrFloat("THINjetPrunedEn");
	      	      //	Int_t* THINjetPassID = inner_data_for_data.GetPtrInt("THINjetPassID");
	      	      //	Float_t* THINjetTau1 = inner_data_for_data.GetPtrFloat("THINjetTau1");
			//	Float_t* THINjetTau2 = inner_data_for_data.GetPtrFloat("THINjetTau2");
			Float_t* THINjetCISVV2 = inner_data_for_data.GetPtrFloat("THINjetCISVV2");
			Int_t * THINjetHadronFlavor=inner_data_for_data.GetPtrInt("THINjetHadronFlavor"); 			
			
			vector<bool> &THINjetPassIDTight = *((vector<bool>*) inner_data_for_data.GetPtr("THINjetPassIDTight"));			
			vector<bool> &THINjetPassIDLoose  = *((vector<bool>*) inner_data_for_data.GetPtr("THINjetPassIDLoose"));


			Float_t *THINjetCEmEF=inner_data_for_data.GetPtrFloat("THINjetCEmEF");
			Float_t *THINjetCHadEF=inner_data_for_data.GetPtrFloat("THINjetCHadEF");	
			Float_t *THINjetNEmEF=inner_data_for_data.GetPtrFloat("THINjetNEmEF");	
			Float_t *THINjetNHadEF=inner_data_for_data.GetPtrFloat("THINjetNHadEF");
			Float_t *THINjetPhoEF=inner_data_for_data.GetPtrFloat("THINjetPhoEF");
	      	      		

			// muon var
	      		Int_t nMu =inner_data_for_data.GetInt("nMu");
	      		TClonesArray* muP4 = (TClonesArray*) inner_data_for_data.GetPtrTObject("muP4");
	      		
	    		
	    		vector<bool> &isGlobalMuon  = *((vector<bool>*) inner_data_for_data.GetPtr("isGlobalMuon"));
	    		vector<bool> &isTrackerMuon  = *((vector<bool>*) inner_data_for_data.GetPtr("isTrackerMuon"));
	    		vector<bool> &isTightMuon  = *((vector<bool>*) inner_data_for_data.GetPtr("isTightMuon"));
	    		vector<bool> &isMediumMuon  = *((vector<bool>*) inner_data_for_data.GetPtr("isMediumMuon"));
	    		vector<bool> &isLooseMuon = *((vector<bool>*) inner_data_for_data.GetPtr("isLooseMuon"));
	    		vector<bool> &isHighPtMuon  = *((vector<bool>*) inner_data_for_data.GetPtr("isHighPtMuon"));
	    		vector<bool> &isSoftMuon  = *((vector<bool>*) inner_data_for_data.GetPtr("isSoftMuon"));
	    		
	    		Float_t* muChHadIso = inner_data_for_data.GetPtrFloat("muChHadIso");
	    		Float_t* muNeHadIso = inner_data_for_data.GetPtrFloat("muNeHadIso");
	    		Float_t* muGamIso = inner_data_for_data.GetPtrFloat("muGamIso");
	    		Float_t* muPUPt = inner_data_for_data.GetPtrFloat("muPUPt");
	    		
/*
	    		Int_t*   muTrkLayers  = inner_data_for_data.GetPtrInt("muTrkLayers");
	    		Int_t*   muPixelHits  = inner_data_for_data.GetPtrInt("muPixelHits");
	    		Int_t*   muHits       = inner_data_for_data.GetPtrInt("muHits");
	    		Int_t*   muMatches    = inner_data_for_data.GetPtrInt("muMatches");
	    		Float_t* mudxy        = inner_data_for_data.GetPtrFloat("mudxy");
	    		Float_t* mudz         = inner_data_for_data.GetPtrFloat("mudz");

	   		//Float_t *muCorrPfIso=inner_data_for_data.GetPtrFloat("muCorrPfIso");				
			//Int_t* muPassID = inner_data_for_data.GetPtrInt("muPassID");
			//Float_t* muCorrTrkIso = inner_data_for_data.GetPtrFloat("muCorrTrkIso");
*/
	      		// ele var
	      		Int_t nEle =inner_data_for_data.GetInt("nEle");
	      		TClonesArray* eleP4 = (TClonesArray*) inner_data_for_data.GetPtrTObject("eleP4");
	      		
	      		vector<bool> &eleIsPassTight  = *((vector<bool>*) inner_data_for_data.GetPtr("eleIsPassTight"));
	      		vector<bool> &eleIsPassMedium  = *((vector<bool>*) inner_data_for_data.GetPtr("eleIsPassMedium"));
	      		vector<bool> &eleIsPassLoose  = *((vector<bool>*) inner_data_for_data.GetPtr("eleIsPassLoose"));
	      		
	      		
/*
			//Float_t *eleDelEtaIn=inner_data_for_data.GetPtrFloat("eleDelEtaIn");
			//Float_t *eleDelPhiIn=inner_data_for_data.GetPtrFloat("eleDelPhiIn");
			//Float_t *eleSigIhIh=inner_data_for_data.GetPtrFloat("eleSigIhIh");
			Float_t *eleHoverE=inner_data_for_data.GetPtrFloat("eleHoverE");
			//Float_t *eleDxy=inner_data_for_data.GetPtrFloat("eleDxy");
			Float_t *eleDz=inner_data_for_data.GetPtrFloat("eleDz");
			Float_t *eleEoverP=inner_data_for_data.GetPtrFloat("eleEoverP");
			//Float_t *eleCorrPfIso=inner_data_for_data.GetPtrFloat("eleCorrPfIso");
			//Int_t *elePassConv=inner_data_for_data.GetPtrInt("elePassConv");
			//Float_t *eleMissingHits=inner_data_for_data.GetPtrFloat("eleMissingHits");

			//Int_t *elePassID=inner_data_for_data.GetPtrInt("elePassID");
			//Float_t *eleUserTrkIso=inner_data_for_data.GetPtrFloat("eleUserTrkIso");
			//Float_t *eleUserCalIso=inner_data_for_data.GetPtrFloat("eleUserCalIso");
			Float_t eleRho=inner_data_for_data.GetFloat("eleRho");
*/

			//met
			Float_t pfMetCorrPt=inner_data_for_data.GetFloat("pfMetCorrPt");
			Float_t pfMetCorrPhi=inner_data_for_data.GetFloat("pfMetCorrPhi");


			event this_ev;

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
					(abs(genParId[i])==5  && abs(genMomParId[i])==6)||
					( abs(genParId[i])<=5  && abs(genParId[i])>=1  && abs(genMomParId[i])==24)||
					( abs(genParId[i])==11  && abs(genMomParId[i])==24)|| 
					( abs(genParId[i])==13  && abs(genMomParId[i])==24)
				)
				{
					/*if ( (abs(genParId[i])>=1 && abs(genParId[i])<=8) || 
					(abs(genParId[i])>=11 && abs(genParId[i])<=18) || 
					(abs(genParId[i])>=21 && abs(genParId[i])<=22)){}
					else{continue;} */
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
				
				for (int elei=0;elei<this_ev.good_eles.size();elei++)
				{	
					ele this_ele = this_ev.good_eles[elei];
					if (delta_R(this_j, this_ele)<0.4){overlap_with_lepton=true;break;}
				}
				
				/*
				for (int mui=0;mui<this_ev.good_mus.size();mui++)
				{
					mu this_mu = this_ev.good_mus[mui];
					if (delta_R(this_j, this_mu)<0.4){overlap_with_lepton=true;break;}
				}
				*/
				
				if (overlap_with_lepton){ continue;}
			
				this_ev.THINjets_deoverlap_with_lepton.push_back(this_j);
				
			}
			
			events_data.push_back(this_ev);
			
			// draw insignificant var (these insignificant var did not saved in class event, excluding pt & eta )
			// // ele (only 1 ele decayed from ttbar)
			
			Float_t *eledEtaAtVtx=inner_data_for_data.GetPtrFloat("eledEtaAtVtx");
			Float_t *eledPhiAtVtx=inner_data_for_data.GetPtrFloat("eledPhiAtVtx");
			Float_t *eleD0=inner_data_for_data.GetPtrFloat("eleD0");
			Float_t *eleEtaseedAtVtx=inner_data_for_data.GetPtrFloat("eleEtaseedAtVtx");
			Float_t *eleHoverE=inner_data_for_data.GetPtrFloat("eleHoverE");
			Int_t *eleMissHits=inner_data_for_data.GetPtrInt("eleMissHits");
			//Float_t *eleMiniIso=inner_data_for_data.GetPtrFloat("eleMiniIso");
			Float_t *eleSigmaIEtaIEta=inner_data_for_data.GetPtrFloat("eleSigmaIEtaIEta");

			if (is_pass_ev_selection(this_ev))
			{
				bool find_candidate_ele=false;
				for (int i=0;i<nEle;i++)
				{
					if (eleIsPassLoose[i] && (*((TLorentzVector*)eleP4->At(i))).Pt()>0 && abs( (*((TLorentzVector*)eleP4->At(i))).Eta())<2.4)
					{find_candidate_ele=true;}
					else{continue;}

					h_ele_pt_a_evcut->Fill( (*((TLorentzVector*)eleP4->At(i))).Pt());
			 		h_ele_eta_a_evcut->Fill( (*((TLorentzVector*)eleP4->At(i))).Eta()); 
			 		h_ele_dEtaAtVtx_a_evcut->Fill(eledEtaAtVtx[i]); 
			 		h_ele_dPhiAtVtx_a_evcut->Fill(eledPhiAtVtx[i]); 
		
					h_ele_D0_a_evcut->Fill(eleD0[i]);  
			 		h_ele_EtaseedAtVtx_a_evcut ->Fill(eleEtaseedAtVtx[i]); 
			 		h_ele_HoverE_a_evcut->Fill(eleHoverE[i]); 
			 		h_ele_MissHits_a_evcut->Fill(eleMissHits[i]); 
			 		//h_ele_MiniIso_a_evcut->Fill(eleMiniIso[i]); 
					h_ele_SigmaIEtaIEta_a_evcut->Fill(eleSigmaIEtaIEta[i]); 

					if (find_candidate_ele){break;}
				}	
			}

	
			// // mu (only 1 mu decayed from ttbar)
/*			
			Float_t *mudxy=inner_data_for_data.GetPtrFloat("mudxy");	
			Float_t *mudz=inner_data_for_data.GetPtrFloat("mudz");		
			Int_t *muHits=inner_data_for_data.GetPtrInt("muHits");		
			Int_t *muMatches=inner_data_for_data.GetPtrInt("muMatches");		
			//Float_t *muMiniIso=inner_data_for_data.GetPtrFloat("muMiniIso");		
			Int_t *muPixelHits=inner_data_for_data.GetPtrInt("muPixelHits");		
			Int_t *muTrkLayers=inner_data_for_data.GetPtrInt("muTrkLayers");		
			Float_t *muTrkPt=inner_data_for_data.GetPtrFloat("muTrkPt");

			if (is_pass_ev_selection(this_ev))
			{
				bool find_candidate_mu=false;
				for (int i=0;i<nMu;i++)
				{
					if (isLooseMuon[i] && (*((TLorentzVector*)muP4->At(i))).Pt()>50 && abs( (*((TLorentzVector*)muP4->At(i))).Eta())<2.1)
					{find_candidate_mu=true;}
					else{continue;}
					
					h_mu_pt_a_evcut->Fill( (*((TLorentzVector*)muP4->At(i))).Pt());	
					h_mu_eta_a_evcut->Fill( (*((TLorentzVector*)muP4->At(i))).Eta());
					h_mu_dxy_a_evcut->Fill(mudxy[i]); 
					h_mu_dz_a_evcut->Fill(mudz[i]); 	
					h_mu_Hits_a_evcut->Fill(muHits[i]); 
					h_mu_Matches_a_evcut->Fill(muMatches[i]); 
					//h_mu_MiniIso_a_evcut->Fill(muMiniIso[i]); 
					h_mu_PixelHits_a_evcut->Fill(muPixelHits[i]); 
					h_mu_TrkLayers_a_evcut->Fill(muTrkLayers[i]); 
					h_mu_TrkPt_a_evcut->Fill(muTrkPt[i]); 
					
					if (find_candidate_mu){break;}
				}	
			}
*/

			// THINjet var (only 4 good jet decayed from ttbar)
			
			if (is_pass_ev_selection(this_ev))
			{
				int N_candidate_jet=0;
				for(int i=0;i<this_ev.THINjets_deoverlap_with_lepton.size();i++)
				{
					THINjet this_j = this_ev.THINjets_deoverlap_with_lepton[i];
					if (!is_good_jet(this_j)){continue;}
					else {N_candidate_jet++;}				
					
					h_THINjet_CEmEF_a_evcut->Fill(this_j.CEmEF);
					h_THINjet_CHadEF_a_evcut->Fill(this_j.CHadEF);
					h_THINjet_NEmEF_a_evcut->Fill(this_j.NEmEF);
					h_THINjet_NHadEF_a_evcut->Fill(this_j.NHadEF);
					h_THINjet_PhoEF_a_evcut->Fill(this_j.PhoEF);
									
					h_THINjet_cisvv2_a_evcut->Fill(this_j.cisvv2);
					h_THINjet_pt_a_evcut->Fill(this_j.p4.Pt());
					h_THINjet_eta_a_evcut->Fill(this_j.p4.Eta());
				
					if (btagging(this_j))
					{
						h_THINBjet_CEmEF_a_evcut->Fill(this_j.CEmEF);
						h_THINBjet_CHadEF_a_evcut->Fill(this_j.CHadEF);
						h_THINBjet_NEmEF_a_evcut->Fill(this_j.NEmEF);
						h_THINBjet_NHadEF_a_evcut->Fill(this_j.NHadEF);
						h_THINBjet_PhoEF_a_evcut->Fill(this_j.PhoEF);
					
						h_THINBtagjet_cisvv2_a_evcut->Fill(this_j.cisvv2);
						h_THINBtagjet_pt_a_evcut->Fill(this_j.p4.Pt());
						h_THINBtagjet_eta_a_evcut->Fill(this_j.p4.Eta());
					}
					
					if (N_candidate_jet>4){break;}

				}
			}			
			
	
			
		}

	
		Nk_from_rec = vector<Int_t>(50,0);
		Nt_from_rec = vector<Int_t>(50,0);
		
		//vector< vector<double> > mask(max_N_of_candidate_bjet+1, vector<double>(max_N_of_candidate_nbjet+1,0));
		//fraction_ratio=mask;

		// event loop : get fraction ratio matrix
		vector< vector<double> > mask(7, vector<double>(20,0));
		fraction_ratio = mask;

//fr

fraction_ratio[0][0]=0; fraction_ratio[0][1]=0; fraction_ratio[0][2]=0; fraction_ratio[0][3]=0; fraction_ratio[0][4]=0.0144474; fraction_ratio[0][5]=0.0699306; fraction_ratio[0][6]=0.0919902; fraction_ratio[0][7]=0.0575296; fraction_ratio[0][8]=0.0388094; fraction_ratio[0][9]=0.0225145; fraction_ratio[0][10]=0.012381; fraction_ratio[0][11]=0.00539195; fraction_ratio[0][12]=0.00452473; fraction_ratio[0][13]=0.00250965; fraction_ratio[0][14]=0.0016047; 
fraction_ratio[1][0]=0; fraction_ratio[1][1]=0; fraction_ratio[1][2]=0; fraction_ratio[1][3]=0.0269914; fraction_ratio[1][4]=0.0515879; fraction_ratio[1][5]=0.0464255; fraction_ratio[1][6]=0.0261715; fraction_ratio[1][7]=0.013363; fraction_ratio[1][8]=0.00310383; fraction_ratio[1][9]=0.00296305; fraction_ratio[1][10]=0.00151024; fraction_ratio[1][11]=0.000225092; fraction_ratio[1][12]=0.000101922; fraction_ratio[1][13]=0.00012902; fraction_ratio[1][14]=4.35834e-05; 
fraction_ratio[2][0]=0; fraction_ratio[2][1]=0; fraction_ratio[2][2]=0.0598157; fraction_ratio[2][3]=0.13332; fraction_ratio[2][4]=0.11862; fraction_ratio[2][5]=0.0777534; fraction_ratio[2][6]=0.0357877; fraction_ratio[2][7]=0.0172737; fraction_ratio[2][8]=0.00566128; fraction_ratio[2][9]=0.00268512; fraction_ratio[2][10]=0.0277489; fraction_ratio[2][11]=0.000286117; fraction_ratio[2][12]=8.40537e-05; fraction_ratio[2][13]=5.48169e-05; fraction_ratio[2][14]=2.32807e-05; 
fraction_ratio[3][0]=0; fraction_ratio[3][1]=0.00048398; fraction_ratio[3][2]=0.00283577; fraction_ratio[3][3]=0.00621093; fraction_ratio[3][4]=0.00502482; fraction_ratio[3][5]=0.00158842; fraction_ratio[3][6]=0.00188543; fraction_ratio[3][7]=0.00033812; fraction_ratio[3][8]=0.000521473; fraction_ratio[3][9]=0.000181353; fraction_ratio[3][10]=4.46707e-06; fraction_ratio[3][11]=1.48902e-06; fraction_ratio[3][12]=0; fraction_ratio[3][13]=0; fraction_ratio[3][14]=0; 
fraction_ratio[4][0]=-3.41404e-06; fraction_ratio[4][1]=3.41404e-06; fraction_ratio[4][2]=0.00193592; fraction_ratio[4][3]=0.000622117; fraction_ratio[4][4]=0.00179028; fraction_ratio[4][5]=0.000969449; fraction_ratio[4][6]=0.000274906; fraction_ratio[4][7]=2.32807e-05; fraction_ratio[4][8]=0; fraction_ratio[4][9]=0.000542496; fraction_ratio[4][10]=1.48902e-06; fraction_ratio[4][11]=0; fraction_ratio[4][12]=0; fraction_ratio[4][13]=0; fraction_ratio[4][14]=0; 
fraction_ratio[5][0]=0; fraction_ratio[5][1]=0; fraction_ratio[5][2]=0.00024199; fraction_ratio[5][3]=0.00024199; fraction_ratio[5][4]=0; fraction_ratio[5][5]=0; fraction_ratio[5][6]=0; fraction_ratio[5][7]=0; fraction_ratio[5][8]=0; fraction_ratio[5][9]=0; fraction_ratio[5][10]=0; fraction_ratio[5][11]=0; fraction_ratio[5][12]=0; fraction_ratio[5][13]=0; fraction_ratio[5][14]=0; 
fraction_ratio[6][0]=0; fraction_ratio[6][1]=0; fraction_ratio[6][2]=0; fraction_ratio[6][3]=0; fraction_ratio[6][4]=0; fraction_ratio[6][5]=3.41404e-06; fraction_ratio[6][6]=0; fraction_ratio[6][7]=0; fraction_ratio[6][8]=0; fraction_ratio[6][9]=0; fraction_ratio[6][10]=0; fraction_ratio[6][11]=0; fraction_ratio[6][12]=0; fraction_ratio[6][13]=0; fraction_ratio[6][14]=0; 










		// event loop : counting b-tagging for data
	  	for (Long64_t Nth = 0; Nth < Nev_data; Nth++)
	    	{
			//cout<<endl<<"=====================event:"<<Nth<<"==============================="<<endl;


			// draw hist before event selection
			
			int gen_b_this_ev=0;


			if(!is_pass_ev_selection(events_data[Nth])){continue;}
			N_available_event_data++;

			// draw rec Tmass
			double this_Tmass = get_Tmass(events_data[Nth]);
			if (this_Tmass!=-1){ h_rec_Tmass_a_evcut->Fill(this_Tmass);}

			
			//counting
			for(int i=0;i<events_data[Nth].THINjets.size();i++){ N_THINjet_before_deoverlap++;}
			for(int i=0;i<events_data[Nth].THINjets_deoverlap_with_lepton.size();i++){ N_THINjet_after_deoverlap++;}		


			vector<THINjet> this_THINjets=events_data[Nth].THINjets_deoverlap_with_lepton;

			// proccess for rec jet 
		
			// counter
			Int_t N_btagging=0;
		
			int N_THINjet_this_ev=0;
			for (int i=0;i<this_THINjets.size();i++)
			{
					THINjet this_j=this_THINjets[i];

					total_kk_rec++;
					N_THINjet_this_ev++;
					total_k_rec++;

					if ( btagging(this_j))
					{ 
						N_btagging++;
						total_btag_pass_cut_rec++;
					}
					else {total_nbtag_pass_cut_rec++;}
				
				
			}
			Nt_from_rec[N_btagging]++;
			Nk_from_rec[N_THINjet_this_ev]++;



	
		} 
	    	

		cout<<"N_available_event_data="<<N_available_event_data<<endl;
		cout<<"=========================  gen level  ==========================================="<<endl<<endl;
		cout<<"total_kk_gen="<<gen_j<<"(b="<<gen_b<<", nb="<<gen_nb<<"),"<<endl;
		cout<<"total_k_gen="<<gen_j_pass_cut<<"(b="<<gen_b_pass_cut<<", nb="<<gen_nb_pass_cut<<")"<<endl;
		cout<<"total_t="<<gen_j_pass_cut_btagging<<"(b="<<gen_b_pass_cut_btagging<<", nb="<<gen_nb_pass_cut_btagging<<")"<<endl;
		cout<<endl;

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
		cout<<"============================  fraction ratio  ========================================="<<endl<<endl;
		double sum_of_fraction_ratio=0;
		double N_of_jet_from_fraction_ratio=0;
		for (int i=0;i<fraction_ratio.size();i++)
		{
		
			for (int j=0;j<fraction_ratio[i].size();j++)
			{ 
				cout<<fraction_ratio[i][j]<<",    ";
				sum_of_fraction_ratio+=fraction_ratio[i][j];	
				N_of_jet_from_fraction_ratio+=(i+j)*fraction_ratio[i][j]*N_available_event_data;
			}
			cout<<endl;
		
		}
		cout<<"sum_of_fraction_ratio="<<sum_of_fraction_ratio<<endl;
		cout<<"N_of_jet_from_fraction_ratio="<<N_of_jet_from_fraction_ratio<<endl;
		cout<<endl;



		//get btagging 2p for rec
		// par_for_minuit(par_number,par_name, par_starting_value, step, min, max)
		double this_step =0.001;	
		par_for_minuit eff_rec(0,"eff_rec", 0.8, this_step ,0+this_step,1-this_step);
		par_for_minuit mr_rec(1,"mr_rec", 0.2, this_step ,0+this_step,1-this_step);	
		vector<par_for_minuit> pars_btagging2p_rec;
		pars_btagging2p_rec.push_back(eff_rec);
		pars_btagging2p_rec.push_back(mr_rec);
		vector<par_err_pair_for_minuit> results_btagging2p_rec=get_min_par_btagging2p_rec(pars_btagging2p_rec);
		best_eff_rec=results_btagging2p_rec[0].v;
		best_mr_rec=results_btagging2p_rec[1].v;



	
	
		cout<<"=========================  rec result  ==========================================="<<endl<<endl;
		cout<<",best_eff_rec="<<best_eff_rec<<"+-"<<results_btagging2p_rec[0].err<<endl;
		cout<<",best_mr_rec="<<best_mr_rec<<"+-"<<results_btagging2p_rec[1].err<<endl;
		cout<<endl;


		
		zprime_sample_result.add_result(Nev_data, N_available_event_data, true_eff, true_mr, best_eff_rec, results_btagging2p_rec[0].err, best_mr_rec,results_btagging2p_rec[1].err );
		
		
				
		
		
	//*/

	/*
	  	//event loop : analyzer
	  	for (Long64_t Nth = 0; Nth < Nev_data_pass_gen_selection; Nth++)
	    	{
			cout<<"=====================event:"<<Nth<<"==============================="<<endl;
		


	      		inner_data_for_data.GetEntry(Nth);

			// genpar var
			Int_t *genParId=inner_data_for_data.GetPtrInt("genParId");
			Int_t *genMomParId=inner_data_for_data.GetPtrInt("genMomParId");
			Int_t nGenPar=inner_data_for_data.GetInt("nGenPar");
			Int_t *genParSt=inner_data_for_data.GetPtrInt("genParSt");
			Int_t nGenJet=inner_data_for_data.GetInt("nGenJet");
			Int_t *genMo1=inner_data_for_data.GetPtrInt("genMo1");

		

	      		// THINjet var
	      		Int_t THINnJet=inner_data_for_data.GetInt("THINnJet");
	      		Float_t *THINjetPt=inner_data_for_data.GetPtrFloat("THINjetPt");  
	      		Float_t *THINjetEta=inner_data_for_data.GetPtrFloat("THINjetEta"); 
	      		Float_t *THINjetPhi=inner_data_for_data.GetPtrFloat("THINjetPhi"); 
	     		Float_t *THINjetEn=inner_data_for_data.GetPtrFloat("THINjetEn");
	     		Float_t *THINjetMass=inner_data_for_data.GetPtrFloat("THINjetMass");
	      		Float_t *THINjetCSV=inner_data_for_data.GetPtrFloat("THINjetCSV");

	  		Float_t *THINjetNHadEF=inner_data_for_data.GetPtrFloat("THINjetNHadEF");
	  		Float_t *THINjetNEmEF=inner_data_for_data.GetPtrFloat("THINjetNEmEF");
	  		Float_t *THINjetCHadEF=inner_data_for_data.GetPtrFloat("THINjetCHadEF");
	  		Float_t *THINjetCEmEF=inner_data_for_data.GetPtrFloat("THINjetCEmEF");
	  		Float_t *THINjetCMulti=inner_data_for_data.GetPtrFloat("THINjetCMulti");;                    
	      		// ELE var
	      		Int_t nEle =inner_data_for_data.GetInt("nEle");
	      		Float_t *elePt=inner_data_for_data.GetPtrFloat("elePt");
	   		Float_t *eleEta=inner_data_for_data.GetPtrFloat("eleEta");
	   		Float_t *elePhi=inner_data_for_data.GetPtrFloat("elePhi");
	   		Float_t *eleEnergy=inner_data_for_data.GetPtrFloat("eleEnergy");
	   		Float_t *eleM=inner_data_for_data.GetPtrFloat("eleM");

			Float_t *eleDelEtaIn=inner_data_for_data.GetPtrFloat("eleDelEtaIn");
			Float_t *eleDelPhiIn=inner_data_for_data.GetPtrFloat("eleDelPhiIn");
			Float_t *eleSigIhIh=inner_data_for_data.GetPtrFloat("eleSigIhIh");
			Float_t *eleHoE=inner_data_for_data.GetPtrFloat("eleHoE");
			Float_t *eleDxy=inner_data_for_data.GetPtrFloat("eleDxy");
			Float_t *eleDz=inner_data_for_data.GetPtrFloat("eleDz");
			Float_t *eleEoverP=inner_data_for_data.GetPtrFloat("eleEoverP");
			Float_t *eleCorrPfIso=inner_data_for_data.GetPtrFloat("eleCorrPfIso");
			Int_t *elePassConv=inner_data_for_data.GetPtrInt("elePassConv");
			Float_t *eleMissingHits=inner_data_for_data.GetPtrFloat("eleMissingHits");

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
	
	
		
	// save hist as png file		

	// insignificant var
	// ele var (insignificant var)
		
	save_hist( h_ele_pt_a_evcut, "pic_13Tev_data", "ele_pt_a_evcut");		
	save_hist( h_ele_eta_a_evcut, "pic_13Tev_data", "ele_eta_a_evcut"); 
	save_hist( h_ele_dEtaAtVtx_a_evcut, "pic_13Tev_data", "ele_dEtaAtVtx_a_evcut"); 
	save_hist( h_ele_dPhiAtVtx_a_evcut, "pic_13Tev_data", "ele_dPhiAtVtx_a_evcut"); 
	save_hist( h_ele_D0_a_evcut, "pic_13Tev_data", "ele_D0_a_evcut"); 
	save_hist( h_ele_EtaseedAtVtx_a_evcut, "pic_13Tev_data", "ele_EtaseedAtVtx_a_evcut"); 
	save_hist( h_ele_HoverE_a_evcut, "pic_13Tev_data", "ele_HoverE_a_evcut");
	save_hist( h_ele_MissHits_a_evcut, "pic_13Tev_data", "ele_MissHits_a_evcut");
	//save_hist( h_ele_MiniIso_a_evcut, "pic_13Tev_data", "ele_MiniIso_a_evcut");
	save_hist( h_ele_SigmaIEtaIEta_a_evcut, "pic_13Tev_data", "ele_SigmaIEtaIEta_a_evcut");
		
		//mu var (insignificant var)
		
	save_hist( h_mu_pt_a_evcut, "pic_13Tev_data", "mu_pt_a_evcut");	
	save_hist( h_mu_eta_a_evcut, "pic_13Tev_data", "mu_eta_a_evcut"); 
	save_hist( h_mu_dxy_a_evcut, "pic_13Tev_data", "mu_dxy_a_evcut"); 
	save_hist( h_mu_dz_a_evcut, "pic_13Tev_data", "mu_dz_a_evcut"); 
	save_hist( h_mu_Hits_a_evcut, "pic_13Tev_data", "mu_Hits_a_evcut"); 
	save_hist( h_mu_Matches_a_evcut, "pic_13Tev_data", "mu_Matches_a_evcut");
	//save_hist( h_mu_MiniIso_a_evcut, "pic_13Tev_data", "mu_MiniIso_a_evcut"); 
	save_hist( h_mu_PixelHits_a_evcut, "pic_13Tev_data", "mu_PixelHits_a_evcut"); 
	save_hist( h_mu_TrkLayers_a_evcut, "pic_13Tev_data", "mu_TrkLayers_a_evcut");
	save_hist( h_mu_TrkPt_a_evcut, "pic_13Tev_data", "mu_TrkPt_a_evcut"); 
		
	// THINjet var	
	
	save_hist(h_THINjet_CEmEF_a_evcut , "pic_13Tev_data" , "THINjet_CEmEF_a_evcut");
	save_hist(h_THINjet_CHadEF_a_evcut , "pic_13Tev_data" , "THINjet_CHadEF_a_evcut");
	save_hist(h_THINjet_NEmEF_a_evcut , "pic_13Tev_data" , "THINjet_NEmEF_a_evcut");
	save_hist(h_THINjet_NHadEF_a_evcut , "pic_13Tev_data" , "THINjet_NHadEF_a_evcut");
	save_hist(h_THINjet_PhoEF_a_evcut , "pic_13Tev_data" , "THINjet_PhoEF_a_evcut");	
	save_hist(h_THINjet_cisvv2_a_evcut , "pic_13Tev_data", "THINjet_cisvv2_a_evcut");
	save_hist(h_THINjet_pt_a_evcut, "pic_13Tev_data", "THINjet_pt_a_evcut");
	save_hist(h_THINjet_eta_a_evcut, "pic_13Tev_data", "THINjet_eta_a_evcut");

	save_hist(h_THINBjet_CEmEF_a_evcut , "pic_13Tev_data" , "THINBjet_CEmEF_a_evcut");
	save_hist(h_THINBjet_CHadEF_a_evcut , "pic_13Tev_data" , "THINBjet_CHadEF_a_evcut");
	save_hist(h_THINBjet_NEmEF_a_evcut , "pic_13Tev_data" , "THINBjet_NEmEF_a_evcut");
	save_hist(h_THINBjet_NHadEF_a_evcut , "pic_13Tev_data" , "THINBjet_NHadEF_a_evcut");
	save_hist(h_THINBjet_PhoEF_a_evcut , "pic_13Tev_data" , "THINBjet_PhoEF_a_evcut");
	save_hist(h_THINBtagjet_cisvv2_a_evcut , "pic_13Tev_data", "THINBtagjet_cisvv2_a_evcut");
	save_hist(h_THINBtagjet_pt_a_evcut, "pic_13Tev_data", "THINBtagjet_pt_a_evcut");
	save_hist(h_THINBtagjet_eta_a_evcut, "pic_13Tev_data", "THINBtagjet_eta_a_evcut");
	
	save_hist(h_rec_Tmass_a_evcut, "pic_13Tev_data", "rec_Tmass_a_evcut");


		
			
		//H_IM->Write();
		//H_IM_csvm->Write();
		//H_IM_csvl->Write();
		//Zee_IM->Write();
	/*
		THINCSV_b_evcut->Write(); 
		THINCSV_a_evcut->Write();
		
		leading_pruned_jet_mass_b_evcut->Write();
		leading_pruned_jet_mass_a_evcut->Write();
	*/	
	}
	//gerrors(zprime_sample_result, "eff");
	//gerrors(zprime_sample_result, "mr");
 	    
}
