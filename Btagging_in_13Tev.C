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

using namespace std;

// global variable================================================
  	//gereral
	Long64_t Nev_MC_before_gen_selection;
	Long64_t Nev_MC_pass_gen_selection =0 ;
	int N_available_event_MC=0;

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
	
	// tcone15 te3st
	float N_in_tcone15=0;
	float N_out_tcone15=0;



//class============================================================================
class PO
{
	public:
	TLorentzVector p4;
	PO(){}
	PO(TLorentzVector p4_)	
	{ p4=p4_;}
	

};
	

class gen_par : public PO
{
	public:
	Int_t id;
	Int_t index;
	Int_t mid;
	Int_t mo1;
	Int_t mo2;
	Int_t status;
	
	gen_par():PO(){}
    	gen_par(TLorentzVector p4_,
		Int_t id_,Int_t index_, Int_t mid_, Int_t mo1_, Int_t mo2_, Int_t status_)
		:PO(p4_)
	{ id=id_; index=index_; mid=mid_; mo1=mo1_; mo2=mo2_; status=status_;}
};

class ele : public PO
{
	public:
	/*
	Float_t del_eta_in;
	Float_t del_phi_in;
	Float_t sig_ieie;
	Float_t hoe;
	Float_t dxy;
	Float_t dz;
	Float_t eop;
	Float_t corr_pf_iso;
	Int_t pass_conv;
	Float_t missing_hits;

	Int_t pass_ID;
	Float_t user_trk_iso;
	Float_t user_cal_iso;
	Float_t rho;
	*/

	bool is_tight;
	bool is_medium;
	bool is_loose;

	ele():PO(){}
    	ele(TLorentzVector p4_ ,/*
			Float_t del_eta_in_, Float_t del_phi_in_, Float_t sig_ieie_, Float_t hoe_, 
			Float_t dxy_,Float_t dz_, Float_t eop_, Float_t corr_pf_iso_, Int_t pass_conv_, Float_t missing_hits_,
			Int_t pass_ID_, Float_t user_trk_iso_, Float_t user_cal_iso_, Float_t rho_ */
			bool is_tight_, bool is_medium_, bool is_loose_)
		:PO(p4_)
    	{
	/*	del_eta_in=del_eta_in_;
		del_phi_in=del_phi_in_;
		sig_ieie=sig_ieie_;
		hoe=hoe_;
		dxy=dxy_;
		dz=dz_;
		eop=eop_;
		corr_pf_iso=corr_pf_iso_;
		pass_conv=pass_conv_;
		missing_hits=missing_hits_;
		pass_ID=pass_ID_;
		user_trk_iso=user_trk_iso_;
		user_cal_iso=user_cal_iso_;
		rho=rho_;
	*/
		is_tight=is_tight_;
		is_medium=is_medium_;		
		is_loose=is_loose_;
	}

};

class mu : public PO
{
	public:

/*
    	Int_t is_global_mu;
    	Int_t is_tracker_mu;
    	Int_t trk_layers;
    	Int_t pixel_hits;
    	Int_t hits;
    	Int_t matches;
    	Float_t dxy;
    	Float_t dz;

	Float_t corr_pf_iso;
	Int_t pass_ID;
	Float_t corr_trk_iso;
*/	
	    		
	    		
	bool is_global;
	bool is_tracker;
	bool is_tight;
	bool is_medium;
	bool is_loose;
	bool is_high_pt;
	bool is_soft;
	
	mu():PO(){}
    	mu(TLorentzVector p4_ , /* //note: mu has no energy leafe
	   		Int_t is_global_mu_, Int_t is_tracker_mu_, Int_t trk_layers_, Int_t pixel_hits_, Int_t hits_, Int_t matches_, Float_t dxy_, Float_t dz_, 
			Float_t corr_pf_iso_, Int_t pass_ID_, Float_t corr_trk_iso_ */
			bool is_global_, bool is_tracker_, bool is_tight_, bool is_medium_, bool is_loose_, bool is_high_pt_, bool is_soft_)
		:PO(p4_)
    	{/*
    		is_global_mu=is_global_mu_;
    		is_tracker_mu=is_tracker_mu_;
    		trk_layers=trk_layers_;
    		pixel_hits=pixel_hits_;
    		hits=hits_;
    		matches=matches_;
    		dxy=dxy_;
    		dz=dz_;	
		corr_pf_iso=corr_pf_iso_;
		pass_ID=pass_ID_;
		corr_trk_iso=corr_trk_iso_;
	*/
		is_global=is_global_;
		is_tracker=is_tracker_;
		is_tight=is_tight_;
		is_medium=is_medium_;
		is_loose=is_loose_;
		is_high_pt=is_high_pt_;
		is_soft=is_soft_;
		
	}

};

  


class jet : public PO
{
  	public:
	Float_t csv;
	Float_t cisvv2;

	//cut var
	/*
	Float_t n_had_ef;
	Float_t n_em_ef;
	Float_t c_had_ef;
	Float_t c_em_ef;
	Float_t c_multi;
	*/

	jet():PO(){}	
    	jet(TLorentzVector p4_,
		Float_t csv_,Float_t cisvv2_/*,
		Float_t n_had_ef_,Float_t n_em_ef_,Float_t c_had_ef_,Float_t c_em_ef_,Float_t c_multi_*/)
		:PO(p4_)
    	{ 
		csv=csv_; cisvv2=cisvv2_;
	 	/*n_had_ef=n_had_ef_; n_em_ef=n_em_ef_; c_had_ef=c_had_ef_; c_em_ef=c_em_ef_; c_multi=c_multi_;*/
	}
};




class THINjet : public jet
{
  	public:
	//Int_t pass_id;
	
	bool is_tight;
	bool is_loose;

	
	THINjet():jet(){}	
    	THINjet(TLorentzVector p4_,
		Float_t csv_,Float_t cisvv2_,  /* Int_t pass_id_,*/  bool is_tight_ ,bool is_loose_ /*, 
		Float_t n_had_ef_,Float_t n_em_ef_,Float_t c_had_ef_,Float_t c_em_ef_,Float_t c_multi_*/)
		:jet(p4_,
		csv_ ,cisvv2_ /*, 
		n_had_ef_, n_em_ef_, c_had_ef_, c_em_ef_, c_multi_*/)
    	{  /*pass_id=pass_id_;*/ is_tight=is_tight_; is_loose=is_loose_; }
};



class met
{
	public:
	Float_t pt;
	Float_t phi;


	met(){}
	met( Float_t pt_, Float_t phi_)
	{ pt=pt_; phi=phi_;}
};



class event
{
	public:
	vector<gen_par> gen_pars;
	vector<THINjet> THINjets;
	vector<THINjet> THINjets_deoverlap_with_lepton;
	vector<THINjet> good_THINjets;
	vector<ele> eles;
	vector<ele> good_eles;
	vector<mu> mus;
	vector<mu> good_mus;
	met missing_et;
	event(){}
};




// counter

class jet_cone
{
	public:
	int gen_b;
	int gen_nb;
	int gen_ele;
	int gen_mu;
	jet_cone(){ gen_b=0;gen_nb=0;gen_ele=0;gen_mu=0;}
	

};

// for minuit
class par_for_minuit
{
	public:
	string name;
	double fit_start;
	double number;
	double fit_step;
	double min;
	double max;
	par_for_minuit(){}
	par_for_minuit(double number_, string name_, double fit_start_,  double fit_step_, double min_, double max_ )
	{  number=number_; name=name_; fit_start=fit_start_;  fit_step=fit_step_; min=min_; max=max_;}
};

class par_err_pair_for_minuit
{
	public:
	Double_t v;
	Double_t err;
	par_err_pair_for_minuit(){}
	par_err_pair_for_minuit(Double_t v_, Double_t err_)
	{ v=v_;err=err_;}
	
};

//sample result

class sample_result
{
	public:
	int N_sample;
	const char* name;
	vector<Long64_t> N_events;
	vector<Long64_t> N_available_events;
	vector<Float_t> gen_effs;
	vector<Float_t> gen_mrs;	
	vector<Float_t> rec_effs;
	vector<Float_t> rec_effs_err;
	vector<Float_t> rec_mrs;
	vector<Float_t> rec_mrs_err;
	
	sample_result(const char* name_){name=name_; N_sample=0;}
	void add_result(Long64_t N_event_, Long64_t N_available_event_,  Float_t gen_eff_,  Float_t gen_mr_, Float_t rec_eff_, Float_t rec_eff_err_, Float_t rec_mr_, Float_t rec_mr_err_)
	{ 
		N_sample++;
		N_events.push_back(N_event_);
		N_available_events.push_back(N_available_event_);
		gen_effs.push_back(gen_eff_);
		gen_mrs.push_back(gen_mr_);	
		rec_effs.push_back(rec_eff_);
		rec_effs_err.push_back(rec_eff_err_);
		rec_mrs.push_back(rec_mr_);
		rec_mrs_err.push_back(rec_mr_err_);
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

void save_hist(TH1D* h_, TString dirname, TString filename) // type: eff mr
{
	   
	TCanvas *c1 = new TCanvas(filename,filename,200,10,700,500);

	//c1->SetFillColor(42);
	c1->SetGrid();
	c1->GetFrame()->SetFillColor(21);
	c1->GetFrame()->SetBorderSize(12);

	h_->Draw();  // overlay
	
	TImage *img1 = TImage::Create();
	img1->FromPad(c1);
	img1->WriteImage("./"+dirname+"/"+filename+".png");
	delete c1;
	delete img1;

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

bool is_pass_csv_13Tev(jet j_)
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

bool is_pass_pteta(PO O1_)
{
	if(O1_.p4.Pt()<20){return false;}
	if(abs(O1_.p4.Eta())>2.4){return false;}
	return true;

}


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

bool is_pass_mu_tight_id(mu m_)
{
	if (m_.is_tight==true){return true;}
	else {return false;}

}

bool is_pass_mu_loose_id(mu m_)
{
	if (m_.is_loose==true){return true;}
	else {return false;}

}

//event selection
//evse gen level
// to check this event is TTbar semiletonic channel
bool is_ttbar_semileptonic_sample(event ev)
{
	int gen_q_from_W =0;
	int gen_l_from_W =0;
	for (int i=0;i<ev.gen_pars.size(); i++)
	{
		gen_par this_g = ev.gen_pars[i];
		if ( abs(this_g.id)<=4 && abs(this_g.id)>=1 && abs(this_g.mid)==24 )
		{gen_q_from_W++;}
		else if ((abs(this_g.id)==11 || abs(this_g.id)==13 )&& abs(this_g.mid)==24)
		{gen_l_from_W++;}
	}
	return (gen_q_from_W==2 && gen_l_from_W==1);

}


//evse rec level
bool is_pass_ev_selection_ttbar_semileptonic(event ev)
{
	int channel=0; //1 for ele, 2 for mu
	int ev_flag=0;

	// 1: event must have 1 good ele or 1 good muon
	vector<ele> good_eles=ev.good_eles;
	vector<mu> good_mus=ev.good_mus;
	
	//if (good_eles.size()<2 && good_mus.size()<2){ev_flag=1; return false;}
	if (good_eles.size()>=1 && good_mus.size()<1){channel=1;}
	else if (good_eles.size()<1 && good_mus.size()>=1){channel=2;}
	else {ev_flag=2; return false;}
	

	// 3: to check if have at least 2 good jet
	// using THINjet
	vector<THINjet> good_jets;
	
	
	for (int i=0;i<ev.THINjets_deoverlap_with_lepton.size();i++)
	{
	
		THINjet this_j = ev.THINjets_deoverlap_with_lepton[i];
		
		if(!(this_j.p4.Pt()>30 && abs(this_j.p4.Eta()<2.5)  &&  //is_pass_csv_13Tev(this_j) &&
		/*this_j.sd_mass>150 && this_j.sd_mass<200 &&  this_j.tau2/this_j.tau1<0.5 && */ is_pass_THINjet_loose_id(this_j)))
		{continue;}
						
		good_jets.push_back(this_j);

	}
	if (good_jets.size()<2){ ev_flag=3; return false;}

	// 4: to check if have 2 good leading jet
	// using THINjet
	bool have_good_leading_jets=false;
	for (int i=0;i<good_jets.size();i++)
	{
		for (int j=0;j<good_jets.size();j++)
		{
			if(i==j){continue;}
			if(good_jets[i].p4.Pt()>=70 && good_jets[j].p4.Pt()>=50){have_good_leading_jets=true;}
		}
	}
	if(!have_good_leading_jets){ev_flag=4;return false;}
	
	
	// 5: to check if missing Et > 20
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





// calculate likelihood



// calculate expection value of Nt
double get_expvalue_of_Nt(  double eff, double mr, int t)
{
        double this_expvalue = 0;
        for (int i = 0; i < fraction_ratio.size() ; i++)
        {
                for (int j = 0; j < fraction_ratio[i].size() ; j++)
                {
                    	double term1 = N_available_event_MC * fraction_ratio[i][j];

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
void Btagging_in_13Tev()
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
	ev_se is_pass_ev_selection = is_pass_ev_selection_ttbar_semileptonic;	
	// Btagging function
	typedef bool (*jet_alg)(jet);
	//jet_alg btagging = is_pass_csvl;
	//jet_alg btagging = is_pass_csvm;
	jet_alg btagging = is_pass_csv_13Tev;
	// is wanted b
	typedef bool (*chk_gen_b)(gen_par);
	chk_gen_b is_wanted_b = is_wanted_b_ttbar;
	//chk_gen_b is_wanted_b = is_wanted_b_Zprime;

		
		
	
	
	sample_result zprime_sample_result("zprime_sample_result");
	
	for (int sample_i=0;sample_i<1; sample_i++) //collect all small samples in a big sample
	{
		// initialization variables
		N_available_event_MC=0;

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
		TH1D* h_THINcisvv2_b_evcut = new TH1D("THINcisvv2_b_evcut","THINcisvv2_b_evcut",100,0,1);
		TH1D* h_THINcisvv2_a_evcut = new TH1D("THINcisvv2_a_evcut","THINcisvv2_a_evcut",100,0,1);

		TH1D* h_genbw_deltaR_b_evcut = new TH1D("genbw_deltaR_b_evcut","genbw_deltaR_b_evcut",100,0,3);		
		TH1D* h_THINnjet = new TH1D("THINnjet","THINnjet",20,0,20);
		TH1D* h_THINnjet_deoverlap = new TH1D("THINnjet_deoverlap","THINnjet_deoverlap",20,0,20);
		TH1D* h_N_THIN_b_jet = new TH1D("N_THIN_b_jet","N_THIN_b_jet",5,0,5);
		TH1D* h_Ngenbjet = new TH1D("h_Ngenbjet","h_Ngenbjet",5,0,5);
		
		
		// gen level hist set
		TH1D* h_dr_tb = new TH1D("h_dr_tb","h_dr_tb",700,0,7);
		TH1D* h_dr_tq = new TH1D("h_dr_tq","h_dr_tq",700,0,7);
		TH1D* h_dr_tqbar = new TH1D("h_dr_tqbar","h_dr_tqbar",700,0,7);
		
		TH1D* h_dr_tP4b = new TH1D("h_dr_tP4b","h_dr_tP4b",700,0,7);
		TH1D* h_dr_tP4q = new TH1D("h_dr_tP4q","h_dr_tP4q",700,0,7);
		TH1D* h_dr_tP4qbar = new TH1D("h_dr_tP4qbar","h_dr_tP4qbar",700,0,7);
		
		
	
	
		// root file
	  	//TChain *MC = new TChain("tree/treeMaker");
	  	//MC->Add("../flattuple.root");// Zprime sample
	  	//MC->Add("../TTbar_sample/flattuple_1.root");// TTbar sample from yu
	  	//MC->Add(samples[sample_i]);
	  	
	  	vector<string> input_files;
	  	for (int i=1;i<=447;i++)
	  	{
	  		if (i%5!=0){continue;};
	  		if (i==78 || i==161 || i==184 || i==424){continue;}
	  		input_files.push_back("../13TevTTbar/NCUGlobalTuples_"+int_to_string(i)+".root");
	  	}
	  	
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_100.root");
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_200.root");
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_300.root");
	  	//input_files.push_back("../13TevTTbar/NCUGlobalTuples_400.root");
	  	
	  				
	  	//TreeReader inner_data_for_MC(MC);
	  	TreeReader inner_data_for_MC(input_files);
	  	Nev_MC_before_gen_selection = inner_data_for_MC.GetEntriesFast();


	  	cout<< "Nev_MC_before_gen_selection : "<< Nev_MC_before_gen_selection << endl;
	
	  	
	  	
	  	
		//configure
		vector<event> events_MC;

	
		int n_THINj=0;
	



		//matching

		// event loop : set up MC
	  	//for (Long64_t Nth = 0; Nth < Nev_MC_before_gen_selection; Nth++)
	  	for (Long64_t Nth = 0; Nth < Nev_MC_before_gen_selection; Nth++)
	    	{
	 		inner_data_for_MC.GetEntry(Nth);

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
			
			
			vector<bool> &THINjetPassIDTight = *((vector<bool>*) inner_data_for_MC.GetPtr("THINjetPassIDTight"));			
			vector<bool> &THINjetPassIDLoose  = *((vector<bool>*) inner_data_for_MC.GetPtr("THINjetPassIDLoose"));

	      	      		
/*      	      	
	  		Float_t *THINjetNHadEF=inner_data_for_MC.GetPtrFloat("THINjetNHadEF");
	  		Float_t *THINjetNEmEF=inner_data_for_MC.GetPtrFloat("THINjetNEmEF");
	  		Float_t *THINjetCHadEF=inner_data_for_MC.GetPtrFloat("THINjetCHadEF");
	  		Float_t *THINjetCEmEF=inner_data_for_MC.GetPtrFloat("THINjetCEmEF");
	  		Float_t *THINjetCMulti=inner_data_for_MC.GetPtrFloat("THINjetCMulti");
*/

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


			event this_ev;

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
	    		
	    			mu mu_i( *((TLorentzVector*)muP4->At(i)), isGlobalMuon[i], isTrackerMuon[i], isTightMuon[i], isMediumMuon[i], isLooseMuon[i],isHighPtMuon[i], isSoftMuon[i]);
				this_ev.mus.push_back(mu_i);
			}
		
	      		for (Int_t i=0;i<THINnJet;i++)
			{	
				//if (THINjetCSV[i]<0 || THINjetCSV[i]>1){continue;}
				if (THINjetCISVV2[i]<0 || THINjetCISVV2[i]>1){continue;}

								
				
				THINjet jet_i( *((TLorentzVector*)THINjetP4->At(i)), THINjetCSV[i], THINjetCISVV2[i],  THINjetPassIDTight[i],THINjetPassIDLoose[i]);
				this_ev.THINjets.push_back(jet_i);

			}
			for(int i=0;i<nGenPar;i++)  //gen par loop for TTbar sample
			{
				// only include qqbar ttbar bbar ele mu 
				if 
				( 
					(abs(genParId[i])==6)  ||
					(abs(genParId[i])==5  && abs(genMomParId[i])==6) ||
					( abs(genParId[i])<=5  && abs(genParId[i])>=1  && abs(genMomParId[i])==24) ||
					( abs(genParId[i])==11  && abs(genMomParId[i])==24) || 
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
				if ( !is_pass_mu_loose_id(p) ){continue;}
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

				for (int mui=0;mui<this_ev.good_mus.size();mui++)
				{
					mu this_mu = this_ev.good_mus[mui];
					if (delta_R(this_j, this_mu)<0.4){overlap_with_lepton=true;break;}
				}
				if (overlap_with_lepton){ continue;}
			
				this_ev.THINjets_deoverlap_with_lepton.push_back(this_j);
				
			}
			

			// to draw hist in all channel and before event selection
			
			int gen_b_index=-1;
			int gen_bbar_index=-1;
			int gen_t_index=-1;
			int gen_tbar_index=-1;
			int gen_q_from_wplus_index=-1;
			int gen_qbar_from_wplus_index=-1;
			int gen_q_from_wminus_index=-1;
			int gen_qbar_from_wminus_index=-1;
			for (int i=0;i< this_ev.gen_pars.size();i++)
			{
				gen_par this_g = this_ev.gen_pars[i];
				if (this_g.id==6){gen_t_index=i;}      // t
				else if (this_g.id==-6){gen_tbar_index=i;}	// tbar
				else if (this_g.id==5 && this_g.mid == 6){gen_b_index=i;}	//b from t
				else if (this_g.id==-5 && this_g.mid == -6){gen_bbar_index=i;} //bbar from tbar
				else if (this_g.id>=0 && this_g.id<=5 && this_g.mid == 24){gen_q_from_wplus_index=i;} //q from W+	
				else if (this_g.id<=0 && this_g.id>=-5 && this_g.mid == 24){gen_qbar_from_wplus_index=i;} //qbar from W+
				else if (this_g.id>=0 && this_g.id<=5 && this_g.mid == -24){gen_q_from_wminus_index=i;}  //q from W-
				else if (this_g.id<=0 && this_g.id>=-5 && this_g.mid == -24){gen_qbar_from_wminus_index=i;} //qbar from W-				
			}
		
		
			// h_dr_tb  :  delta R(t or tbar, daughter b or bbar )
			if ( gen_t_index!=-1 && gen_b_index!=-1 )
			{h_dr_tb->Fill(delta_R(this_ev.gen_pars[gen_t_index],this_ev.gen_pars[gen_b_index])); }
			if ( gen_tbar_index!=-1 && gen_bbar_index!=-1)
			{h_dr_tb->Fill(delta_R(this_ev.gen_pars[gen_tbar_index],this_ev.gen_pars[gen_bbar_index])); }
			// h_dr_tq : delta R(t or tbar, daughter q )
			if ( gen_t_index!=-1 && gen_q_from_wplus_index!=-1 )
			{h_dr_tq->Fill(delta_R(this_ev.gen_pars[gen_t_index],this_ev.gen_pars[gen_q_from_wplus_index])); }
			if ( gen_tbar_index!=-1 && gen_q_from_wminus_index!=-1)
			{h_dr_tq->Fill(delta_R(this_ev.gen_pars[gen_tbar_index],this_ev.gen_pars[gen_q_from_wminus_index])); }
			// h_dr_tqbar : delta R(t or tbar, daughter qbar )
			if ( gen_t_index!=-1 && gen_qbar_from_wplus_index!=-1 )
			{h_dr_tqbar->Fill(delta_R(this_ev.gen_pars[gen_t_index],this_ev.gen_pars[gen_qbar_from_wplus_index])); }
			if ( gen_tbar_index!=-1 && gen_qbar_from_wminus_index!=-1)
			{h_dr_tqbar->Fill(delta_R(this_ev.gen_pars[gen_tbar_index],this_ev.gen_pars[gen_qbar_from_wminus_index])); }
			

			// h_dr_tP4b, h_dr_tP4q,  h_dr_tP4qbar
			if (gen_t_index!=-1 && gen_b_index!=-1 && gen_q_from_wplus_index!=-1 && gen_qbar_from_wplus_index!=-1)
			{
				TLorentzVector bl =this_ev.gen_pars[gen_b_index].p4;
				TLorentzVector ql =this_ev.gen_pars[gen_q_from_wplus_index].p4;
				TLorentzVector qbarl =this_ev.gen_pars[gen_qbar_from_wplus_index].p4;
				TLorentzVector tl = bl+ql+qbarl;
				h_dr_tP4b->Fill(delta_R(tl, bl ));
				h_dr_tP4q->Fill(delta_R(tl, ql ));
				h_dr_tP4qbar->Fill(delta_R(tl, qbarl ));
				if (delta_R(tl, bl )<=1.5 && delta_R(tl, ql )<1.5 && delta_R(tl, qbarl )<1.5)
				{ N_in_tcone15++;}
				else 
				{N_out_tcone15++;}
			}
			if (gen_tbar_index!=-1 && gen_bbar_index!=-1 && gen_q_from_wminus_index!=-1 && gen_qbar_from_wminus_index!=-1)
			{
				TLorentzVector bl =this_ev.gen_pars[gen_bbar_index].p4;
				TLorentzVector ql =this_ev.gen_pars[gen_q_from_wminus_index].p4;
				TLorentzVector qbarl =this_ev.gen_pars[gen_qbar_from_wminus_index].p4;
				TLorentzVector tl = bl+ql+qbarl;
				h_dr_tP4b->Fill(delta_R(tl, bl ));
				h_dr_tP4q->Fill(delta_R(tl, ql ));
				h_dr_tP4qbar->Fill(delta_R(tl, qbarl ));
				if (delta_R(tl, bl )<=1.5 && delta_R(tl, ql )<1.5 && delta_R(tl, qbarl )<1.5)
				{ N_in_tcone15++;}
				else 
				{N_out_tcone15++;}
			}
			
			
			

			// to check this ttbar sample is semiletonic channel
			if (!is_ttbar_semileptonic_sample(this_ev)){continue;}
			Nev_MC_pass_gen_selection++;
			events_MC.push_back(this_ev);

		}
		cout<< "Nev_MC_pass_gen_selection : "<< Nev_MC_pass_gen_selection << endl;
	
		// event loop : get size of fraction ratio matrix
		int max_N_of_candidate_bjet=0;
		int max_N_of_candidate_nbjet=0;
	  	for (Long64_t Nth = 0; Nth < Nev_MC_pass_gen_selection; Nth++)
	    	{
			//cout<<endl<<"=====================event:"<<Nth<<"==============================="<<endl;
	
			// matching gen par and THINjet

			if(!is_pass_ev_selection(events_MC[Nth])){continue;}

	 		//to check this event has 2 b-THINjet, to check deltaR of this 2 subjet > 0.3
			int N_of_THINjet=0;
			int N_of_candidate_bjet=0;
			int N_of_candidate_nbjet=0;
			double R=0.4;

			vector<THINjet>  this_THINjets=events_MC[Nth].THINjets_deoverlap_with_lepton;
			for (int i=0;i<this_THINjets.size();i++)
			{	
			
					THINjet this_j=this_THINjets[i];

					if(!is_pass_pteta(this_j)){continue;}

					N_of_THINjet++;
					bool matching_with_b=false;
					for(int ii=0;ii<events_MC[Nth].gen_pars.size();ii++)  //gen par loop
					{
						gen_par this_g = events_MC[Nth].gen_pars[ii];
						if (delta_R(this_j,this_g)<=R  && is_wanted_b(this_g))
						{
							matching_with_b=true; 
						}
						
					}
					if (matching_with_b){ N_of_candidate_bjet++;}

				
			}

			// draw hist h_N_THIN_b_jet
			h_N_THIN_b_jet->Fill(N_of_candidate_bjet);
			
			N_available_event_MC++;
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
		Nk_from_rec = vector<Int_t>(50,0);
		Nt_from_rec = vector<Int_t>(50,0);
		
	
		cout<<"max_true_b="<<max_true_b<<", max_true_nb="<<max_true_nb<<";    ";
		vector< vector<double> > mask(max_N_of_candidate_bjet+1, vector<double>(max_N_of_candidate_nbjet+1,0));
		fraction_ratio=mask;

		// event loop : get fraction ratio matrix
	  	for (Long64_t Nth = 0; Nth < Nev_MC_pass_gen_selection; Nth++)
	    	{
			//cout<<endl<<"=====================event:"<<Nth<<"==============================="<<endl;
	
			if(!is_pass_ev_selection(events_MC[Nth])){continue;}

			int N_of_THINjet=0;
			int N_of_candidate_bjet=0;
			int N_of_candidate_nbjet=0;

			double R=0.4;

			vector<THINjet>  this_THINjets=events_MC[Nth].THINjets_deoverlap_with_lepton;
			for (int i=0;i<this_THINjets.size();i++)
			{	

					THINjet this_j=this_THINjets[i];

					if(!is_pass_pteta(this_j)){continue;}

					N_of_THINjet++;
					bool matching_with_b=false;
					for(int ii=0;ii<events_MC[Nth].gen_pars.size();ii++)  //gen par loop
					{
						gen_par this_g = events_MC[Nth].gen_pars[ii];
						if (delta_R(this_j,this_g)<=R  && is_wanted_b(this_g))
						{matching_with_b=true;}
						
					}
					if (matching_with_b){ N_of_candidate_bjet++;}

				
			}

			N_of_candidate_nbjet=N_of_THINjet-N_of_candidate_bjet;
			fraction_ratio[N_of_candidate_bjet][N_of_candidate_nbjet]+=(1/double(N_available_event_MC));
		}	






		// event loop : counting b-tagging for MC
	  	for (Long64_t Nth = 0; Nth < Nev_MC_pass_gen_selection; Nth++)
	    	{
			//cout<<endl<<"=====================event:"<<Nth<<"==============================="<<endl;
		




			// draw hist before event selection
			
			h_THINnjet->Fill(events_MC[Nth].THINjets.size());
			h_THINnjet_deoverlap->Fill(events_MC[Nth].THINjets_deoverlap_with_lepton.size()); //note: dont include subjet deoverlap 
			for(int i=0;i<events_MC[Nth].THINjets_deoverlap_with_lepton.size();i++)
			{
				h_THINcisvv2_b_evcut->Fill(events_MC[Nth].THINjets_deoverlap_with_lepton[i].cisvv2);
			}
			
			int gen_b_this_ev=0;
			for(int i=0;i<events_MC[Nth].gen_pars.size();i++)
			{
				if (is_wanted_b(events_MC[Nth].gen_pars[i]))
				{ gen_b_this_ev++;}
			}
			h_Ngenbjet->Fill(gen_b_this_ev);
			
			

	
		
			// draw genbw_deltaR		
			int gen_b_index=-1;
			int gen_wplus_index=-1;
			int gen_bbar_index=-1;
			int gen_wminus_index=-1;
			for (int i=0;i<events_MC[Nth].gen_pars.size();i++)
			{
				gen_par this_g = events_MC[Nth].gen_pars[i];
				if(this_g.id==5 && this_g.mid==6){gen_b_index=i; break;}
			}
			for (int i=0;i<events_MC[Nth].gen_pars.size();i++)
			{
				gen_par this_g = events_MC[Nth].gen_pars[i];
				if(this_g.id==24 && this_g.mid==6){gen_wplus_index=i; break;}
			}
			for (int i=0;i<events_MC[Nth].gen_pars.size();i++)
			{
				gen_par this_g = events_MC[Nth].gen_pars[i];
				if(this_g.id==-5 && this_g.mid==-6){gen_bbar_index=i; break;}
			}
			for (int i=0;i<events_MC[Nth].gen_pars.size();i++)
			{
				gen_par this_g = events_MC[Nth].gen_pars[i];
				if(this_g.id==-24 && this_g.mid==-6){gen_wminus_index=i; break;}
			}			

			
			if (gen_b_index!=-1 && gen_wplus_index!=-1 && gen_bbar_index!=-1 && gen_wminus_index!=-1)
			{
				h_genbw_deltaR_b_evcut->Fill(delta_R(events_MC[Nth].gen_pars[gen_b_index], events_MC[Nth].gen_pars[gen_wplus_index]) );
				h_genbw_deltaR_b_evcut->Fill(delta_R(events_MC[Nth].gen_pars[gen_bbar_index], events_MC[Nth].gen_pars[gen_wminus_index]) );
			}

			
			
		

			if(!is_pass_ev_selection(events_MC[Nth])){continue;}

			// draw hist in event passed event selection
			for(int i=0;i<events_MC[Nth].THINjets_deoverlap_with_lepton.size();i++)
			{
				h_THINcisvv2_a_evcut->Fill(events_MC[Nth].THINjets_deoverlap_with_lepton[i].cisvv2);
			}

			
			//counting
			for(int i=0;i<events_MC[Nth].THINjets.size();i++){ N_THINjet_before_deoverlap++;}
			for(int i=0;i<events_MC[Nth].THINjets_deoverlap_with_lepton.size();i++){ N_THINjet_after_deoverlap++;}		


			vector<THINjet> this_THINjets=events_MC[Nth].THINjets_deoverlap_with_lepton;
		
		
	 		//calculate gen jets
			for (int i=0;i<this_THINjets.size();i++)
			{

					THINjet this_j=this_THINjets[i];
	
					if(!is_pass_pteta(this_j)){continue;}
		
					n_THINj++;
					string this_j_is="";
							
					// to see if this sj is b or nb
			
					double R=0.4;
					jet_cone this_jet_cone;
			
					for(int ii=0;ii<events_MC[Nth].gen_pars.size();ii++)  //gen par loop
					{
						gen_par this_g = events_MC[Nth].gen_pars[ii];
						if (delta_R(this_j,this_g)<=R)
						{
							int this_g_id=abs(this_g.id);
							if(is_wanted_b(this_g)){this_jet_cone.gen_b++;}
							if(this_g_id<=4 && this_g_id>=0){this_jet_cone.gen_nb++;}
							if(this_g_id==11){this_jet_cone.gen_ele++;}
							if(this_g_id==13){this_jet_cone.gen_mu++;}				
						}	
					}
					if (this_jet_cone.gen_b>0)
					{
						gen_j++; gen_b++;
						if(this_jet_cone.gen_b==1){this_j_is="b";}
						else if (this_jet_cone.gen_b==2){this_j_is="2b";}				
					}
					else{this_j_is="nb"; gen_j++; gen_nb++;} //actually, this_j matching with nbq, ele, or mu, or nothing. 


			
					//get 2p
					if (this_jet_cone.gen_b>0)
					{
						if ( is_pass_pteta(this_j))
						{
							gen_j_pass_cut++;
							gen_b_pass_cut++;
				
							if(btagging(this_j))
							{ gen_j_pass_cut_btagging++; gen_b_pass_cut_btagging++; }
						}
					}
					else if (this_j_is=="nb")
					{
						if ( is_pass_pteta(this_j))
						{
							gen_j_pass_cut++;
							gen_nb_pass_cut++;
				
							if(btagging(this_j))
							{ gen_j_pass_cut_btagging++; gen_nb_pass_cut_btagging++; }
						}
					}

				
			}

			// proccess for rec jet 
		
			// counter
			Int_t N_btagging=0;
			Int_t N_jet_pass_cut=0;
		
			int N_THINjet_this_ev=0;
			for (int i=0;i<this_THINjets.size();i++)
			{


					THINjet this_j=this_THINjets[i];

					if(!is_pass_pteta(this_j)){continue;}

					total_kk_rec++;
					N_THINjet_this_ev++;
					if ( is_pass_pteta(this_j))
					{N_jet_pass_cut++;total_k_rec++;}
					else
					{continue;}
					if ( btagging(this_j))
					{ N_btagging++;total_btag_pass_cut_rec++;}
					else {total_nbtag_pass_cut_rec++;}
				
				
			}
			Nt_from_rec[N_btagging]++;
			Nk_from_rec[N_jet_pass_cut]++;
		
			//----------------------------k----------------------t-------------------------------
			//cout<<"this event have "<<N_THINsubjet_this_ev<<"jet, and "<<N_jet_pass_cut<< "jet pass cut, and "<<N_btagging<<" btagging."<<endl;


	
		} 
	    	

		cout<<"N_available_event_MC="<<N_available_event_MC<<endl;
		cout<<"N_available_event_MC/Nev_MC_pass_gen_selection ="<<double(N_available_event_MC)/double(Nev_MC_pass_gen_selection)<<endl;
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

		cout<<"=========================  2p in gen level  ==========================================="<<endl<<endl;
		true_eff=(double)gen_b_pass_cut_btagging/(double)gen_b_pass_cut;
		true_mr=(double)gen_nb_pass_cut_btagging/(double)gen_nb_pass_cut;
		cout<<"eff="<<true_eff<<endl;
		cout<<"mr="<<true_mr<<endl;
		cout<<endl;

/*
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
	
		cout<<"=========================  Nk_from_data  ==========================================="<<endl<<endl;
		double sum_of_Nk_from_data=0;
		double N_of_jet_from_Nk_from_data=0;
		for (int i =0;i<Nk_from_data.size();i++)
		{ 
			cout<<"Nk_from_data["<<i<<"]:"<<Nk_from_data[i]<<endl;
			sum_of_Nk_from_data+=Nk_from_data[i];
			N_of_jet_from_Nk_from_data+=Nk_from_data[i]*i;
		}
		cout<<"sum_of_Nk_from_data="<<sum_of_Nk_from_data<<endl;
		cout<<"N_of_jet_from_Nk_from_data="<<N_of_jet_from_Nk_from_data<<endl;
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
	
		cout<<"=========================  Nt_from_data  ==========================================="<<endl<<endl;
		double sum_of_Nt_from_data=0;
		double N_of_jet_from_Nt_from_data=0;
		for (int i =0;i<Nt_from_data.size();i++)
		{ 
			cout<<"Nt_from_data["<<i<<"]:"<<Nt_from_data[i]<<endl;
			sum_of_Nt_from_data+=Nt_from_data[i];
			N_of_jet_from_Nt_from_data+=Nt_from_data[i]*i;
		}
		cout<<"sum_of_Nt_from_data="<<sum_of_Nt_from_data<<endl;
		cout<<"N_of_jet_from_Nt_from_data="<<N_of_jet_from_Nt_from_data<<endl;
		cout<<endl;

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

		cout<<"============================  fraction ratio  ========================================="<<endl<<endl;
		double sum_of_fraction_ratio=0;
		double N_of_jet_from_fraction_ratio=0;
		for (int i=0;i<fraction_ratio.size();i++)
		{
		
			for (int j=0;j<fraction_ratio[i].size();j++)
			{ 
				cout<<fraction_ratio[i][j]<<",    ";
				sum_of_fraction_ratio+=fraction_ratio[i][j];	
				N_of_jet_from_fraction_ratio+=(i+j)*fraction_ratio[i][j]*N_available_event_MC;
			}
			cout<<endl;
		
		}
		cout<<"sum_of_fraction_ratio="<<sum_of_fraction_ratio<<endl;
		cout<<"N_of_jet_from_fraction_ratio="<<N_of_jet_from_fraction_ratio<<endl;
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



	
	
		cout<<"=========================  rec result  ==========================================="<<endl<<endl;
		cout<<",best_eff_rec="<<best_eff_rec<<"+-"<<results_btagging2p_rec[0].err<<endl;
		cout<<",best_mr_rec="<<best_mr_rec<<"+-"<<results_btagging2p_rec[1].err<<endl;
		cout<<endl;

		cout<<"=========================  ratio rec/gen  ==========================================="<<endl<<endl;
		cout<<",best_eff_rec/true_eff="<<best_eff_rec/true_eff<<endl;
		cout<<",best_mr_rec/true_mr="<<best_mr_rec/true_mr<<endl;
		cout<<endl;
		
		cout<<"=========================  tcone15 test  ==========================================="<<endl<<endl;
		cout<<"N_in_tcone15/(N_in_tcone15+N_out_tcone15)="<<N_in_tcone15/(N_in_tcone15+N_out_tcone15)<<endl;
		cout<<endl;
				
		zprime_sample_result.add_result(Nev_MC_pass_gen_selection, N_available_event_MC, true_eff, true_mr, best_eff_rec, results_btagging2p_rec[0].err, best_mr_rec,results_btagging2p_rec[1].err );
		
		
		
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
	
	
		
		
	save_hist(h_THINcisvv2_b_evcut, "pic_13Tev_TTbar", "THINcisvv2_b_evcut");
	save_hist(h_THINcisvv2_a_evcut, "pic_13Tev_TTbar", "THINcisvv2_a_evcut");
	save_hist(h_genbw_deltaR_b_evcut, "pic_13Tev_TTbar", "genbw_deltaR_b_evcut");
	save_hist(h_THINnjet, "pic_13Tev_TTbar", "THINnjet");
	save_hist(h_THINnjet_deoverlap, "pic_13Tev_TTbar", "THINnjet_deoverlap");
	save_hist(h_N_THIN_b_jet, "pic_13Tev_TTbar", "N_THIN_b_jet");
	save_hist(h_Ngenbjet, "pic_13Tev_TTbar", "Ngenbjet");

	save_hist(h_dr_tb, "pic_13Tev_TTbar", "h_dr_tb");
	save_hist(h_dr_tq, "pic_13Tev_TTbar", "h_dr_tq");
	save_hist(h_dr_tqbar, "pic_13Tev_TTbar", "h_dr_tqbar");
	save_hist(h_dr_tP4b, "pic_13Tev_TTbar", "h_dr_tP4b");
	save_hist(h_dr_tP4q, "pic_13Tev_TTbar", "h_dr_tP4q");
	save_hist(h_dr_tP4qbar, "pic_13Tev_TTbar", "h_dr_tP4qbar");
		
			
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
