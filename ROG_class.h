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
	
	Float_t ChHadIso;
	Float_t NeHadIso;
	Float_t GamIso;
	Float_t PUPt;
	
	mu():PO(){}
    	mu(TLorentzVector p4_ , /* //note: mu has no energy leafe
	   		Int_t is_global_mu_, Int_t is_tracker_mu_, Int_t trk_layers_, Int_t pixel_hits_, Int_t hits_, Int_t matches_, Float_t dxy_, Float_t dz_, 
			Float_t corr_pf_iso_, Int_t pass_ID_, Float_t corr_trk_iso_ */
			bool is_global_, bool is_tracker_, bool is_tight_, bool is_medium_, bool is_loose_, bool is_high_pt_, bool is_soft_,
			Float_t ChHadIso_, Float_t NeHadIso_, Float_t GamIso_, Float_t PUPt_)
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
	
		ChHadIso=ChHadIso_;
		NeHadIso=NeHadIso_;
		GamIso=GamIso_;
		PUPt=PUPt_;
		
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
	
	Int_t hadron_flavor;
	bool is_tight;
	bool is_loose;

	Float_t CEmEF;
	Float_t CHadEF;
	Float_t NEmEF;
	Float_t NHadEF;
	Float_t PhoEF;
	
	THINjet():jet(){}	
    	THINjet(TLorentzVector p4_,
		Float_t csv_,Float_t cisvv2_, Int_t hadron_flavor_, /* Int_t pass_id_,*/  bool is_tight_ ,bool is_loose_ , 
		Float_t CEmEF_,Float_t CHadEF_,Float_t NEmEF_,Float_t NHadEF_,Float_t PhoEF_)
		:jet(p4_,
		csv_ ,cisvv2_ /*, 
		n_had_ef_, n_em_ef_, c_had_ef_, c_em_ef_, c_multi_*/)
    	{  	
    		hadron_flavor=hadron_flavor_; is_tight=is_tight_; is_loose=is_loose_;
    		CEmEF=CEmEF_; CHadEF=CHadEF_; NEmEF=NEmEF_; NHadEF=NHadEF_; PhoEF=PhoEF_;	
    	}
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
	bool is_pass_ele_trigger = false;
	bool is_pass_mu_trigger = false;
	Float_t event_weight;
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




