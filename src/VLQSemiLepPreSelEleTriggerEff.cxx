#include <algorithm>
#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/PartonHT.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/CollectionProducer.h"
#include "UHH2/common/include/TriggerSelection.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/common/include/NSelections.h"

#include "UHH2/VLQSemiLepPreSel/include/VLQCommonModules.h"
#include "UHH2/VLQSemiLepPreSel/include/VLQSemiLepPreSelHists.h"
#include "UHH2/VLQSemiLepPreSel/include/VLQGenHists.h"
#include "UHH2/VLQSemiLepPreSel/include/EventHists.h"
#include "UHH2/VLQSemiLepPreSel/include/SelectionHists.h"
#include "UHH2/VLQSemiLepPreSel/include/VLQTrigger_selectionItems.h"




using namespace std;
using namespace uhh2;


template <typename TYPE>
static bool is_true(const TYPE &, const Event &) {
    return true;
}





class VLQSemiLepPreSelEleTriggerEff: public AnalysisModule {
public:
    explicit VLQSemiLepPreSelEleTriggerEff(Context & ctx);
    virtual bool process(Event & event) override;

private:
    std::string version;
    std::string type;
    // modules for setting up collections and cleaning
    std::vector<std::unique_ptr<AnalysisModule>> v_pre_modules;
    unique_ptr<SelectionProducer> sel_module;
    unique_ptr<AnalysisModule> common_modules_with_lumi_sel;
    std::unique_ptr<Selection> trigger_sel_Ele50;
    std::unique_ptr<Selection> trigger_sel_Photo;
    std::unique_ptr<Selection> trigger_sel_Mu50;
    std::unique_ptr<Selection> trigger_sel_TkMu50;
    std::unique_ptr<Selection> mu_sel;
    std::unique_ptr<Selection> ele_sel;

    unique_ptr<AnalysisModule> Mu_SF_BtoF;
    unique_ptr<AnalysisModule> Mu_SF_GtoH;

    std::unique_ptr<Hists> n_input_events_hist;
    std::vector<std::unique_ptr<Hists>> v_hists;
    std::vector<std::unique_ptr<Hists>> v_hists_post;
    std::vector<std::unique_ptr<Hists>> v_hists_post_trig_Ele50_Photo;


    double lumi_tot;
    double lumi1;
    double lumi2;

};


VLQSemiLepPreSelEleTriggerEff::VLQSemiLepPreSelEleTriggerEff(Context & ctx) {
    version = ctx.get("dataset_version", "");
    type = ctx.get("dataset_type", "");
    auto data_dir_path = ctx.get("data_dir_path");

    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }

    ctx.set("lumi_file", "/nfs/dust/cms/user/schumas/ANALYSIS/80X_v4/CMSSW_8_0_24_patch1/src/UHH2/common/data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.root");

    n_input_events_hist.reset(new NInputEventsHist(ctx));

    // use centrally managed PU reweighting, jet corrections, jet lepton cleaning, jet smearing ....
    CommonModules* commonObjectCleaning = new CommonModules();
    commonObjectCleaning->set_jet_id(AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(25.,7.)));
    commonObjectCleaning->set_electron_id(AndId<Electron>(ElectronID_Spring16_tight_noIso,PtEtaCut(55.0, 2.5)));
    commonObjectCleaning->set_muon_id(AndId<Muon>(MuonIDTight(),PtEtaCut(55., 2.4)));
    commonObjectCleaning->switch_jetlepcleaner(true);
    commonObjectCleaning->switch_jetPtSorter(true);
    commonObjectCleaning->disable_jersmear();
    commonObjectCleaning->init(ctx);
    common_modules_with_lumi_sel.reset(commonObjectCleaning);

    if (type == "MC") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_L123_AK8PFchs_MC, 
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
    } else {
      if (version == "DataSingleMuB" || version == "DataSingleMuC" || version == "DataSingleMuD") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
      }
      if (version == "DataSingleMuE" || version == "DataSingleMuF") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
      }
      if (version == "DataSingleMuG") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
      }
      if (version == "DataSingleMuH") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
      }
      }

     v_pre_modules.emplace_back(new MCPileupReweight(ctx));
     v_pre_modules.emplace_back(new PrimaryLepton(ctx, "PrimaryLepton", prim_lep_min_ele_pt, prim_lep_min_mu_pt));

     v_pre_modules.emplace_back(new MCMuonScaleFactor(ctx,
        data_dir_path + "MuonID_EfficienciesAndSF_average_RunBtoH.root",
        "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1., "id", "nominal", "prim_mu_coll"));

     v_pre_modules.emplace_back(new MCElecScaleFactor(ctx,
       data_dir_path + "egammaEffi.txt_EGM2D_CutBased_Tight_ID.root",1., "eleid", "nominal"));

     Mu_SF_BtoF.reset(new MCMuonScaleFactor(ctx,
        data_dir_path + "MuonTrigger_EfficienciesAndSF_RunBtoF.root",
        "Mu50_OR_TkMu50_PtEtaBins", 1., "trg", "nominal", "prim_mu_coll"));
     Mu_SF_GtoH.reset(new MCMuonScaleFactor(ctx,
        data_dir_path + "MuonTrigger_EfficienciesAndSF_RunGtoH.root",
        "Mu50_OR_TkMu50_PtEtaBins", 1., "trg", "nominal", "prim_mu_coll"));

    lumi_tot = string2double(ctx.get("target_lumi"));
    lumi1 = 19870.;//0.622/fb in 2016 data                                                                                                      
    lumi2 = lumi_tot - lumi1;


    v_pre_modules.emplace_back(new EventWeightOutputHandle(ctx, "weight"));
    v_pre_modules.emplace_back(new LeadingJetPtProducer(ctx));
    v_pre_modules.emplace_back(new STCalculator(ctx, "ST", JetId(PtEtaCut(30., 2.4))));
    v_pre_modules.emplace_back(new PrimaryLeptonInfoProducer(ctx));
    v_pre_modules.emplace_back(new TwoDCutProducer(ctx));



    SelItemsHelper sel_helper(SEL_ITEMS_PRESEL, ctx);
    // sel_helper.write_cuts_to_texfile();
    sel_module.reset(new SelectionProducer(ctx, sel_helper));

    trigger_sel_Mu50.reset(new TriggerSelection("HLT_Mu50_v*"));
    trigger_sel_TkMu50.reset(new TriggerSelection("HLT_TkMu50_v*"));
    trigger_sel_Ele50.reset(new TriggerSelection("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*"));
    trigger_sel_Photo.reset(new TriggerSelection("HLT_Photon175_v*"));

    ele_sel.reset(new NElectronSelection(1));
    mu_sel.reset(new NMuonSelection(1));

    // 2. setup histograms
    auto * nm1_hists = new Nm1SelHists(ctx, "Nm1Selection", sel_helper);
    auto * cf_hists = new VLQ2HTCutflow(ctx, "Cutflow", sel_helper);

    v_hists.emplace_back(nm1_hists);
    v_hists.emplace_back(cf_hists);
    v_hists_post.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPost")); 
    v_hists_post_trig_Ele50_Photo.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPostTrig_Ele50_Photo"));


    // append 2D cut
    unsigned pos_2d_cut = 0;
    sel_module->insert_selection(pos_2d_cut, new TwoDCutSel(ctx, DR_2D_CUT_PRESEL, DPT_2D_CUT_PRESEL));
    nm1_hists->insert_hists(pos_2d_cut, new TwoDCutHist(ctx, "Nm1Selection"));
    cf_hists->insert_step(pos_2d_cut, "2D cut");

    // general histograms
    v_hists_post.emplace_back(new HistCollector(ctx, "EventHistsPost", type == "MC"));
    v_hists_post_trig_Ele50_Photo.emplace_back(new HistCollector(ctx, "EventHistsPostTrig_Ele50_Photo", type == "MC"));


}


bool VLQSemiLepPreSelEleTriggerEff::process(Event & event) {

    n_input_events_hist->fill(event);

 
    // common modules & check lumi selection on data
    if (!common_modules_with_lumi_sel->process(event)) {
        return false;
    }

    // run all event modules
    for (auto & mod : v_pre_modules) {
        mod->process(event);
    }

    double w_wo_HLT = event.weight;
    // muon-HLT eff                                                                                                                             
    Mu_SF_BtoF->process(event);
    double w1 = event.weight;
    event.weight = w_wo_HLT;
    Mu_SF_GtoH->process(event);
    double w2 = event.weight;
    double w = (lumi1*w1+lumi2*w2)/(lumi_tot);
    event.weight = w;


    // run selection
    bool all_accepted = sel_module->process(event);
    bool trigger_Mu50 = trigger_sel_Mu50->passes(event);
    bool trigger_TkMu50 = false;
    if (version != "DataSingleMuB"){
      trigger_TkMu50 = trigger_sel_TkMu50->passes(event); }
    bool trigger_Ele50 = trigger_sel_Ele50->passes(event);
    bool trigger_Photo = trigger_sel_Photo->passes(event);

    bool sel_lep = (ele_sel->passes(event) && mu_sel->passes(event));
    //bool sel_lep = (mu_sel->passes(event));
    if (trigger_TkMu50){
      std::cout << "all_accepted=" << all_accepted << "  trigger_Mu50=" << trigger_Mu50 << "  trigger_TkMu50=" << trigger_TkMu50 << "  trigger_Ele50=" << trigger_Ele50 << "  trigger_Photo=" << trigger_Photo << "  sel_lep=" << sel_lep << std::endl;
    }
    // fill ctrl hists
    for (auto & h : v_hists) {
        h->fill(event);
    }
    if(sel_lep && all_accepted){
      if (version == "DataSingleMuB")
	{ 
	if (trigger_Mu50){
	  for (auto & h : v_hists_post) {
	    h->fill(event);
	  }
	  if (trigger_Photo || trigger_Ele50) {
	    for (auto & h : v_hists_post_trig_Ele50_Photo) { 
	      h->fill(event);
	    }
	  }
	}
      }
      else
	{
	  if (trigger_TkMu50 || trigger_Mu50){
	    for (auto & h : v_hists_post) {
	      h->fill(event);
	    }
	    if (trigger_Photo || trigger_Ele50) {
	      for (auto & h : v_hists_post_trig_Ele50_Photo) { 
		h->fill(event);
	      }
	    }
	  }
	}
    }	


    return all_accepted;
}

UHH2_REGISTER_ANALYSIS_MODULE(VLQSemiLepPreSelEleTriggerEff)
