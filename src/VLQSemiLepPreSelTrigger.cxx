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

class HiggsJetBuilder: public AnalysisModule {
public:
    explicit HiggsJetBuilder(Context & ctx):
        h_topjets       (ctx.get_handle<vector<TopJet>>("topjets")),
        h_softdrop_jets (ctx.get_handle<vector<TopJet>>("patJetsAk8CHSJetsSoftDropPacked_daughters")),
        h_output        (ctx.get_handle<vector<TopJet>>("combined_ak8_jets")) {}

    virtual bool process(Event & event) override {
        const auto & topjets = event.get(h_topjets);
        const auto & sd_jets = event.get(h_softdrop_jets);

        vector<TopJet> comb_jets;
        for (auto top_j : topjets) {  // making a copy here
            for (const auto & sd_j : sd_jets) {
                if (deltaR(top_j.v4(), sd_j.v4()) < 0.5) {
                    auto tj_P2 = top_j.v4().P2();
                    auto sd_mass = sd_j.v4().mass();
                    auto nrg = sqrt(sd_mass*sd_mass + tj_P2);
                    top_j.set_energy(nrg);  // set mass through energy
                    comb_jets.push_back(top_j);
                }
            }
        }
        event.set(h_output, comb_jets);
        return true;
    }

private:
    Event::Handle<vector<TopJet>> h_topjets;
    Event::Handle<vector<TopJet>> h_softdrop_jets;
    Event::Handle<vector<TopJet>> h_output;
};  // HiggsJetBuilder



class VLQSemiLepPreSelTrigger: public AnalysisModule {
public:
    explicit VLQSemiLepPreSelTrigger(Context & ctx);
    virtual bool process(Event & event) override;

private:
    std::string version;
    std::string type;
    // modules for setting up collections and cleaning
    std::vector<std::unique_ptr<AnalysisModule>> v_pre_modules;
    unique_ptr<SelectionProducer> sel_module;
    unique_ptr<AnalysisModule> leptonic_decay_checker;
    unique_ptr<AnalysisModule> common_modules_with_lumi_sel;
    std::unique_ptr<Selection> trigger_sel_Ele50;
    std::unique_ptr<Selection> trigger_sel_Mu50;
    std::unique_ptr<Selection> trigger_sel_TkMu50;
    std::unique_ptr<Selection> trigger_sel_Photo;
    unique_ptr<AnalysisModule> Mu_SF_BtoF;
    unique_ptr<AnalysisModule> Mu_SF_GtoH;

    std::unique_ptr<Hists> n_input_events_hist;
    std::vector<std::unique_ptr<Hists>> v_hists;
    std::vector<std::unique_ptr<Hists>> v_hists_post;
    std::vector<std::unique_ptr<Hists>> v_hists_post_trig_Mu50;
    std::vector<std::unique_ptr<Hists>> v_hists_post_trig_TkMu50;
    std::vector<std::unique_ptr<Hists>> v_hists_post_trig_Ele50;
    std::vector<std::unique_ptr<Hists>> v_hists_post_trig_Photo;

    double lumi_tot;
    double lumi1;
    double lumi2;

};


VLQSemiLepPreSelTrigger::VLQSemiLepPreSelTrigger(Context & ctx) {
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
    commonObjectCleaning->set_muon_id(AndId<Muon>(MuonIDMedium(),PtEtaCut(55., 2.4)));
    commonObjectCleaning->switch_jetlepcleaner(false);
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
      if (version == "DataSingleEleB1" || version == "DataSingleEleB2" || version == "DataSingleEleC" || version == "DataSingleEleD" || version == "DataSingleMuB1" || version == "DataSingleMuB2" || version == "DataSingleMuC" || version == "DataSingleMuD") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
      }
      if (version == "DataSingleEleE" || version == "DataSingleEleF" || version == "DataSingleMuE" || version == "DataSingleMuF") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
      }
      if (version == "DataSingleEleG" || version == "DataSingleMuG") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
      }
      if (version == "DataSingleEleH" || version == "DataSingleMuH") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
      }
      }

     v_pre_modules.emplace_back(new MCPileupReweight(ctx));

     v_pre_modules.emplace_back(new MCMuonScaleFactor(ctx,
        data_dir_path + "MuonID_EfficienciesAndSF_average_RunBtoH.root",
        "MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta", 1., "id", "nominal", "prim_mu_coll"));

     v_pre_modules.emplace_back(new MCElecScaleFactor(ctx,
        data_dir_path + "egammaEffi.txt_EGM2D_CutBased_Tight_ID.root",1., "eleid", "nominal", "prim_ele_coll"));

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
    v_pre_modules.emplace_back(new HiggsJetBuilder(ctx));
    v_pre_modules.emplace_back(new TopJetCleaner(ctx, PtEtaCut(200., 2.4), "combined_ak8_jets"));

    v_pre_modules.emplace_back(new PrimaryLepton(ctx, "PrimaryLepton", prim_lep_min_ele_pt, prim_lep_min_mu_pt));
    v_pre_modules.emplace_back(new STCalculator(ctx, "ST", JetId(PtEtaCut(30., 2.4))));
    v_pre_modules.emplace_back(new CollectionSizeProducer<Jet>(ctx, "jets", "n_btags", JetId(CSVBTag(CSVBTag::WP_LOOSE))));
    v_pre_modules.emplace_back(new CollectionSizeProducer<TopJet>(ctx, 
        "patJetsAk8CHSJetsSoftDropPacked_daughters", "n_hcands", 
            TopJetId(HiggsTag(40., 99999., is_true<Jet>))));
    v_pre_modules.emplace_back(new LeadingJetPtProducer(ctx));
    v_pre_modules.emplace_back(new SubleadingJetPtProducer(ctx));
    v_pre_modules.emplace_back(new PrimaryLeptonInfoProducer(ctx));
    v_pre_modules.emplace_back(new TwoDCutProducer(ctx));
    v_pre_modules.emplace_back(new CollectionProducer<Jet>(ctx, "jets", "b_jets", JetId(AndId<Jet>(PtEtaCut(0., 2.4), CSVBTag(CSVBTag::WP_LOOSE)))));




    if (type == "MC") {
        v_pre_modules.emplace_back(new PartonHT(ctx.get_handle<double>("parton_ht")));
    }

    SelItemsHelper sel_helper(SEL_ITEMS_PRESEL, ctx);
    // sel_helper.write_cuts_to_texfile();
    sel_module.reset(new SelectionProducer(ctx, sel_helper));

    trigger_sel_Mu50.reset(new TriggerSelection("HLT_Mu50_v*"));
    trigger_sel_TkMu50.reset(new TriggerSelection("HLT_TkMu50_v*"));
    trigger_sel_Ele50.reset(new TriggerSelection("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v*"));
    trigger_sel_Photo.reset(new TriggerSelection("HLT_Photon175_v*"));


    // 2. setup histograms
    auto * nm1_hists = new Nm1SelHists(ctx, "Nm1Selection", sel_helper);
    auto * cf_hists = new VLQ2HTCutflow(ctx, "Cutflow", sel_helper);

    v_hists.emplace_back(nm1_hists);
    v_hists.emplace_back(cf_hists);
    v_hists_post.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPost"));
    v_hists_post_trig_Mu50.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPostTrig_Mu50"));
    v_hists_post_trig_TkMu50.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPostTrig_TkMu50"));
    v_hists_post_trig_Ele50.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPostTrig_Ele50"));
    v_hists_post_trig_Photo.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPostTrig_Photo"));

    // append 2D cut
    unsigned pos_2d_cut = 2;
    sel_module->insert_selection(pos_2d_cut, new TwoDCutSel(ctx, DR_2D_CUT_PRESEL, DPT_2D_CUT_PRESEL));
    nm1_hists->insert_hists(pos_2d_cut, new TwoDCutHist(ctx, "Nm1Selection"));
    cf_hists->insert_step(pos_2d_cut, "2D cut");

    // general histograms
    v_hists_post.emplace_back(new HistCollector(ctx, "EventHistsPost", type == "MC"));
    v_hists_post_trig_Mu50.emplace_back(new HistCollector(ctx, "EventHistsPostTrig_Mu50", type == "MC"));
    v_hists_post_trig_TkMu50.emplace_back(new HistCollector(ctx, "EventHistsPostTrig_TkMu50", type == "MC"));
    v_hists_post_trig_Ele50.emplace_back(new HistCollector(ctx, "EventHistsPostTrig_Ele50", type == "MC"));
    v_hists_post_trig_Photo.emplace_back(new HistCollector(ctx, "EventHistsPostTrig_Photo", type == "MC"));

    if (version.size() > 10 && version.substr(version.size() - 9, 100) == "_lepDecay") {
        cout << "USING leptonic_decay_checker!" << endl;
        leptonic_decay_checker.reset(new LeptonicDecayVLQ());
    }
}


bool VLQSemiLepPreSelTrigger::process(Event & event) {

    n_input_events_hist->fill(event);



    // signal: check for leptonic decay mode
    if (leptonic_decay_checker.get()
        && !leptonic_decay_checker->process(event)) {
        return false;
    }
 
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
    bool trigger_TkMu50 = trigger_sel_TkMu50->passes(event);
    bool trigger_Ele50 = trigger_sel_Ele50->passes(event);
    bool trigger_Photo = trigger_sel_Photo->passes(event);
    // fill ctrl hists
    for (auto & h : v_hists) {
        h->fill(event);
    }

    // fill ctrl hists after selection
    if (all_accepted) {
        for (auto & h : v_hists_post) {
            h->fill(event);
        }
    }

    if (all_accepted && trigger_Mu50) {
        for (auto & h : v_hists_post_trig_Mu50) {
            h->fill(event);
        }
    }

    if (all_accepted && trigger_TkMu50) {
        for (auto & h : v_hists_post_trig_TkMu50) {
            h->fill(event);
        }
    }

    if (all_accepted && trigger_Photo) {
      for (auto & h : v_hists_post_trig_Photo) { 
	    h->fill(event);
        }
    }

    if (all_accepted && trigger_Ele50) {
        for (auto & h : v_hists_post_trig_Ele50) {
            h->fill(event);
        }
    }


    return all_accepted;
}

UHH2_REGISTER_ANALYSIS_MODULE(VLQSemiLepPreSelTrigger)
