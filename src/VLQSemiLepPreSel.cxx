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

#include "UHH2/VLQSemiLepPreSel/include/VLQCommonModules.h"
#include "UHH2/VLQSemiLepPreSel/include/VLQSemiLepPreSelHists.h"
#include "UHH2/VLQSemiLepPreSel/include/VLQGenHists.h"
#include "UHH2/VLQSemiLepPreSel/include/EventHists.h"
#include "UHH2/VLQSemiLepPreSel/include/SelectionHists.h"
#include "UHH2/VLQSemiLepPreSel/include/VLQSLPS_selectionItems.h"


using namespace std;
using namespace uhh2;


template <typename TYPE>
static bool is_true(const TYPE &, const Event &) {
    return true;
}


class VLQSemiLepPreSel: public AnalysisModule {
public:
    explicit VLQSemiLepPreSel(Context & ctx);
    virtual bool process(Event & event) override;

private:
    std::string version;
    std::string type;
    // modules for setting up collections and cleaning
    std::vector<std::unique_ptr<AnalysisModule>> v_pre_modules;
    unique_ptr<SelectionProducer> sel_module;
    unique_ptr<AnalysisModule> leptonic_decay_checker;
    unique_ptr<AnalysisModule> common_modules_with_lumi_sel;

    std::unique_ptr<Hists> n_input_events_hist;
    std::vector<std::unique_ptr<Hists>> v_hists;
    std::vector<std::unique_ptr<Hists>> v_hists_post;
};


VLQSemiLepPreSel::VLQSemiLepPreSel(Context & ctx) {
    version = ctx.get("dataset_version", "");
    type = ctx.get("dataset_type", "");
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }

    if (version == "Run2015D_Ele") {
        ctx.set("lumi_file", "/afs/desy.de/user/t/tholenhe/xxl-af-cms/CMSSW_7_4_15_patch1/src/UHH2/common/data/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_NoBadBSRuns.root");
    } else if (version == "Run2015D_Mu") {
        ctx.set("lumi_file", "/afs/desy.de/user/t/tholenhe/xxl-af-cms/CMSSW_7_4_15_patch1/src/UHH2/common/data/Latest_2015_Golden_JSON.root");
    }

    n_input_events_hist.reset(new NInputEventsHist(ctx));

    // use centrally managed PU reweighting, jet corrections, jet lepton cleaning, jet smearing ....
    CommonModules* commonObjectCleaning = new CommonModules();
    commonObjectCleaning->set_jet_id(AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(25.,7.)));
    commonObjectCleaning->set_electron_id(AndId<Electron>(ElectronID_Spring15_25ns_tight_noIso,PtEtaCut(50.0, 2.5)));
    commonObjectCleaning->set_muon_id(AndId<Muon>(MuonIDMedium(),PtEtaCut(47., 2.1)));
    commonObjectCleaning->switch_jetlepcleaner(true);
    commonObjectCleaning->switch_jetPtSorter(true);
    commonObjectCleaning->disable_jersmear();
    commonObjectCleaning->init(ctx);
    common_modules_with_lumi_sel.reset(commonObjectCleaning);

    if (type == "MC") {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer15_25ns_L123_AK8PFchs_MC, 
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer15_25ns_L123_AK4PFchs_MC,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
    } else {
        v_pre_modules.emplace_back(new GenericTopJetCorrector(ctx,
            JERFiles::Summer15_25ns_L123_AK8PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
        v_pre_modules.emplace_back(new GenericSubJetCorrector(ctx,
            JERFiles::Summer15_25ns_L123_AK4PFchs_DATA,
                "patJetsAk8CHSJetsSoftDropPacked_daughters"));
    }

    v_pre_modules.emplace_back(new EventWeightOutputHandle(ctx, "weight"));
    v_pre_modules.emplace_back(new PrimaryLepton(ctx, "PrimaryLepton", prim_lep_min_ele_pt, prim_lep_min_mu_pt));
    v_pre_modules.emplace_back(new STCalculator(ctx, "ST", JetId(PtEtaCut(30., 2.4))));
    v_pre_modules.emplace_back(new CollectionSizeProducer<Jet>(ctx, "jets", "n_btags", JetId(CSVBTag(CSVBTag::WP_LOOSE))));
    v_pre_modules.emplace_back(new CollectionSizeProducer<TopJet>(ctx, 
        "patJetsAk8CHSJetsSoftDropPacked_daughters", "n_hcands", 
            TopJetId(HiggsTag(40., 99999., is_true<Jet>))));
    v_pre_modules.emplace_back(new LeadingJetPtProducer(ctx));
    v_pre_modules.emplace_back(new PrimaryLeptonInfoProducer(ctx));
    v_pre_modules.emplace_back(new TwoDCutProducer(ctx));

    if (version == "TTbar") {
        v_pre_modules.emplace_back(new TTbarGenProducer(ctx, "ttbargen", false));
        v_pre_modules.emplace_back(new TopPtWeight(ctx, "ttbargen", 0.156, -0.00137, "weight_ttbar", false));
        v_hists.emplace_back(new TopPtWeightHist(ctx, "TTbarReweight", "weight_ttbar"));
    }

    // if (version == "Run2015D_Mu") {
    //     v_pre_modules.emplace_back(new TriggerAcceptProducer(ctx,
    //         PRESEL_TRIGGER_PATHS_DATA, "trigger_accept"));
    // } else if (version == "Run2015D_Ele") {
    //     v_pre_modules.emplace_back(new TriggerAcceptProducer(ctx,
    //         PRESEL_TRIGGER_PATHS_DATA, PRESEL_TRIGGER_PATHS_DATA_ELE_VETO, "trigger_accept"));
    // } else if (version == "Run2015D_Had") {
    //     v_pre_modules.emplace_back(new TriggerAcceptProducer(ctx,
    //         PRESEL_TRIGGER_PATHS_DATA, PRESEL_TRIGGER_PATHS_DATA_HAD_VETO, "trigger_accept"));
    // } else {
    //     v_pre_modules.emplace_back(new TriggerAcceptProducer(ctx,
    //         PRESEL_TRIGGER_PATHS, "trigger_accept"));
    // }

    if (type == "MC") {
        v_pre_modules.emplace_back(new PartonHT(ctx.get_handle<double>("parton_ht")));
    }

    SelItemsHelper sel_helper(SEL_ITEMS_PRESEL, ctx);
    // sel_helper.write_cuts_to_texfile();
    sel_module.reset(new SelectionProducer(ctx, sel_helper));

    // 2. setup histograms
    auto * nm1_hists = new Nm1SelHists(ctx, "Nm1Selection", sel_helper);
    auto * cf_hists = new VLQ2HTCutflow(ctx, "Cutflow", sel_helper);

    v_hists.emplace_back(nm1_hists);
    v_hists.emplace_back(cf_hists);
    v_hists_post.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPost"));

    // append 2D cut
    // unsigned pos_2d_cut = 2;
    // sel_module->insert_selection(pos_2d_cut, new TwoDCutSel(ctx, DR_2D_CUT_PRESEL, DPT_2D_CUT_PRESEL));
    // nm1_hists->insert_hists(pos_2d_cut, new TwoDCutHist(ctx, "Nm1Selection"));
    // cf_hists->insert_step(pos_2d_cut, "2D cut");

    // general histograms
    v_hists_post.emplace_back(new HistCollector(ctx, "EventHistsPost", type == "MC"));

    if (version.size() > 10 && version.substr(version.size() - 9, 100) == "_lepDecay") {
        cout << "USING leptonic_decay_checker!" << endl;
        leptonic_decay_checker.reset(new LeptonicDecayVLQ());
    }
}


bool VLQSemiLepPreSel::process(Event & event) {

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

    // run selection
    bool all_accepted = sel_module->process(event);

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

    return all_accepted;
}

UHH2_REGISTER_ANALYSIS_MODULE(VLQSemiLepPreSel)
