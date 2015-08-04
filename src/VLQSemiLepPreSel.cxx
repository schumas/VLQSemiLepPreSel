#include <algorithm>
#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/JetIds.h"
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


class VLQSemiLepPreSel: public AnalysisModule {
public:
    explicit VLQSemiLepPreSel(Context & ctx);
    virtual bool process(Event & event) override;

private:
    std::string version;
    // modules for setting up collections and cleaning
    std::vector<std::unique_ptr<AnalysisModule>> v_pre_modules;
    unique_ptr<SelectionProducer> sel_module;
    unique_ptr<AnalysisModule> leptonic_decay_checker;

    std::vector<std::unique_ptr<Hists>> v_hists;
    std::vector<std::unique_ptr<Hists>> v_hists_post;
};


VLQSemiLepPreSel::VLQSemiLepPreSel(Context & ctx) {
    version = ctx.get("dataset_version", "");
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }

    // use centrally managed PU reweighting, jet corrections, jet lepton cleaning, jet smearing ....
    CommonModules* commonOjectCleaning = new CommonModules();
    commonOjectCleaning->set_jet_id(PtEtaCut(30.0,7.0));
    commonOjectCleaning->set_electron_id(AndId<Electron>(ElectronID_PHYS14_25ns_medium_noIso,PtEtaCut(20.0, 2.4)));
    commonOjectCleaning->set_muon_id(AndId<Muon>(MuonIDTight(),PtEtaCut(20.0, 2.4)));
    commonOjectCleaning->switch_jetlepcleaner(true);
    commonOjectCleaning->switch_jetPtSorter(true);
    commonOjectCleaning->init(ctx);
    v_pre_modules.emplace_back(commonOjectCleaning);

    v_pre_modules.emplace_back(new PrimaryLepton(ctx));
    v_pre_modules.emplace_back(new PartonHT(ctx.get_handle<double>("parton_ht")));
    v_pre_modules.emplace_back(new STCalculator(ctx));
    v_pre_modules.emplace_back(new CollectionSizeProducer<Jet>(ctx, "jets", "n_btags", JetId(CSVBTag(CSVBTag::WP_LOOSE))));
    v_pre_modules.emplace_back(new CollectionSizeProducer<Jet>(ctx, "jets", "n_btags", JetId(CSVBTag(CSVBTag::WP_LOOSE))));
    v_pre_modules.emplace_back(new CollectionSizeProducer<TopJet>(ctx, "patJetsAk8CHSJetsSoftDropPacked_daughters", "n_higgstags", TopJetId(HiggsTag())));
    v_pre_modules.emplace_back(new CollectionSizeProducer<TopJet>(ctx, "topjets", "n_toptags", TopJetId(CMSTopTag())));
    v_pre_modules.emplace_back(new LeadingJetPtProducer(ctx));
    v_pre_modules.emplace_back(new LeptonPtProducer(ctx));
    v_pre_modules.emplace_back(new TwoDCutProducer(ctx));
    v_pre_modules.emplace_back(new TriggerAcceptProducer(ctx, PRESEL_TRIGGER_PATHS, "trigger_accept"));

    SelItemsHelper sel_helper(SEL_ITEMS_PRESEL, ctx);
    // sel_helper.write_cuts_to_texfile();
    sel_module.reset(new SelectionProducer(ctx, sel_helper));

    // 2. setup histograms
    auto * nm1_hists = new Nm1SelHists(ctx, "Nm1Selection", sel_helper);
    auto * cf_hists = new VLQ2HTCutflow(ctx, "Cutflow", sel_helper);

    v_hists.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPre"));
    v_hists.emplace_back(new HistCollector(ctx, "EventHistsPre"));
    // v_hists.emplace_back(new VLQGenHists(ctx, "GenHistsPre"));
    v_hists.emplace_back(nm1_hists);
    v_hists.emplace_back(cf_hists);
    sel_helper.fill_hists_vector(v_hists, "NoSelection");

    v_hists_post.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPost"));
    v_hists_post.emplace_back(new HistCollector(ctx, "EventHistsPost"));
    // v_hists_post.emplace_back(new VLQGenHists(ctx, "GenHistsPost"));

    // append 2D cut
    unsigned pos_2d_cut = 4;
    sel_module->insert_selection(pos_2d_cut, new TwoDCutSel(ctx, DR_2D_CUT_PRESEL, DPT_2D_CUT_PRESEL));
    nm1_hists->insert_hists(pos_2d_cut, new TwoDCutHist(ctx, "Nm1Selection"));
    cf_hists->insert_step(pos_2d_cut, "2D cut");
    v_hists.insert(v_hists.begin() + pos_2d_cut, move(unique_ptr<Hists>(new TwoDCutHist(ctx, "NoSelection"))));

    if (version.substr(version.size() - 4, 100) == "_lepDecay") {
        leptonic_decay_checker.reset(new LeptonicDecayVLQ());
    }
}


bool VLQSemiLepPreSel::process(Event & event) {

    // only for signal in leptonic decay mode
    if (leptonic_decay_checker.get()
        && !leptonic_decay_checker->process(event)) {
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
