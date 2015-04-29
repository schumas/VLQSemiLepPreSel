#include <algorithm>
#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"

#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronIds.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/JetCorrections.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/PartonHT.h"
#include "UHH2/common/include/EventVariables.h"

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
    unique_ptr<AnalysisModule> sel_module;

    std::vector<std::unique_ptr<Hists>> v_hists;
    std::vector<std::unique_ptr<Hists>> v_hists_post;
};


VLQSemiLepPreSel::VLQSemiLepPreSel(Context & ctx) {

    // If needed, access the configuration of the module here, e.g.:
    // string testvalue = ctx.get("TestKey", "<not set>");
    // cout << "TestKey in the configuration was: " << testvalue << endl;
    version = ctx.get("dataset_version", "");
    
    // If running in SFrame, the keys "dataset_version", "dataset_type" and "dataset_lumi"
    // are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
    
    // 1. setup modules to prepare the event.
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new ElectronCleaner(
        AndId<Electron>(
            ElectronID_PHYS14_25ns_medium_noIso,
            PtEtaCut(20.0, 2.4)
        )
    )));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new MuonCleaner(
        AndId<Muon>(
            MuonIDTight(),
            PtEtaCut(20.0, 2.4)
        )
    )));
    v_pre_modules.emplace_back(new JetCorrector(JERFiles::PHYS14_L123_MC));
    v_pre_modules.emplace_back(new JetLeptonCleaner(JERFiles::PHYS14_L123_MC));
    v_pre_modules.emplace_back(new JetCleaner(30.0, 7.0));
    v_pre_modules.emplace_back(new JetPtSorter());
    v_pre_modules.emplace_back(new PrimaryLepton(ctx));
    v_pre_modules.emplace_back(new PartonHT(ctx.get_handle<double>("parton_ht")));
    v_pre_modules.emplace_back(new HTCalculator(ctx));
    v_pre_modules.emplace_back(new STCalculator(ctx));
    v_pre_modules.emplace_back(new NBTagProducer(ctx));
    v_pre_modules.emplace_back(new NTaggedTopJetProducer(ctx, TopJetId(HiggsTag()), "n_higgstags", "patJetsCa15CHSJetsFilteredPacked"));
    v_pre_modules.emplace_back(new NTaggedTopJetProducer(ctx, TopJetId(CMSTopTag()), "n_toptags"));
    v_pre_modules.emplace_back(new LeadingJetPtProducer(ctx));
    v_pre_modules.emplace_back(new LeptonPtProducer(ctx));

    SelItemsHelper sel_helper(SEL_ITEMS_PRESEL, ctx);
    sel_module.reset(new SelectionProducer(ctx, sel_helper));

    // 2. setup histograms
    sel_helper.fill_hists_vector(v_hists, "NoSelection");
    v_hists.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPre"));
    v_hists.emplace_back(new HistCollector(ctx, "EventHistsPre"));
    v_hists.emplace_back(new VLQGenHists(ctx, "GenHistsPre"));
    v_hists.emplace_back(new Nm1SelHists(ctx, "Nm1Selection", sel_helper));
    v_hists.emplace_back(new VLQ2HTCutflow(ctx, "Cutflow", sel_helper));

    v_hists_post.emplace_back(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPost"));
    v_hists_post.emplace_back(new HistCollector(ctx, "EventHistsPost"));
    v_hists_post.emplace_back(new VLQGenHists(ctx, "GenHistsPost"));
}


bool VLQSemiLepPreSel::process(Event & event) {

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
