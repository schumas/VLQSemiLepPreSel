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
//#include "UHH2/VLQSemiLepPreSel/include/VLQToHiggsPairProdHists.h"
#include "UHH2/VLQSemiLepPreSel/include/VLQGenHists.h"


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

    Event::Handle<FlavorParticle> h_primlep;
    Event::Handle<double> h_st;

    std::vector<std::unique_ptr<Hists>> v_hists;
    std::vector<std::unique_ptr<Hists>> v_hists_post;
};


VLQSemiLepPreSel::VLQSemiLepPreSel(Context & ctx):
    h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")),
    h_st(ctx.get_handle<double>("ST"))
{

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
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new JetCorrector(JERFiles::PHYS14_L123_MC)));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new JetLeptonCleaner(JERFiles::PHYS14_L123_MC)));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new JetCleaner(30.0, 7.0)));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new JetPtSorter()));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new PrimaryLepton(ctx)));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new PartonHT(ctx.get_handle<double>("parton_ht"))));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new STCalculator(ctx)));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new NBTagProducer(ctx)));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new NHTagProducer(ctx, "patJetsCa15CHSJetsFilteredPacked")));
    v_pre_modules.push_back(std::unique_ptr<AnalysisModule>(new TopTagCalculator(ctx.get_handle<int>("n_toptags"))));

    v_hists.push_back(std::unique_ptr<Hists>(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPre")));
    //v_hists.push_back(std::unique_ptr<Hists>(new HistCollector(ctx, "EventsHistsPre")));
    v_hists.push_back(std::unique_ptr<Hists>(new VLQGenHists(ctx, "GenHistsPre")));

    v_hists_post.push_back(std::unique_ptr<Hists>(new VLQSemiLepPreSelHists(ctx, "PreSelCtrlPost")));
    //v_hists_post.push_back(std::unique_ptr<Hists>(new HistCollector(ctx, "EventsHistsPost")));
    //v_hists_post.push_back(std::unique_ptr<Hists>(new VLQGenHists(ctx, "GenHistsPost")));
}


bool VLQSemiLepPreSel::process(Event & event) {

    // reject event if no jets or lepton
    if (!event.jets->size()) {
        return false;
    }
    if (!(event.muons->size() || event.electrons->size())) {
        return false;
    }

    // run all event modules
    for (auto & mod : v_pre_modules) {
         mod->process(event);
    }

    // test again + good lepton
    if (!event.jets->size()) {
        return false;
    }
    if (!(event.muons->size() || event.electrons->size())) {
        return false;
    }
    if (!event.is_valid(h_primlep)) {
        return false;
    }

    // fill ctrl hists
    for (auto & h : v_hists) {
        h->fill(event);
    }

    // decide
    bool accept = event.get(h_primlep).pt() >= 50.
                  && event.jets->at(0).pt() >= 200.
                  && event.get(h_st) >= 500.;

    // fill ctrl hists after selection
    if (accept) {
        for (auto & h : v_hists_post) {
            h->fill(event);
        }
    }

    return accept;
}

UHH2_REGISTER_ANALYSIS_MODULE(VLQSemiLepPreSel)
