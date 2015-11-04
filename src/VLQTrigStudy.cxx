#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>
#include <map>

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
#include "UHH2/VLQSemiLepPreSel/include/HandleHist.h"


using namespace std;
using namespace uhh2;


class VLQTrigHists: public Hists {
public:
    explicit VLQTrigHists(Context & ctx, const string & dir, float min_st):
        Hists(ctx, dir),
        st_denom(HandleHist<double>(ctx, dir+"/denom", "ST", ";ST;events / 50 GeV", 50, 0, 2500)),
        lep_pt_denom(HandleHist<float>(ctx, dir+"/denom", "primary_lepton_pt", ";primary lepton p_{T}; events / 20 GeV", 80, 0, 800)),
        lead_jet_pt_denom(HandleHist<float>(ctx, dir+"/denom", "leading_jet_pt", ";leading jet p_{T}; events / 20 GeV", 80, 0, 800)),
        sublead_jet_pt_denom(HandleHist<float>(ctx, dir+"/denom", "subleading_jet_pt", ";subleading jet p_{T}; events / 20 GeV", 80, 0, 800)),
        st_passing(HandleHist<double>(ctx, dir+"/passing", "ST", ";ST;events / 50 GeV", 50, 0, 2500)),
        lep_pt_passing(HandleHist<float>(ctx, dir+"/passing", "primary_lepton_pt", ";primary lepton p_{T}; events / 20 GeV", 80, 0, 800)),
        lead_jet_pt_passing(HandleHist<float>(ctx, dir+"/passing", "leading_jet_pt", ";leading jet p_{T}; events / 20 GeV", 80, 0, 800)),
        sublead_jet_pt_passing(HandleHist<float>(ctx, dir+"/passing", "subleading_jet_pt", ";subleading jet p_{T}; events / 20 GeV", 80, 0, 800)),
        hndl_st(ctx.get_handle<double>("ST")),
        hndl_lep_pt(ctx.get_handle<float>("primary_lepton_pt")),
        hndl_trg(ctx.get_handle<int>(dir)),
        st_cut(min_st)
    {}

    virtual void fill(const Event & e) override {
        auto st = e.get(hndl_st);

        // if (e.is_valid(hndl_lep_pt) && e.get(hndl_lep_pt) < 50.) {
        //     return;
        // }

        st_denom.fill(e);
        if (st > st_cut) {
            lep_pt_denom.fill(e);
            lead_jet_pt_denom.fill(e);
            sublead_jet_pt_denom.fill(e);            
        }
        if (e.get(hndl_trg)) {
            st_passing.fill(e);
            if (st > st_cut) {
                lep_pt_passing.fill(e);
                lead_jet_pt_passing.fill(e);
                sublead_jet_pt_passing.fill(e);            
            }
        }
    }

private:
    HandleHist<double> st_denom;
    HandleHist<float> lep_pt_denom;
    HandleHist<float> lead_jet_pt_denom;
    HandleHist<float> sublead_jet_pt_denom;
    HandleHist<double> st_passing;
    HandleHist<float> lep_pt_passing;
    HandleHist<float> lead_jet_pt_passing;
    HandleHist<float> sublead_jet_pt_passing;
    Event::Handle<double> hndl_st;
    Event::Handle<float> hndl_lep_pt;
    Event::Handle<int> hndl_trg;
    float st_cut;
};


class VLQTrigStudy: public AnalysisModule {
public:
    explicit VLQTrigStudy(Context & ctx);
    virtual bool process(Event & event) override;

private:
    std::string version;
    std::string type;
    // modules for setting up collections and cleaning
    std::vector<std::unique_ptr<AnalysisModule>> v_pre_modules;
    unique_ptr<AnalysisModule> leptonic_decay_checker_ele;
    unique_ptr<AnalysisModule> leptonic_decay_checker_mu;
    unique_ptr<AnalysisModule> common_modules_with_lumi_sel;
    unique_ptr<Selection> event_sel;

    std::vector<std::unique_ptr<Hists>> v_hists_ele;
    std::vector<std::unique_ptr<Hists>> v_hists_mu;

    Event::Handle<float> hndl_lep_pt;
};


VLQTrigStudy::VLQTrigStudy(Context & ctx) {
    version = ctx.get("dataset_version", "");
    type = ctx.get("dataset_type", "");
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }

    leptonic_decay_checker_ele.reset(new LeptonicDecayVLQ({11, -11}));
    leptonic_decay_checker_mu.reset(new LeptonicDecayVLQ({13, -13}));

    // use centrally managed PU reweighting, jet corrections, jet lepton cleaning, jet smearing ....
    // CommonModules* commonObjectCleaning = new CommonModules();
    // commonObjectCleaning->set_jet_id(AndId<Jet>(JetPFID(JetPFID::WP_LOOSE), PtEtaCut(30.0,3.6)));
    // commonObjectCleaning->set_electron_id(AndId<Electron>(ElectronID_Spring15_25ns_medium_noIso,PtEtaCut(20.0, 2.1)));
    // commonObjectCleaning->set_muon_id(AndId<Muon>(MuonIDTight(),PtEtaCut(20.0, 2.1)));
    // commonObjectCleaning->switch_jetlepcleaner(true);
    // commonObjectCleaning->switch_jetPtSorter(true);
    // commonObjectCleaning->disable_mcpileupreweight();
    // commonObjectCleaning->init(ctx);
    // common_modules_with_lumi_sel.reset(commonObjectCleaning);

    v_pre_modules.emplace_back(new PrimaryLepton(ctx, "PrimaryLepton", 20., 20.));
    v_pre_modules.emplace_back(new STCalculator(ctx, "ST", JetId(PtEtaCut(40., 2.4))));
    v_pre_modules.emplace_back(new SubleadingJetPtProducer(ctx));
    v_pre_modules.emplace_back(new LeadingJetPtProducer(ctx));
    v_pre_modules.emplace_back(new TwoDCutProducer(ctx));
    // v_pre_modules.emplace_back(new LeptonPtProducer(ctx));

    float min_st = version.substr(0,9) == "MC_TpTp_M" ? 700. : 400.;
    cout << "MINIMUM ST CUT:" << min_st << endl;

    map<string, vector<string>> items_ele = {
        {"trg_Ele32_eta2p1_WP75_Gsf",
            {"HLT_Ele32_eta2p1_WP75_Gsf_v*",}
        },
        // {"trg_Ele32_eta2p1_WP75_Gsf_OR_PFHT800",
        //     {"HLT_Ele32_eta2p1_WP75_Gsf_v*", "HLT_PFHT800Emu_v*",}
        // },
        {"trg_Ele105_CaloIdVT_GsfTrkIdT",
            {"HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",}
        },
        // {"trg_Ele105_CaloIdVT_GsfTrkIdT_OR_PFHT800",
        //     {"HLT_Ele105_CaloIdVT_GsfTrkIdT_v*", "HLT_PFHT800Emu_v*",}
        // },
        {"trg_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50",
            {"HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*",}
        },
        {"trg_Ele15_IsoVVVL_PFHT600_v",
            {"HLT_Ele15_IsoVVVL_PFHT600_v*",}
        },
    };
    map<string, vector<string>> items_mu = {
        {"trg_IsoMu24_eta2p1",
            {"HLT_IsoMu24_eta2p1_v*",}
        },
        // {"trg_IsoMu24_eta2p1_OR_PFHT800",
        //     {"HLT_IsoMu24_eta2p1_v*", "HLT_PFHT800Emu_v*",}
        // },
        {"trg_Mu45_eta2p1",
            {"HLT_Mu45_eta2p1_v*",}
        },
        // {"trg_Mu45_eta2p1_OR_PFHT800",
        //     {"HLT_Mu45_eta2p1_v*", "HLT_PFHT800Emu_v*",}
        // },
        {"trg_Mu40_eta2p1_PFJet200_PFJet50",
            {"HLT_Mu40_eta2p1_PFJet200_PFJet50_v*",}
        },
        {"trg_Mu15_IsoVVVL_PFHT600",
            {"HLT_Mu15_IsoVVVL_PFHT600_v*",}
        },
        {"trg_PFHT800",
            {"HLT_PFHT800Emu_v*",}
        },
    };
    for (auto & kv: items_ele) {
        v_pre_modules.emplace_back(new TriggerAcceptProducer(ctx, kv.second, kv.first));
        v_hists_ele.emplace_back(new VLQTrigHists(ctx, kv.first, min_st));
    }
    for (auto & kv: items_mu) {
        v_pre_modules.emplace_back(new TriggerAcceptProducer(ctx, kv.second, kv.first));
        v_hists_mu.emplace_back(new VLQTrigHists(ctx, kv.first, min_st));
    }

    event_sel.reset(new TwoDCutSel(ctx, 0.4, 40.));
    hndl_lep_pt = ctx.get_handle<float>("primary_lepton_pt");
}


bool VLQTrigStudy::process(Event & event) {

    bool is_ele = leptonic_decay_checker_ele->process(event);
    bool is_mu = leptonic_decay_checker_mu->process(event);

    static int all = 0;
    static int sub = 0;


    if (!(all % 1000)) {
        cout << "rate of events: " << sub/ (float) all << ", sub: " << sub << ", all: " << all << endl;
    }


    all += 1;
    // signal: check for leptonic decay mode
    if (is_ele == is_mu) {  // this is "not a xor b"
        return false;
    }

    // common modules & check lumi selection on data
    if (common_modules_with_lumi_sel.get() && !common_modules_with_lumi_sel->process(event)) {
        return false;
    }

    // run all event modules
    for (auto & mod : v_pre_modules) {
        mod->process(event);
    }

    if (event_sel.get() && !event_sel->passes(event)) {
        return false;
    }

    if (is_ele && event.electrons->size()) {
        event.set(hndl_lep_pt, event.electrons->at(0).pt()); 

        sub += 1;
        for (auto & hist: v_hists_ele) {
            hist->fill(event);
        }
    } else if (is_mu && event.muons->size()) {
        event.set(hndl_lep_pt, event.muons->at(0).pt());

        sub += 1;
        for (auto & hist: v_hists_mu) {
            hist->fill(event);
        }
    }

    return false;  // don't want output
}

UHH2_REGISTER_ANALYSIS_MODULE(VLQTrigStudy)
