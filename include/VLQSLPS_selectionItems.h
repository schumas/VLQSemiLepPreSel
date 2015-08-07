#pragma once

#include <vector>
#include "TH2F.h"

#include "UHH2/common/include/Utils.h"
#include "UHH2/VLQSemiLepPreSel/include/SelectionItem.h"


typedef SelectionItemData<int>      SelDatI;
typedef SelectionItemData<float>    SelDatF;
typedef SelectionItemData<double>   SelDatD;


static const float DR_2D_CUT_PRESEL = 2.5;
static const float DPT_2D_CUT_PRESEL = 40.0;


static const vector<shared_ptr<SelectionItem>> SEL_ITEMS_PRESEL {
    shared_ptr<SelectionItem>(new SelDatI("trigger_accept",    "trigger accept",                    2, -.5, 1.5       ,1      )),
    shared_ptr<SelectionItem>(new SelDatF("leading_jet_pt",    "leading jet p_{T} / 20 GeV",        75, 0, 1500       ,100    )),
    shared_ptr<SelectionItem>(new SelDatF("primary_lepton_pt", "primary lepton p_{T} / 20 GeV",     75, 0, 1500       ,25     )),
    shared_ptr<SelectionItem>(new SelDatD("ST",                "ST / 50 GeV",                       100, 0, 5000      ,400    )),
    shared_ptr<SelectionItem>(new SelDatI("n_btags",           "number of loose btags",             11, -.5, 10.5     ,1      )),
};


static const vector<std::string> PRESEL_TRIGGER_PATHS {
    // "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*",
    // "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
    // "HLT_Ele32_eta2p1_WP75_Gsf_v*",

    // "HLT_Mu40_eta2p1_PFJet200_PFJet50_v*",
    "HLT_Mu45_eta2p1_v*",
    "HLT_Mu50_v*",
    "HLT_IsoMu24_eta2p1_v*",
    "HLT_IsoMu27_v*",

    // "HLT_PFHT800Emu_v*",
};

static const vector<std::string> PRESEL_TRIGGER_PATHS_DATA {
    // "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*",
    // "HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
    // "HLT_Ele32_eta2p1_WP75_Gsf_v*",

    // "HLT_Mu40_eta2p1_PFJet200_PFJet50_v*",
    "HLT_Mu45_eta2p1_v*",
    "HLT_Mu50_v*",
    "HLT_IsoMu24_eta2p1_v*",
    "HLT_IsoMu27_v*",

    // "HLT_PFHT800_v*",
};

