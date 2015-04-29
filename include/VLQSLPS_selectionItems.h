#pragma once

#include <vector>
#include "UHH2/VLQSemiLepPreSel/include/SelectionItem.h"


typedef SelectionItemData<int>      SelDatI;
typedef SelectionItemData<float>    SelDatF;
typedef SelectionItemData<double>   SelDatD;


static const vector<shared_ptr<SelectionItem>> SEL_ITEMS_PRESEL {
    shared_ptr<SelectionItem>(new SelDatF("leading_jet_pt",    "leading jet p_{T}",        50, 100, 1500       ,200    )),
    shared_ptr<SelectionItem>(new SelDatF("primary_lepton_pt", "primary lepton p_{T}",     50, 100, 1500       ,50     )),
    shared_ptr<SelectionItem>(new SelDatD("ST",                "ST",                       100, 0, 5000        ,500    )),
};