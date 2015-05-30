#pragma once

#include <vector>
#include "TH2F.h"

#include "UHH2/common/include/Utils.h"
#include "UHH2/VLQSemiLepPreSel/include/SelectionItem.h"


typedef SelectionItemData<int>      SelDatI;
typedef SelectionItemData<float>    SelDatF;
typedef SelectionItemData<double>   SelDatD;


static const vector<shared_ptr<SelectionItem>> SEL_ITEMS_PRESEL {
    shared_ptr<SelectionItem>(new SelDatF("leading_jet_pt",    "leading jet p_{T}",        50, 0, 1500       ,100    )),
    shared_ptr<SelectionItem>(new SelDatF("primary_lepton_pt", "primary lepton p_{T}",     50, 0, 1500       ,50     )),
    shared_ptr<SelectionItem>(new SelDatD("ST",                "ST",                       100, 0, 5000      ,400    )),
    shared_ptr<SelectionItem>(new SelDatI("n_btags",           "number of loose btags",    11, -.5, 10.5     ,1      )),
};


namespace {

class TwoDCutSel: public Selection {
public:
    explicit TwoDCutSel(Context & ctx,
                        float min_dr,
                        float min_ptrel,
                        const string & dr_name = "TwoDCut_dr",
                        const string & pt_name = "TwoDCut_ptrel"):
        h_dr(ctx.get_handle<float>(dr_name)),
        h_pt(ctx.get_handle<float>(pt_name)),
        min_dr_(min_dr),
        min_ptrel_(min_ptrel) {}

    virtual bool passes(const Event & e) override {
        if (!e.is_valid(h_dr) || !e.is_valid(h_pt)) {
            return false;
        }
        return (e.get(h_dr) > min_dr_) || (e.get(h_pt) > min_ptrel_);
    }

private:
    Event::Handle<float> h_dr;
    Event::Handle<float> h_pt;
    float min_dr_;
    float min_ptrel_;
};


class TwoDCutHist: public Hists {
public:
    explicit TwoDCutHist(Context & ctx,
                         const string & dirname,
                         const string & dr_name = "TwoDCut_dr",
                         const string & pt_name = "TwoDCut_ptrel"):
        Hists(ctx, dirname),
        h_dr(ctx.get_handle<float>(dr_name)),
        h_pt(ctx.get_handle<float>(pt_name)),
        hist(book<TH2F>("TwoDCut",
                        ";min #DeltaR(lep., jet);min p_{T,rel}(lep., jet)",
                        200, 0., 1., 200, 0., 500.)) {}

    virtual void fill(const Event & e) override {
        if (!e.is_valid(h_dr) || !e.is_valid(h_pt)) {
            return;
        }
        hist->Fill(e.get(h_dr), e.get(h_pt), e.weight);
    }

private:
    Event::Handle<float> h_dr;
    Event::Handle<float> h_pt;
    TH2F * hist;
};

}