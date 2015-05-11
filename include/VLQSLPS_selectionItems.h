#pragma once

#include <vector>
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


class TwoDCutSel: public Selection {
public:
    explicit TwoDCutSel(Context & ctx,
                        float min_dr,
                        float min_ptrel):
        h_prim_lep(ctx.get_handle<FlavorParticle>("PrimaryLepton")),
        min_dr_(min_dr),
        min_ptrel_(min_ptrel) {}

    virtual bool passes(const Event & e) override {
        if (!e.is_valid(h_prim_lep)) {
            return false;
        }
        auto prim_lep = e.get(h_prim_lep);
        float dr, pt;
        std::tie(dr, pt) = drmin_pTrel(prim_lep, *e.jets);
        return (dr > min_dr_) || (pt > min_ptrel_);
    }

private:
    Event::Handle<FlavorParticle> h_prim_lep;
    float min_dr_;
    float min_ptrel_;
};


class TwoDCutHist: public Hists {
public:
    explicit TwoDCutHist(Context & ctx,
                         const string & dirname):
        Hists(ctx, dirname),
        h_prim_lep(ctx.get_handle<FlavorParticle>("PrimaryLepton")),
        hist(book<TH2F>("TwoDCut",
                        ";min #DeltaR(lep., jet);min p_{T,rel}(lep., jet)",
                        100, 0., 5., 100, 0., 500.)) {}

    virtual void fill(const Event & e) override {
        if (!e.is_valid(h_prim_lep)) {
            return;
        }
        auto prim_lep = e.get(h_prim_lep);
        float dr, pt;
        std::tie(dr, pt) = drmin_pTrel(prim_lep, *e.jets);
        hist->Fill(dr, pt, e.weight);
    }

private:
    Event::Handle<FlavorParticle> h_prim_lep;
    TH2F * hist;
};