#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "TH1F.h"
#include "TH1D.h"

#include "UHH2/VLQSemiLepPreSel/include/HandleHist.h"
#include "UHH2/VLQSemiLepPreSel/include/SelectionItem.h"


class Nm1SelHists: public Hists {
public:
    explicit Nm1SelHists(Context & ctx,
                         const string & dir,
                         const SelItemsHelper & sel_helper):
        Hists(ctx, dir),
        h_sel_res(ctx.get_handle<vector<bool>>("sel_accept"))
    {
        sel_helper.fill_hists_vector(v_hists, dir);
    }

    void add_hists(Hists * hists) {
        v_hists.emplace_back(hists);
    }

    virtual void fill(const Event & event) override {
        const auto & v_accept = event.get(h_sel_res);
        for (unsigned i=0; i<v_hists.size(); ++i) {
            bool accept_nm1 = true;
            for (unsigned j=0; j<v_accept.size(); ++j) {
                if (i==j) {
                    continue;
                }
                if (!v_accept[j]) {
                    accept_nm1 = false;
                    break;
                }
            }
            if (accept_nm1) {
                v_hists[i]->fill(event);
            }
        }
    }

private:
    Event::Handle<vector<bool>> h_sel_res;
    vector<unique_ptr<Hists>> v_hists;
};


class VLQ2HTCutflow: public Hists {
public:
    VLQ2HTCutflow(Context & ctx,
                  const string & dir,
                  const SelItemsHelper & sel_helper):
        Hists(ctx, dir),
        h_sel(ctx.get_handle<vector<bool>>("sel_accept")),
        v_names(sel_helper.get_item_names()),
        h(book<TH1D>("cutflow", "", 1, 0, -1))
    {
        h->SetBit(TH1::kCanRebin);
        h->Fill("input", 1e-7);
        for (const string & name : v_names) {
            h->Fill(name.c_str(), 1e-7);
        }
    }

    void add_step(const string & name) {
        v_names.push_back(name);
        h->Fill(name.c_str(), 1e-7);
    }

    virtual void fill(const uhh2::Event & e) override {
        float w = e.weight;
        const auto & sel = e.get(h_sel);
        h->Fill("input", w);
        for (unsigned i = 0; i < v_names.size(); ++i) {
            if (sel[i]) {
                h->Fill(v_names[i].c_str(), w);
            } else {
                break;
            }
        }
    }

private:
    Event::Handle<vector<bool>> h_sel;
    vector<string> v_names;
    TH1D * h;
};
