#pragma once

#include <assert.h>
#include <memory>

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "TH1F.h"
#include "TH1D.h"

#include "UHH2/VLQSemiLepPreSel/include/HandleHist.h"
#include "UHH2/VLQSemiLepPreSel/include/SelectionItem.h"


using namespace std;
using namespace uhh2;


class Nm1SelHists: public Hists {
public:
    explicit Nm1SelHists(Context & ctx,
                         const string & dir,
                         const SelItemsHelper & sel_helper):
        Hists(ctx, dir),
        h_sel_res(ctx.get_handle<vector<bool>>(sel_helper.s_vec_bool()))
    {
        sel_helper.fill_hists_vector(v_hists, dir);
    }

    void insert_hists(unsigned pos, Hists * hists) {
        v_hists.insert(v_hists.begin() + pos, move(unique_ptr<Hists>(hists)));
    }

    virtual void fill(const Event & event) override {
        const auto & v_accept = event.get(h_sel_res);
        assert(v_accept.size() == v_hists.size());
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


class SelectedSelHists: public Hists {
public:
    explicit SelectedSelHists(Context & ctx,
                         const string & dir,
                         const SelItemsHelper & sel_helper,
                         const vector<string> & selections,
                         const vector<string> & veto_selections = {}):
        Hists(ctx, dir),
        h_sel_res(ctx.get_handle<vector<bool>>(sel_helper.s_vec_bool())),
        v_all_items(sel_helper.get_item_names()),
        v_selections(selections),
        v_veto_selections(veto_selections)
    {
        sel_helper.fill_hists_vector(v_hists, dir);
    }

    void insert_hist_and_sel(unsigned pos, Hists * hists, const string & sel) {
        v_hists.insert(v_hists.begin() + pos, move(unique_ptr<Hists>(hists)));
        v_all_items.insert(v_all_items.begin() + pos, sel);
    }

    void insert_additional_hist(Hists * hists) {
        v_hists.push_back(move(unique_ptr<Hists>(hists)));
    }

    virtual void fill(const Event & event) override {
        const auto & v_accept = event.get(h_sel_res);
        assert(v_accept.size() == v_all_items.size());
        for (unsigned i=0; i<v_hists.size(); ++i) {
            bool accept_sels = true;
            for (unsigned j=0; j<v_accept.size(); ++j) {
                for (const string & sel : v_selections) {
                    if (sel == v_all_items[j] && !v_accept[j]) {
                        accept_sels = false;
                        break;    
                    }
                }

                bool ignore_sel = v_veto_selections.size() ? false : true;
                for (const string & sel : v_veto_selections) {
                    if (sel == v_all_items[j]) {
                        ignore_sel = true;   
                    }
                }
                if (!ignore_sel) {
                    if (!v_accept[j]) {
                        accept_sels = false;
                        break;
                    }
                }
            }
            if (accept_sels) {
                v_hists[i]->fill(event);
            }
        }
    }

private:
    Event::Handle<vector<bool>> h_sel_res;
    vector<unique_ptr<Hists>> v_hists;
    vector<string> v_all_items;
    const vector<string> v_selections;
    const vector<string> v_veto_selections;
};


class VLQ2HTCutflow: public Hists {
public:
    explicit VLQ2HTCutflow(Context & ctx,
                  const string & dir,
                  const SelItemsHelper & sel_helper):
        Hists(ctx, dir),
        h_sel(ctx.get_handle<vector<bool>>(sel_helper.s_vec_bool())),
        v_names(sel_helper.get_item_names()),
        h(book<TH1D>("cutflow", "", 1, 0, -1)),
        init_done(false) {}

    void insert_step(unsigned pos, const string & name) {
        assert(!init_done);  // no changing names after init
        v_names.insert(v_names.begin() + pos, name);
    }

    virtual void fill(const uhh2::Event & e) override {
        if (!init_done) {
            initialize_histo();
        }
        const float w = e.weight;
        const auto & sel = e.get(h_sel);
        assert(sel.size() == v_names.size());
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
    void initialize_histo() {
        h->SetBit(TH1::kCanRebin);
        h->Fill("input", 1e-7);
        for (const string & name : v_names) {
            h->Fill(name.c_str(), 1e-7);
        }
        init_done = true;
    }

    Event::Handle<vector<bool>> h_sel;
    vector<string> v_names;
    TH1D * h;
    bool init_done;
};
