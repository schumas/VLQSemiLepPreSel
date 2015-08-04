#pragma once

#include <map>
#include <memory>
#include <set>
#include <fstream>
#include <sstream>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/VLQSemiLepPreSel/include/HandleSelection.h"
#include "UHH2/VLQSemiLepPreSel/include/HandleHist.h"
#include "UHH2/VLQSemiLepPreSel/include/HandleToBranchWriter.h"


using namespace std;
using namespace uhh2;


class SelectionItem {
public:
    SelectionItem(const string & name) : name_(name) {}
    virtual Selection       * make_selection(Context & ctx) const = 0;
    virtual Hists           * make_hists(Context & ctx, const string & dir) const = 0;
    virtual AnalysisModule  * make_branch_writer(Context & ctx, TTree * tree) const = 0;
    virtual void              declare_for_output(Context & ctx) const = 0;
    const string & name() const {return name_;}
    virtual float cutvalue_min() const = 0;
    virtual float cutvalue_max() const = 0;

protected:
    const string name_;
};


template<typename DATATYPE>
class SelectionItemData : public SelectionItem {
public:
    SelectionItemData(const string & name, const string & title,
                      int n_bins, float x_min, float x_max,
                      DATATYPE min_val=-99999.0, DATATYPE max_val=99999.0):
        SelectionItem(name), title_(title),
        n_bins_(n_bins), x_min_(x_min), x_max_(x_max),
        min_val_(min_val), max_val_(max_val) {}

    virtual Selection * make_selection(Context & ctx) const override {
        return new HandleSelection<DATATYPE>(ctx, name_, min_val_, max_val_);
    }

    virtual Hists * make_hists(Context & ctx, const string & dir) const override {
        // TODO pass title with _%g(min)_to_%g(max)
        return new HandleHist<DATATYPE>(ctx, dir, name_, title_.c_str(), n_bins_, x_min_, x_max_);
    }

    virtual AnalysisModule * make_branch_writer(Context & ctx, TTree * tree) const override {
        return new HandleToBranchWriter<DATATYPE>(ctx, name_, tree);
    }

    virtual void declare_for_output(Context & ctx) const override {
        ctx.declare_event_output<DATATYPE>(name_);
    }

    virtual float cutvalue_min() const override {
        return min_val_;
    }

    virtual float cutvalue_max() const override {
        return max_val_;
    }

private:
    string title_;
    int n_bins_;
    float x_min_;
    float x_max_;
    DATATYPE min_val_;
    DATATYPE max_val_;
};


class SelItemsHelper {
public:
    SelItemsHelper(const vector<shared_ptr<SelectionItem>> & sel_items,
                   Context & context,
                   const vector<string> & names = vector<string>()):
        items(sel_items),
        ctx(context),
        item_names(names.size() ? names : all_item_names())
    {
        set<string> name_set;
        for (const auto & it : items) {
            string n = it->name();
            if (name_set.count(n)) {
                throw runtime_error(
                    "SelectionItem available more then once: " + n);
            }
            name_set.insert(n);
        }
    }

    const vector<string> & get_item_names() const {
        return item_names;
    }

    const vector<string> & all_item_names() const {
        static vector<string> v;
        if (!v.size()) {
            for (const auto & it: items) {
                v.push_back(it->name());
            }
        }
        return v;
    }

    const shared_ptr<SelectionItem> get_sel_item(const string & name) const {
        static map<string, shared_ptr<SelectionItem>> m;
        if (!m.size()) {
            for (const auto & it: items) {
                m[it->name()] = it;
            }
        }
        if (m.count(name)) {
            return m[name];
        } else {
            return NULL;
        }
    }

    void fill_hists_vector(vector<unique_ptr<Hists>> & target,
                          const string & dir) const {
        for (const auto & name: item_names) {
            target.emplace_back(get_sel_item(name)->make_hists(ctx, dir));
        }
    }

    void fill_sel_vector(vector<unique_ptr<Selection>> & target) const {
        for (const auto & name: item_names) {
            target.emplace_back(get_sel_item(name)->make_selection(ctx));
        }
    }

    void fill_wrtr_vector(vector<unique_ptr<AnalysisModule>> & target,
                          TTree * tree) const {
        for (const auto & name: item_names) {
            target.emplace_back(get_sel_item(name)->make_branch_writer(ctx, tree));
        }
    }

    TreeWriter * make_tree_writer(const string & filename) const {
        auto * wrtr = new TreeWriter(ctx, filename);
        fill_wrtr_vector(wrtr->writers(), wrtr->tree());
        return wrtr;
    }

    void declare_items_for_output() const {
        for (const auto & name: item_names) {
            get_sel_item(name)->declare_for_output(ctx);
        }
    }

    void write_cuts_to_texfile(const string & filename = "cuts.tex") const {
        ofstream ofs(filename);
        for (const auto & it : items) {
            float min = it->cutvalue_min();
            float max = it->cutvalue_max();
            ofs << it->name() << " & ";
            if (min > -99998.) {
                ofs << min << " & ";
            } else {
                ofs << "--- & ";
            }
            if (max < 99998.) {
                ofs << max;
            } else {
                ofs << "---";
            }
            ofs << " \\\\ \n";
        }
        ofs.close();
    }

private:
    const vector<shared_ptr<SelectionItem>> & items;
    Context & ctx;
    const vector<string> & item_names;
};


class SelectionProducer : public AnalysisModule {
public:
    explicit SelectionProducer(Context & ctx,
                               const SelItemsHelper & sel_helper):
        h_sel_res(ctx.get_handle<vector<bool>>("sel_accept")),
        h_all_acc(ctx.get_handle<bool>("sel_all_accepted"))
    {
        sel_helper.fill_sel_vector(v_sel);
    }

    void insert_selection(unsigned pos, Selection * sel) {
        v_sel.insert(v_sel.begin() + pos, move(unique_ptr<Selection>(sel)));
    }

    virtual bool process(Event & event) override {
        bool all_accepted = true;
        vector<bool> v_accept(v_sel.size());
        for (unsigned i=0; i<v_sel.size(); ++i) {
            bool accept = v_sel[i]->passes(event);
            v_accept[i] = accept;
            if (!accept) {
                all_accepted = false;
            }
        }
        event.set(h_sel_res, v_accept);
        event.set(h_all_acc, all_accepted);
        return all_accepted;
    }

private:
    vector<unique_ptr<Selection>> v_sel;
    Event::Handle<vector<bool>> h_sel_res;
    Event::Handle<bool> h_all_acc;
};