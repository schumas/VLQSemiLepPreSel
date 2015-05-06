#pragma once

#include <assert.h>
#include <typeinfo>
#include <vector>
#include "TFile.h"
#include "TTree.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"


using namespace std;
using namespace uhh2;


template<typename DATATYPE>
class HandleToBranchWriter : public AnalysisModule {
public:
    HandleToBranchWriter(Context & ctx, const string & handlename, TTree * tree):
        value(0.),
        hndl(ctx.get_handle<DATATYPE>(handlename))
    {
        auto type_idx = type_index(typeid(DATATYPE));
        string typ = type_idx == type_index(typeid(int))    ? "/I" :
                     type_idx == type_index(typeid(float))  ? "/F" :
                     type_idx == type_index(typeid(double)) ? "/D" :
                     "";
        assert(typ.size());

        tree->Branch(handlename.c_str(),
                     &value,
                     (handlename + typ).c_str());
    }

    virtual bool process(Event & e) override {
        value = e.get(hndl);
        return true;
    }

protected:
    DATATYPE value;
    Event::Handle<DATATYPE> hndl;
};


class WeightToBranchWriter : public HandleToBranchWriter<double> {
public:
    WeightToBranchWriter(Context & ctx, TTree * tree):
        HandleToBranchWriter(ctx, "weight", tree) {}

    virtual bool process(Event & e) override {
        value = e.weight;
        return true;
    }
};


class TreeWriter : public AnalysisModule {
public:
    TreeWriter(Context & ctx,
               const string & filename):
        context(ctx),
        out_file(unique_ptr<TFile>(new TFile(("tree_writer."+filename+".root").c_str(),
                                             "RECREATE",
                                             filename.c_str()))),
        out_tree(unique_ptr<TTree>(new TTree(filename.c_str(),
                                             filename.c_str())))
    {
        v_writers.emplace_back(new WeightToBranchWriter(context, out_tree.get()));
    }

    ~TreeWriter() {
        out_file->Write();
        out_file->Close();
        out_tree.release();
    }

    TTree * tree() {
        return out_tree.get();
    }

    vector<unique_ptr<AnalysisModule>> & writers() {
        return v_writers;
    }

    virtual bool process(Event & e) override {
        for (auto &wrtr : v_writers) {
            wrtr->process(e);
        }
        out_tree->Fill();
        return true;
    }

private:
    Context & context;
    unique_ptr<TFile> out_file;
    unique_ptr<TTree> out_tree;
    vector<unique_ptr<AnalysisModule>> v_writers;
};