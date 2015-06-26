#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "TH1F.h"


using namespace std;
using namespace uhh2;


template<class HTYPE>
class HandleHist: public Hists {
public:
    template<class... TARGS>
    explicit HandleHist(Context & ctx,
                        const string & dirname,
                        const string & handlename,
                        TARGS... args):
        Hists(ctx, dirname),
        hndl(ctx.get_handle<HTYPE>(handlename)),
        hist(book<TH1F>(handlename.c_str(), args...)) {}

    virtual void fill(const Event & e) override {
        if (e.is_valid(hndl)) {
            hist->Fill(e.get(hndl), e.weight);
        }
    }

private:
    Event::Handle<HTYPE> hndl;
    TH1F * hist;
};


template<class HTYPE, class CHI2_HTYPE>
class HandleHistChi2Weight: public Hists {
public:
    template<class... TARGS>
    explicit HandleHistChi2Weight(Context & ctx,
                                  const string & dirname,
                                  const string & handlename,
                                  const string & chi2_handlename,
                                  TARGS... args):
        Hists(ctx, dirname),
        hndl(ctx.get_handle<HTYPE>(handlename)),
        chi2_hndl(ctx.get_handle<CHI2_HTYPE>(chi2_handlename)),
        hist(book<TH1F>((handlename + "_chi2_weighted").c_str(), args...)),
        weight_sum(0.),
        n_events(0.) {}

    virtual ~HandleHistChi2Weight() {
        float average_weight = weight_sum / n_events;
        hist->Scale(1./average_weight);
    }

    virtual void fill(const Event & e) override {
        if (e.is_valid(hndl) && e.is_valid(chi2_hndl)) {
            float chi2_weight = 1 / (1 + e.get(chi2_hndl));
            weight_sum += chi2_weight;
            n_events += 1.;
            hist->Fill(
                e.get(hndl),
                e.weight * chi2_weight
            );
        }
    }

private:
    Event::Handle<HTYPE> hndl;
    Event::Handle<CHI2_HTYPE> chi2_hndl;
    TH1F * hist;
    float weight_sum;
    float n_events;
};