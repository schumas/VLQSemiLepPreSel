#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "TH1F.h"


using namespace std;
using namespace uhh2;


template<class HANDLETYPE>
class HandleHist: public Hists {
public:
    template<class... TARGS>
    explicit HandleHist(Context & ctx,
                        const string & dirname,
                        const string & handlename,
                        TARGS... args):
        Hists(ctx, dirname),
        hndl(ctx.get_handle<HANDLETYPE>(handlename)) {
            hist=book<TH1F>(handlename.c_str(), args...);
        }

    virtual void fill(const Event & e) override {
        if (e.is_valid(hndl)) {
            hist->Fill(e.get(hndl), e.weight);
        }
    }

private:
    Event::Handle<HANDLETYPE> hndl;
    TH1F * hist;
};