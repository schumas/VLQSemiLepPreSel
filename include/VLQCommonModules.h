#include <algorithm>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"


using namespace std;
using namespace uhh2;


class FwdJetSwitch: public AnalysisModule {
public:
    explicit FwdJetSwitch(Context & ctx):
        hndl(ctx.get_handle<std::vector<Jet> >("fwd_jets")) {}

    bool process(Event & event){
        std::vector<Jet> fwd;
        std::vector<Jet> cnt;
        for(const auto & jet: *event.jets) {
            if (fabs(jet.eta()) > 2.4) {
                fwd.push_back(jet);
            } else {
                cnt.push_back(jet);
            }
        }
        event.set(hndl, fwd);
        swap(*event.jets, cnt);
        return true;
    }

private:
    Event::Handle<std::vector<Jet> > hndl;
};  // class FwdJetSwitch


class NBTagProducer: public AnalysisModule {
public:
    explicit NBTagProducer(Context & ctx,
                           CSVBTag::wp wp = CSVBTag::WP_LOOSE):
        hndl(ctx.get_handle<int>("n_btags")),
        tagger(CSVBTag(wp)) {}

    bool process(Event & event){
        int nbtag = 0;
        for(const Jet & j : *event.jets){
            if (tagger(j, event)) {
                ++nbtag;
            }
        }
        event.set(hndl, nbtag);
        return true;
    }

private:
    Event::Handle<int> hndl;
    CSVBTag tagger;
};  // class NBTagProducer


class NLeadingBTagProducer: public AnalysisModule {
public:
    explicit NLeadingBTagProducer(Context & ctx,
                                  CSVBTag::wp wp = CSVBTag::WP_LOOSE):
        hndl(ctx.get_handle<int>("n_leading_btags")),
        tagger(CSVBTag(wp)) {}

    bool process(Event & event){
        int ntag = 0;
        for(const Jet & j : *event.jets){
            if (tagger(j, event)) {
                ++ntag;
            } else {
                break;
            }
        }
        event.set(hndl, ntag);
        return true;
    }

private:
    Event::Handle<int> hndl;
    CSVBTag tagger;
};  // class NLeadingBTagProducer


class NHTagProducer: public AnalysisModule {
public:
    explicit NHTagProducer(Context & ctx, const std::string & coll):
        hndl(ctx.get_handle<int>("n_higgs_tags")),
        coll_hndl(ctx.get_handle<std::vector<TopJet>>(coll)),
        tagger(HiggsTag()) {}

    bool process(Event & event){
        int ntag = 0;
        for(const TopJet & j : event.get(coll_hndl)){
            if (tagger(j, event)) {
                ++ntag;
            }
        }
        event.set(hndl, ntag);
        return true;
    }

private:
    Event::Handle<int> hndl;
    Event::Handle<std::vector<TopJet>> coll_hndl;
    TopJetId tagger;
};  // class NHTagProducer


class STCalculator: public uhh2::AnalysisModule {
public:
    explicit STCalculator(uhh2::Context & ctx):
        h_st(ctx.get_handle<double>("ST")),
        h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")) {}

    virtual bool process(uhh2::Event & event) override {
        if (!event.is_valid(h_primlep)) {
            return false;
        }
        float st = event.get(h_primlep).pt();
        st += event.met->pt();
        for (const auto & j : *event.jets) {
            if (fabs(j.eta()) < 2.4) {
                st += j.pt();
            }
        }
        event.set(h_st, st);
        return true;
    }

private:
    uhh2::Event::Handle<double> h_st;
    uhh2::Event::Handle<FlavorParticle> h_primlep;
};  // class STCalculator