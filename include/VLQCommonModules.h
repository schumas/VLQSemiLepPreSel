#pragma once

#include <algorithm>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/GenTools.h"

using namespace std;
using namespace uhh2;

namespace {

class FwdJetSwitch: public AnalysisModule {
public:
    explicit FwdJetSwitch(Context & ctx):
        hndl(ctx.get_handle<vector<Jet> >("fwd_jets")) {}

    bool process(Event & event) override {
        vector<Jet> fwd;
        vector<Jet> cnt;
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
    Event::Handle<vector<Jet> > hndl;
};  // class FwdJetSwitch


class BJetsProducer: public AnalysisModule {
public:
    explicit BJetsProducer(Context & ctx,
                           CSVBTag::wp wp = CSVBTag::WP_MEDIUM):
        hndl(ctx.get_handle<vector<Jet>>("b_jets")),
        tagger(CSVBTag(wp)) {}

    bool process(Event & event) override {
        vector<Jet> b_jets;
        for(const Jet & j : *event.jets){
            if (tagger(j, event)) {
                b_jets.push_back(j);
            }
        }
        event.set(hndl, b_jets);
        return true;
    }

private:
    Event::Handle<vector<Jet>> hndl;
    CSVBTag tagger;
};  // class BJetsProducer


class NBTagProducer: public AnalysisModule {
public:
    explicit NBTagProducer(Context & ctx,
                           CSVBTag::wp wp = CSVBTag::WP_MEDIUM):
        hndl(ctx.get_handle<int>("n_btags")),
        tagger(CSVBTag(wp)) {}

    bool process(Event & event) override {
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
                                  CSVBTag::wp wp = CSVBTag::WP_MEDIUM):
        hndl(ctx.get_handle<int>("n_leading_btags")),
        tagger(CSVBTag(wp)) {}

    bool process(Event & event) override {
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


class LeadingJetPtProducer: public AnalysisModule {
public:
    explicit LeadingJetPtProducer(Context & ctx):
        h(ctx.get_handle<float>("leading_jet_pt")) {}

    virtual bool process(Event & e) override {
        if (e.jets->size() > 0) {
            e.set(h, e.jets->at(0).pt());
            return true;
        } else {
            e.set(h, 0.);
            return false;
        }
    }

private:
    Event::Handle<float> h;
};


class LeptonPtProducer: public AnalysisModule {
public:
    explicit LeptonPtProducer(Context & ctx):
        h(ctx.get_handle<float>("primary_lepton_pt")),
        h_prim_lep(ctx.get_handle<FlavorParticle>("PrimaryLepton")) {}

    virtual bool process(Event & e) override {
        if (e.is_valid(h_prim_lep)) {
            e.set(h, e.get(h_prim_lep).pt());
        } else {
            e.set(h, 0.);
        }
        return true;
    }

private:
    Event::Handle<float> h;
    Event::Handle<FlavorParticle> h_prim_lep;
};


class LargestJetEtaProducer: public AnalysisModule {
public:
    explicit LargestJetEtaProducer(Context & ctx,
                           const string & jets_name = "jets",
                           const string & fwd_jets_name = ""):
        h_largest_jet_eta(ctx.get_handle<float>("largest_jet_eta")),
        h_jets(ctx.get_handle<vector<Jet>>(jets_name)),
        h_fwd_jets(ctx.get_handle<vector<Jet>>(fwd_jets_name)) {}

    virtual float largest_eta(const vector<Jet> & jets) {
        float largest_jet_eta = 0.;
        for (const Jet & j : jets) {
            float eta = j.eta();
            if (fabs(largest_jet_eta) < fabs(eta)) {
                largest_jet_eta = eta;
            }
        }
        return largest_jet_eta;
    }

    virtual bool process(Event & e) override {
        if (e.is_valid(h_fwd_jets) && e.get(h_fwd_jets).size()) {
            e.set(h_largest_jet_eta, largest_eta(e.get(h_fwd_jets)));
            return true;
        }
        if (e.is_valid(h_jets) && e.get(h_jets).size()) {
            e.set(h_largest_jet_eta, largest_eta(e.get(h_jets)));
            return true;
        }
        return false;
    }

private:
    Event::Handle<float>        h_largest_jet_eta;
    Event::Handle<vector<Jet>>  h_jets;
    Event::Handle<vector<Jet>>  h_fwd_jets;
};  // class LargestJetEtaProducer


class STCalculator: public AnalysisModule {
public:
    explicit STCalculator(Context & ctx):
        h_st(ctx.get_handle<double>("ST")),
        h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")) {}

    virtual bool process(Event & event) override {
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
    Event::Handle<double> h_st;
    Event::Handle<FlavorParticle> h_primlep;
};  // class STCalculator


class JetPtSorter : public AnalysisModule {
public:
    explicit JetPtSorter() {}
    virtual bool process(Event & event) override {
        vector<Jet> & ev_jets = *event.jets;
        sort_by_pt(ev_jets);

        return true;
    }
};  // class JetPtSorter


class JetTagCalculator : public AnalysisModule {
public:
    explicit JetTagCalculator(Context & ctx, string hndl_name, JetId const & id = JetId(CSVBTag(CSVBTag::WP_MEDIUM))) :
        tagger_(id), hndl_(ctx.get_handle<int>(hndl_name)) {}

    virtual bool process(Event & event) override {
        int n_btags = 0;
        for (const Jet & jet : *event.jets) {
            if (tagger_(jet, event))
                n_btags++;
        }
        event.set(hndl_, n_btags);
        return true;
    }

private:
    JetId tagger_;
    Event::Handle<int> hndl_;
};  // class JetTagCalculator


class TaggedTopJetProducer: public AnalysisModule {
public:
    explicit TaggedTopJetProducer(Context & ctx, TopJetId const & id, const string & coll_out, const string & coll_in = ""):
        hndl_in(ctx.get_handle<vector<TopJet>>(coll_in)),
        hndl_out(ctx.get_handle<vector<TopJet>>(coll_out)),
        tagger(id) {}

    bool process(Event & event) override {
        const vector<TopJet> & topjets = event.is_valid(hndl_in) ? event.get(hndl_in) : *event.topjets;
        vector<TopJet> out_jets;
        for(const TopJet & j : topjets) {
            if (tagger(j, event)) {
                out_jets.push_back(j);
            }
        }
        event.set(hndl_out, out_jets);
        return true;
    }

private:
    Event::Handle<vector<TopJet>> hndl_in;
    Event::Handle<vector<TopJet>> hndl_out;
    TopJetId tagger;
};  // class TaggedTopJetProducer


class NTaggedTopJetProducer: public AnalysisModule {
public:
    explicit NTaggedTopJetProducer(Context & ctx, TopJetId const & id, const string hndl_name, const string & coll_in = ""):
        hndl_in(ctx.get_handle<vector<TopJet>>(coll_in)),
        hndl_out(ctx.get_handle<int>(hndl_name)),
        tagger(id) {}

    bool process(Event & event) override {
        const vector<TopJet> & topjets = event.is_valid(hndl_in) ? event.get(hndl_in) : *event.topjets;
        int n_topjettags = 0;
        for(const TopJet & j : topjets) {
            if (tagger(j, event)) {
                n_topjettags++;
            }
        }
        event.set(hndl_out, n_topjettags);
        return true;
    }

private:
    Event::Handle<vector<TopJet>> hndl_in;
    Event::Handle<int> hndl_out;
    TopJetId tagger;
};  // class NTaggedTopJetProducer


class EventWeightOutputHandle: public AnalysisModule {
public:
    explicit EventWeightOutputHandle(Context & ctx):
        hndl(ctx.declare_event_output<double>("weight")) {}

    bool process(Event & e) override {
        e.set(hndl, e.weight);
        return true;
    }

private:
    Event::Handle<double> hndl;
};


class GenParticleMotherId
{
public:
    GenParticleMotherId(int mother_id = 0, int veto_mother_id = 0) :
        mother_id_(mother_id), veto_mother_id_(veto_mother_id)
        {}

    bool operator()(const GenParticle & genp, const Event & event) const
    {
        if (mother_id_ > 0 || veto_mother_id_ > 0)
        {
            // cout << "  Looking for particle with mother " << mother_id_ << " and not from " << veto_mother_id_ << endl;
            bool right_mother = mother_id_ > 0 ? false : true;
            GenParticle const * gen_mother = findMother(genp, event.genparticles);
            while (gen_mother)
            {
                // cout << "   Mother id: " << gen_mother->pdgId() << endl;
                if (mother_id_ > 0 && abs(gen_mother->pdgId()) == mother_id_)
                {
                    right_mother = true;
                }
                else if (veto_mother_id_ > 0 && abs(gen_mother->pdgId()) == veto_mother_id_)
                {
                    right_mother = false;
                    break;
                }
                gen_mother = findMother(*gen_mother, event.genparticles);
            }
            if (!right_mother)
            {
                // cout << "  Bad mother, rejected!\n";
                return false;
            }
        }

        // cout << "  Found right mother!\n";
        return true;
    }

private:
    int mother_id_, veto_mother_id_;
};


class GenParticleDaughterId
{
public:
    GenParticleDaughterId(int part_id, int daughter1_id = 0, int daughter2_id = 0) :
        part_id_(part_id), daughter1_id_(daughter1_id), daughter2_id_(daughter2_id)
        {}

    bool operator()(const GenParticle & genp, const Event & event) const
    {
        if (abs(genp.pdgId()) == part_id_)
        {
            GenParticle const * daughter1 = genp.daughter(event.genparticles, 1);
            GenParticle const * daughter2 = genp.daughter(event.genparticles, 2);
            if (!(daughter1 && daughter2))
                return false;
            if ((abs(daughter1->pdgId()) == daughter1_id_ && abs(daughter2->pdgId()) == daughter2_id_)
                || (abs(daughter1->pdgId()) == daughter2_id_ && abs(daughter2->pdgId()) == daughter1_id_))
                return true;
        }

        // cout << "  Found right mother!\n";
        return false;
    }

private:
    int part_id_, daughter1_id_, daughter2_id_;
};


// function to grab the best hypothesis
// note: member HYP.discriminators must be public map<string, float>
template<typename HYP>
const HYP * get_best_hypothesis(
    const vector<HYP> & hyps,
    const string & label,
    float & best_discr)
{
    const HYP * best = nullptr;
    float current_best_disc = numeric_limits<float>::infinity();
    for(const auto & hyp : hyps){
        if(!hyp.discriminators.count(label)) continue;
        float disc = hyp.discriminators.find(label)->second;
        if(disc < current_best_disc){
            best = &hyp;
            current_best_disc = disc;
        }
    }
    best_discr = current_best_disc;
    return best;  // note: might be nullptr
}

}