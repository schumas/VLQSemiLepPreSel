#pragma once

#include <algorithm>
#include <cmath>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"


#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/GenTools.h"

using namespace std;
using namespace uhh2;

namespace {

template<typename T>
class AbsValueProducer: public AnalysisModule {
public:
    AbsValueProducer(Context & ctx,
                     const string & h_name):
        h_in(ctx.get_handle<T>(h_name)),
        h_out(ctx.get_handle<T>("abs_" + h_name)) {}

    virtual bool process(Event & e) override {
        if (e.is_valid(h_in)) {
            // will only compile for int, float...
            e.set(h_out, abs(e.get(h_in)));
            return true;
        }
        return false;
    }

private:
    Event::Handle<T> h_in;
    Event::Handle<T> h_out;
};

template<typename T>
class DiffValueProducer: public AnalysisModule {
public:
    DiffValueProducer(Context & ctx,
                     const string & h_name1,
                     const string & h_name2):
        h_in1(ctx.get_handle<T>(h_name1)),
        h_in2(ctx.get_handle<T>(h_name2)),
        h_out(ctx.get_handle<T>("diff_" + h_name1 + '_' + h_name2)) {}

    virtual bool process(Event & e) override {
        if (e.is_valid(h_in1) && e.is_valid(h_in2)) {
            // will only compile for int, float...
            e.set(h_out, e.get(h_in1) - e.get(h_in2));
            return true;
        }
        return false;
    }

private:
    Event::Handle<T> h_in1;
    Event::Handle<T> h_in2;
    Event::Handle<T> h_out;
};

class LeptonPtProducer: public AnalysisModule {
public:
    explicit LeptonPtProducer(Context & ctx,
                              const string & prim_lep_hndl = "PrimaryLepton",
                              const string & h_name = "primary_lepton_pt"):
        h(ctx.get_handle<float>(h_name)),
        h_prim_lep(ctx.get_handle<FlavorParticle>(prim_lep_hndl)) {}

    virtual bool process(Event & e) override {
        if (e.is_valid(h_prim_lep)) {
            e.set(h, e.get(h_prim_lep).pt());
        }
        return true;
    }

private:
    Event::Handle<float> h;
    Event::Handle<FlavorParticle> h_prim_lep;
};  // LeptonPtProducer


class NLeptonsProducer: public AnalysisModule {
public:
    explicit NLeptonsProducer(Context & ctx,
                              string const & h_name = "n_leptons"):
        h(ctx.get_handle<int>(h_name)) {}

    virtual bool process(Event & e) override {
        e.set(h, e.electrons->size() + e.muons->size());
        return true;
    }

private:
    Event::Handle<int> h;
};  // NLeptonsProducer



// DEPRECATED: replaced by CollectionProducer class; needs to be tested still!
// like this:
// static bool is_fwd_jet(const Jet & j, const Event &) {return fabs(j.eta()) >= 2.4;}
// static bool is_cntrl_jet(const Jet & j, const Event &) {return fabs(j.eta()) < 2.4;}
// CollectionProducer<Jet>(ctx, is_fwd_jet, "jets", "fwd_jets")
// CollectionProducer<Jet>(ctx, is_cntrl_jet, "jets", "jets")
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


class NLeadingBTagProducer: public AnalysisModule {
public:
    explicit NLeadingBTagProducer(Context & ctx,
                                  CSVBTag::wp wp = CSVBTag::WP_MEDIUM,
                                  const string & h_name = "n_leading_btags"):
        hndl(ctx.get_handle<int>(h_name)),
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
    explicit LeadingJetPtProducer(Context & ctx,
                                  const string & h_name = "leading_jet_pt"):
        h(ctx.get_handle<float>(h_name)) {}

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
};  // LeadingJetPtProducer


class SubleadingJetPtProducer: public AnalysisModule {
public:
    explicit SubleadingJetPtProducer(Context & ctx,
                                     const string & h_name = "subleading_jet_pt"):
        h(ctx.get_handle<float>(h_name)) {}

    virtual bool process(Event & e) override {
        if (e.jets->size() > 1) {
            e.set(h, e.jets->at(1).pt());
            return true;
        } else {
            e.set(h, 0.);
            return false;
        }
    }

private:
    Event::Handle<float> h;
};  // SubleadingJetPtProducer


class LargestJetEtaProducer: public AnalysisModule {
public:
    explicit LargestJetEtaProducer(Context & ctx,
                                   const string & output_name = "largest_jet_eta",
                                   const string & jets_name = "jets"):
        h_largest_jet_eta(ctx.get_handle<float>(output_name)),
        h_jets(ctx.get_handle<vector<Jet>>(jets_name)) {}

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
        if (e.is_valid(h_jets) && e.get(h_jets).size()) {
            e.set(h_largest_jet_eta, largest_eta(e.get(h_jets)));
            return true;
        }
        return false;
    }

private:
    Event::Handle<float>        h_largest_jet_eta;
    Event::Handle<vector<Jet>>  h_jets;
};  // class LargestJetEtaProducer


class STCalculator: public AnalysisModule {
public:
    explicit STCalculator(Context & ctx,
                          string const & h_name = "ST"):
        h_st(ctx.get_handle<double>(h_name)),
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

// DEPRECATED, use PtSorter below instead
class JetPtSorter : public AnalysisModule {
public:
    explicit JetPtSorter() {}
    virtual bool process(Event & event) override {
        vector<Jet> & ev_jets = *event.jets;
        sort_by_pt(ev_jets);

        return true;
    }
};  // class JetPtSorter


template<typename T>
class PtSorter : public AnalysisModule {
public:
    explicit PtSorter(Context & ctx,
                        const string & h_coll):
        h_coll_(ctx.get_handle<vector<T>>(h_coll)) {}
    virtual bool process(Event & event) override {
        if (event.is_valid(h_coll_)) {
            vector<T> & coll = event.get(h_coll_);
            sort_by_pt(coll);

            return true;
        } else {
            return false;
        }
    }
private:
    Event::Handle<vector<T>> h_coll_;
};  // class PtSorter


class TriggerAcceptProducer : public AnalysisModule {
public:
    explicit TriggerAcceptProducer(Context & ctx,
                                   const vector<string> & trig_names,
                                   const string & h_name = "trigger_accept"):
        v_trig_names(trig_names),
        h(ctx.get_handle<int>(h_name)) {}

    virtual bool process(Event & e) override {
        vector<Event::TriggerIndex> v_trig_ind;
        for (const auto & trig_name : v_trig_names) {
            v_trig_ind.emplace_back(e.get_trigger_index(trig_name));
        }
        // auto ele_trig = e.get_trigger_index("HLT_Ele95_CaloIdVT_GsfTrkIdT_v*");
        // auto mu_trig = e.get_trigger_index("HLT_Mu40_v*");
        int passes_any = 0;
        for (Event::TriggerIndex & trig_ind : v_trig_ind) {
            if (e.passes_trigger(trig_ind)) {
                passes_any = 1;
                break;
            }
        }
        e.set(h, passes_any);
        return true;
    }

private:
    vector<string> v_trig_names;
    Event::Handle<int> h;
};  // TriggerAcceptProducer


class TwoDCutProducer: public AnalysisModule {
public:
    explicit TwoDCutProducer(Context & ctx,
                             const string & primlep_name = "PrimaryLepton",
                             const string & dr_name = "TwoDCut_dr",
                             const string & pt_name = "TwoDCut_ptrel"):
        h_dr(ctx.get_handle<float>(dr_name)),
        h_pt(ctx.get_handle<float>(pt_name)),
        h_prim_lep(ctx.get_handle<FlavorParticle>(primlep_name)) {}

    bool process(Event & e) override {
        if (e.is_valid(h_prim_lep)) {
            auto prim_lep = e.get(h_prim_lep);
            float dr, pt;
            std::tie(dr, pt) = drmin_pTrel(prim_lep, *e.jets);
            e.set(h_dr, dr);
            e.set(h_pt, pt);
            return true;
        }
        return false;
    }

private:
    Event::Handle<float> h_dr;
    Event::Handle<float> h_pt;
    Event::Handle<FlavorParticle> h_prim_lep;
};  // TwoDCutProducer


class EventWeightOutputHandle: public AnalysisModule {
public:
    explicit EventWeightOutputHandle(Context & ctx,
                                     const string & h_name = "weight"):
        hndl(ctx.declare_event_output<double>(h_name)) {}

    bool process(Event & e) override {
        e.set(hndl, e.weight);
        return true;
    }

private:
    Event::Handle<double> hndl;
};  // EventWeightOutputHandle


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
};  // GenParticleMotherId


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
};  // GenParticleDaughterId


template<typename T>
class LeadingPartMassProducer : public AnalysisModule {
public:
    explicit LeadingPartMassProducer(Context & ctx,
                        const string & h_in,
                        const string & h_out):
        h_in_(ctx.get_handle<vector<T>>(h_in)),
        h_out_(ctx.get_handle<float>(h_out)) {}
    virtual bool process(Event & event) override {
        if (event.is_valid(h_in_)) {
            vector<T> & coll = event.get(h_in_);
            if (coll.size()) {
                event.set(h_out_, coll[0].v4().M());
            } else {
                event.set(h_out_, -1.);
            }

            return true;
        } else {
            event.set(h_out_, -1.);
            return false;
        }
    }
private:
    Event::Handle<vector<T>> h_in_;
    Event::Handle<float> h_out_;
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


class VectorAndSelection: public Selection {
public:
    explicit VectorAndSelection(const vector<Selection*> &sel_vec) {
        for (const auto & sel : sel_vec) {
            sel_vec_.emplace_back(sel);
        }
    }

    bool passes(const Event & e) override {
        for (const auto & sel : sel_vec_) {
            if (!sel->passes(e)) {
                return false;
            }
        }
        return true;
    }

private:
    vector<unique_ptr<Selection>> sel_vec_;
};


template<typename TYPE>
class MinDeltaRId
{
public:
    MinDeltaRId(Context & ctx,
                const string & h_comp_coll,
                double min_dr = 1.0,
                bool only_leading = false) :h_comp_coll_(ctx.get_handle<vector<TYPE>>(h_comp_coll)),min_dr_(min_dr),only_leading_(only_leading) {}

    bool operator()(const Particle & part, const Event & event) const
    {
        // TODO: make assert statement that part (or rather, TYPE1) really inherits from Particle!
        if (event.is_valid(h_comp_coll_)){
            const vector<TYPE> & comp_coll = event.get(h_comp_coll_);
            if (only_leading_){
                if (comp_coll.size()){
                    const TYPE & ld_part = comp_coll[0];
                    if (deltaR(part, ld_part) < min_dr_)
                        return false;
                }
                return true;
            } else {
                const TYPE * closest_part = closestParticle<TYPE>(part, comp_coll);
                if (closest_part){
                    if (deltaR(part, *closest_part) < min_dr_)
                        return false;
                }
                return true;
            }
        }

        std::cout << "WARNING: in MinDeltaRId: handle to h_comp_coll_ is not valid!\n";
        return true;
    }

private:
    Event::Handle<vector<TYPE>> h_comp_coll_;
    double min_dr_;
    bool only_leading_;
};  // MinDeltaRId

}
