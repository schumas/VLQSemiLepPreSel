#pragma once

#include <algorithm>
#include <cmath>

#include "TH2.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Selection.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Utils.h"

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
        if (e.is_valid(h_prim_lep) && e.get(h_prim_lep).pt() > 0.001) {
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
                          string const & h_name = "ST",
                          boost::optional<JetId> const & jet_id = boost::none):
        jet_id_(jet_id),
        h_st(ctx.get_handle<double>(h_name)),
        h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")) {}

    virtual bool process(Event & event) override {
        if (!event.is_valid(h_primlep)) {
            return false;
        }
        float st = event.get(h_primlep).pt();
        st += event.met->pt();
        for (const auto & j : *event.jets) {
            if (jet_id_){
                if ((*jet_id_)(j, event)) {
                    st += j.pt();
                }
            } else {
                st += j.pt();
            }
        }
        event.set(h_st, st);
        return true;
    }

private:
    boost::optional<JetId> jet_id_;
    Event::Handle<double> h_st;
    Event::Handle<FlavorParticle> h_primlep;
};  // class STCalculator


class LepPtPlusMETProducer: public AnalysisModule {
public:
    explicit LepPtPlusMETProducer(
        Context & ctx,
        const string & prim_lep_name = "PrimaryLepton",
        const string & out_name = "lep_plus_met",
        const string & out_name_vec_sum = "lep_plus_met_vec_sum"
    ):
        h_primlep(ctx.get_handle<FlavorParticle>(prim_lep_name)),
        h_out(ctx.get_handle<float>(out_name)),
        h_out_vec_sum(ctx.get_handle<float>(out_name_vec_sum)) {}

    virtual bool process(Event & e) override {
        const auto & prim_lep = e.get(h_primlep);
        e.set(
            h_out, 
            prim_lep.pt() + e.met->pt()
        );
        e.set(
            h_out_vec_sum, 
            (prim_lep.v4() + e.met->v4()).pt()
        );
        return true;
    }

private:
    Event::Handle<FlavorParticle> h_primlep;
    Event::Handle<float> h_out;
    Event::Handle<float> h_out_vec_sum;
};


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

    explicit TriggerAcceptProducer(Context & ctx,
                                   const vector<string> & trig_names,
                                   const vector<string> & trig_veto_names,
                                   const string & h_name = "trigger_accept"):
        v_trig_names(trig_names),
        v_trig_veto_names(trig_veto_names),
        h(ctx.get_handle<int>(h_name)) {}

    virtual bool process(Event & e) override {
        // process vetos first
        if (v_trig_veto_names.size()) {
            vector<Event::TriggerIndex> v_trig_veto_ind;
            for (const auto & trig_name : v_trig_veto_names) {
                v_trig_veto_ind.emplace_back(e.get_trigger_index(trig_name));
            }
            for (Event::TriggerIndex & trig_ind : v_trig_veto_ind) {
                if (e.passes_trigger(trig_ind)) {
                    e.set(h, 0);
                    return false;
                }
            }
        }

        // check for normal triggers
        vector<Event::TriggerIndex> v_trig_ind;
        for (const auto & trig_name : v_trig_names) {
            v_trig_ind.emplace_back(e.get_trigger_index(trig_name));
        }
        for (Event::TriggerIndex & trig_ind : v_trig_ind) {
            if (e.passes_trigger(trig_ind)) {
                e.set(h, 1);
                return true;
            }
        }

        e.set(h, 0);
        return false;
    }

private:
    vector<string> v_trig_names;
    vector<string> v_trig_veto_names;
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


/**
  * veto_mother_id is stronger then mother_id and will always be checked for
  * all mothers.
  */
class GenParticleMotherId
{
public:
    GenParticleMotherId(int mother_id = 0, int veto_mother_id = 0) :
        mother_id_({mother_id}), veto_mother_id_({veto_mother_id})
        {}

    GenParticleMotherId(const vector<int> & mother_id = {0},
                        const vector<int> & veto_mother_id = {0}) :
        mother_id_(mother_id), veto_mother_id_(veto_mother_id)
        {}

    bool operator()(const GenParticle & genp, const Event & event) const
    {
        if (mother_id_.size() == 1 && mother_id_[0] == 0
            && veto_mother_id_.size() == 1 && veto_mother_id_[0] == 0) {
            return true;
        }

        bool right_mother = mother_id_.size() == 1 && mother_id_[0] == 0;
        GenParticle const * gen_mother = findMother(genp, event.genparticles);
        while (gen_mother)
        {
            for (int mom_id : mother_id_) {
                if (mom_id > 0 && abs(gen_mother->pdgId()) == mom_id) {
                    right_mother = true;
                }
            }
            for (int veto_id : veto_mother_id_) {
                if (veto_id > 0 && abs(gen_mother->pdgId()) == veto_id) {
                    return false;
                }
            }
            gen_mother = findMother(*gen_mother, event.genparticles);
        }
        return right_mother;
    }

private:
    vector<int> mother_id_, veto_mother_id_;
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


class GenParticlePdgIdId
{
public:
    GenParticlePdgIdId(const std::vector<int> & pdgids) :
        pdgids_(pdgids) {}

    bool operator()(const GenParticle & genp, const Event &) const
    {
        for (int id : pdgids_) {
            if (genp.pdgId() == id)  // ATTENTION: NO ABS! NEED TO SPECIFY ALL
            {
                return true;
            }
        }

        // cout << "  Found right mother!\n";
        return false;
    }

private:
    std::vector<int> pdgids_;
};  // GenParticlePdgIdId


class LeptonicDecayVLQ : public AnalysisModule {
public:
    explicit LeptonicDecayVLQ(const vector<int> &lep_ids = {11, -11, 13, -13}):
        id_(AndId<GenParticle>(
            GenParticlePdgIdId(lep_ids),
            GenParticleMotherId({8000001, 6000006, 6000007, 6000008})
        )) {}

    virtual bool process(Event & event) override {
        for (auto & gp : *event.genparticles) {
            if (id_(gp, event)) {
                return true;
            }
        }
        return false;
    }

private:
    GenParticleId id_;
};


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


class LeadingTopjetLorentzVectorProducer : public AnalysisModule {
public:
    explicit LeadingTopjetLorentzVectorProducer(Context & ctx,
                        const string & h_in,
                        const string & h_out):
        h_in_(ctx.get_handle<vector<TopJet>>(h_in)),
        h_out_(ctx.get_handle<LorentzVector>(h_out)) {}

    virtual bool process(Event & event) override {
        if (event.is_valid(h_in_)) {
            vector<TopJet> & coll = event.get(h_in_);
            if (coll.size()) {
                event.set(h_out_, coll[0].v4());
            }
            return true;
        }
        return false;
    }
private:
    Event::Handle<vector<TopJet>> h_in_;
    Event::Handle<LorentzVector> h_out_;
};


class LeadingTopjetMassProducer : public AnalysisModule {
public:
    explicit LeadingTopjetMassProducer(Context & ctx,
                        const string & h_in,
                        const string & h_out):
        h_in_(ctx.get_handle<vector<TopJet>>(h_in)),
        h_out_(ctx.get_handle<float>(h_out)) {}

    virtual bool process(Event & event) override {
        if (event.is_valid(h_in_)) {
            vector<TopJet> & coll = event.get(h_in_);
            if (coll.size()) {
                if (coll[0].subjets().size()){
                    LorentzVector sum_subjets;
                    for (Jet const & subjet : coll[0].subjets())
                        sum_subjets += subjet.v4();
                    event.set(h_out_, sum_subjets.M());
                } else {
                    event.set(h_out_, -1.);    
                }
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
    Event::Handle<vector<TopJet>> h_in_;
    Event::Handle<float> h_out_;
};


class LeadingTopjetNSubjettinessProducer : public AnalysisModule {
public:
    explicit LeadingTopjetNSubjettinessProducer(Context & ctx,
                        const string & h_in,
                        const string & h_out,
                        bool tau21_not_tau32 = true):
        h_in_(ctx.get_handle<vector<TopJet>>(h_in)),
        h_out_(ctx.get_handle<float>(h_out)),
        tau21_not_tau32_(tau21_not_tau32) {}

    virtual bool process(Event & event) override {
        if (event.is_valid(h_in_)) {
            vector<TopJet> & coll = event.get(h_in_);
            if (coll.size()) {
                if (tau21_not_tau32_)
                    event.set(h_out_, coll[0].tau2()/coll[0].tau1());
                else
                    event.set(h_out_, coll[0].tau3()/coll[0].tau2());
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
    Event::Handle<vector<TopJet>> h_in_;
    Event::Handle<float> h_out_;
    bool tau21_not_tau32_;
};


template<typename TYPE>
class PartPtProducer: public AnalysisModule {
public:
    explicit PartPtProducer(Context & ctx,
                            const string & h_in,
                            const string & h_out,
                            int part_num = -1):
        h_in_(ctx.get_handle<vector<TYPE>>(h_in)),
        h_out_(ctx.get_handle<float>(h_out)),
        part_num_(part_num) {}

    virtual bool process(Event & e) override {
        if (e.is_valid(h_in_)){
            const vector<TYPE> & coll = e.get(h_in_);
            if (part_num_ < 0) {
                if (coll.size() > 0) {
                    e.set(h_out_, coll.back().pt());
                    return true;
                } else {
                    e.set(h_out_, -1.);
                    return false;
                }
            } else if (part_num_ > 0){
                if (int(coll.size()) >= part_num_) {
                    e.set(h_out_, coll[part_num_-1].pt());
                    return true;
                } else {
                    e.set(h_out_, -1.);
                    return false;
                }
            } else {
                std::cout << "In PartPtProducer: to calculate pt of the pt leading particle, give 1 as argument!\n";
                assert(false);
            }
        }

        e.set(h_out_, -1.);
        return false;

    }

private:
    Event::Handle<vector<TYPE>> h_in_;
    Event::Handle<float> h_out_;
    int part_num_;
};  // PartPtProducer


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


class PrimaryLeptonDeltaPhiId {
public :
    explicit PrimaryLeptonDeltaPhiId(
        Context & ctx, 
        float min_dist, 
        const string & prim_lep="PrimaryLepton"
    ):
        min_dist_(min_dist),
        prim_lep_(ctx.get_handle<FlavorParticle>(prim_lep)) {}

    bool operator()(const Particle & p, const Event & e) const {
        const auto & lep = e.get(prim_lep_);
        return deltaPhi(lep, p) > min_dist_;
    }

private:
    float min_dist_;
    Event::Handle<FlavorParticle> prim_lep_;
};  // class PrimaryLeptonDeltaPhiId


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
};  // VectorAndSelection


template<typename TYPE>
class MinMaxDeltaRId
{
public:
    MinMaxDeltaRId(Context & ctx,
                const string & h_comp_coll,
                float min_dr = 1.0,
                bool only_leading = false,
                bool use_min = true) :
        h_comp_coll_(ctx.get_handle<vector<TYPE>>(h_comp_coll)),
        min_dr_(min_dr),
        only_leading_(only_leading),
        use_min_(use_min),
        use_handle_(false)
        {}

    MinMaxDeltaRId(Context & ctx,
                const string & h_comp_coll,
                const string & h_min_dr,
                bool only_leading = false,
                bool use_min = true) :
        h_comp_coll_(ctx.get_handle<vector<TYPE>>(h_comp_coll)),
        h_min_dr_(ctx.get_handle<float>(h_min_dr)),
        only_leading_(only_leading),
        use_min_(use_min),
        use_handle_(true)
        {}

    bool operator()(const Particle & part, const Event & event) const
    {
        // TODO: make assert statement that part (or rather, TYPE1) really inherits from Particle!
        if (event.is_valid(h_comp_coll_)){
            float dyn_mindr = 1.0;
            if (use_handle_) {
                if (event.is_valid(h_min_dr_))
                    dyn_mindr = event.get(h_min_dr_);
            } else {
                dyn_mindr = min_dr_;
            }
            const vector<TYPE> & comp_coll = event.get(h_comp_coll_);
            if (only_leading_){
                if (comp_coll.size()){
                    const TYPE & ld_part = comp_coll[0];
                    if (use_min_ ? deltaR(part, ld_part) <= dyn_mindr : deltaR(part, ld_part) > dyn_mindr)
                        return false;
                }
                return true;
            } else {
                const TYPE * closest_part = closestParticle<TYPE>(part, comp_coll);
                if (closest_part){
                    if (use_min_ ? deltaR(part, *closest_part) <= dyn_mindr : deltaR(part, *closest_part) > dyn_mindr)
                        return false;
                }
                return true;
            }
        }

        std::cout << "WARNING: in MinMaxDeltaRId: handle to h_comp_coll_ is not valid!\n";
        return true;
    }

private:
    Event::Handle<vector<TYPE>> h_comp_coll_;
    Event::Handle<float> h_min_dr_;
    float min_dr_;
    bool only_leading_;
    bool use_min_;
    bool use_handle_;
};  // MinMaxDeltaRId


class TwoDCutSel: public Selection {
public:
    explicit TwoDCutSel(Context & ctx,
                        float min_dr,
                        float min_ptrel,
                        const string & dr_name = "TwoDCut_dr",
                        const string & pt_name = "TwoDCut_ptrel"):
        h_dr(ctx.get_handle<float>(dr_name)),
        h_pt(ctx.get_handle<float>(pt_name)),
        min_dr_(min_dr),
        min_ptrel_(min_ptrel) {}

    virtual bool passes(const Event & e) override {
        if (!e.is_valid(h_dr) || !e.is_valid(h_pt)) {
            return false;
        }
        return (e.get(h_dr) > min_dr_) || (e.get(h_pt) > min_ptrel_);
    }

private:
    Event::Handle<float> h_dr;
    Event::Handle<float> h_pt;
    float min_dr_;
    float min_ptrel_;
};


class TwoDCutHist: public Hists {
public:
    explicit TwoDCutHist(Context & ctx,
                         const string & dirname,
                         const string & dr_name = "TwoDCut_dr",
                         const string & pt_name = "TwoDCut_ptrel"):
        Hists(ctx, dirname),
        h_dr(ctx.get_handle<float>(dr_name)),
        h_pt(ctx.get_handle<float>(pt_name)),
        hist(book<TH2F>("TwoDCut",
                        ";min #DeltaR(lep., jet);min p_{T,rel}(lep., jet)",
                        200, 0., 1., 200, 0., 500.)) {}

    virtual void fill(const Event & e) override {
        if (!e.is_valid(h_dr) || !e.is_valid(h_pt)) {
            return;
        }
        hist->Fill(e.get(h_dr), e.get(h_pt), e.weight);
    }

private:
    Event::Handle<float> h_dr;
    Event::Handle<float> h_pt;
    TH2F * hist;
};


class AntiHiggsBVetoTag {
public:
    explicit AntiHiggsBVetoTag(float minmass = 60.f, float maxmass = std::numeric_limits<float>::infinity(), 
                               JetId const & id = CSVBTag(CSVBTag::WP_MEDIUM)) :
        minmass_(minmass), maxmass_(maxmass), btagid_(id) {}

    bool operator()(TopJet const & topjet, uhh2::Event const & event) const {
        auto subjets = topjet.subjets();
        if(subjets.size() < 2) return false;
        unsigned n_sj_btags = 0;
        unsigned n_sj_btagvetos = 0;
        for (const auto & sj : subjets) {
            if (btagid_(sj, event))
                n_sj_btags++;
            else if (!btagid_(sj, event))
                n_sj_btagvetos++;
        }

        if (!(n_sj_btags == 1 && n_sj_btagvetos >= 1))
            return false;

        LorentzVector firsttwosubjets = subjets[0].v4() + subjets[1].v4();
        if(!firsttwosubjets.isTimelike()) {
            return false;
        }
        auto mjet = firsttwosubjets.M();
        if(mjet < minmass_) return false;
        if(mjet > maxmass_) return false;
        return true;
    }

private:
    float minmass_, maxmass_;
    JetId btagid_;

};


class METProducer: public AnalysisModule {
public:
    explicit METProducer(Context & ctx,
                                  const string & h_name = "met"):
        h(ctx.get_handle<float>(h_name)) {}

    virtual bool process(Event & e) override {
        if (e.met) {
            e.set(h, e.met->pt());
            return true;
        } else {
            e.set(h, 0.);
            return false;
        }
    }

private:
    Event::Handle<float> h;
};  // METProducer


class NInputEventsHist: public Hists {
public:
    explicit NInputEventsHist(Context & ctx):
        Hists(ctx, ""),
        n_events_total_str(ctx.get("n_events_total")),
        hist(book<TH1F>("NInputEventsHist", "NInputEventsHist", 1, -.5, 0.5))
    {
        hist->SetBit(TH1::kCanRebin);
    }

    virtual void fill(const Event &) override {
        // the number of total events should not be added when adding multiple
        // outputs from proof workers. This is why it's used as axis label.
        hist->Fill(n_events_total_str.c_str(), 1.);  
    }

private:
    string n_events_total_str;
    TH1F * hist;
};


} // namespace
