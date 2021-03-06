#pragma once

#include "UHH2/common/include/ObjectIdUtils.h"

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"


/**
 *   Example class for booking and filling histograms, the new version using AnalysisModule mechanisms.
 */

class CustomizableGenHists: public uhh2::Hists {
public:

    struct GenHistColl
    {
        int pdgid;
        unsigned int order_num;
        TH1F *h_pt, *h_eta, *h_phi, *h_n, *h_mass, *h_charge, *h_decay, *h_mother, *h_dRDecay, *h_dPhiDecay, *h_dEtaDecay;
        boost::optional<GenParticleId> genp_id;
    };

    // use the same constructor arguments as Hists for forwarding:
    CustomizableGenHists(uhh2::Context & ctx, const std::string & dirname, const std::string & h_part_ht = "parton_ht");
    virtual void fill(const uhh2::Event & ev) override;
    void add_genhistcoll(int pdgid, unsigned int order_num, const std::vector<std::string> & variables = {"number", "mass", "charge", "decay", "mother", "dRDecay", "dPhiDecay", "dEtaDecay"}, const boost::optional<GenParticleId> & genp_id = boost::none,
                const std::string & suffix = "");

    void fill_genhistcoll(const uhh2::Event & ev, GenHistColl & gen_histcoll);

    std::map<int, std::pair<float, float> > & minmax_pts() {
        return minmax_pts_;
    }
    std::map<int, std::pair<float, float> > & minmax_masses() {
        return minmax_masses_;
    }

    virtual ~CustomizableGenHists();
private:
    template<class T>
    void fill_wrapper(std::vector<const T*> plot_particles, const uhh2::Event & event,
        GenHistColl & gen_histcoll);

    template<class T>
    void fill_hists(const T * ipart, const std::vector<GenParticle> & genparticles,
        GenHistColl & gen_histcoll, double w);

    uhh2::Context & ctx_;
    std::string dirname_;
	uhh2::Event::Handle<double> h_part_ht_;

    static std::map<int, std::pair<float, float> > minmax_pts_;
    static std::map<int, std::pair<float, float> > minmax_masses_;

    std::vector<GenHistColl> all_hists_;

	TH1F *spec_parton_ht, *spec_deltaR_bb_h, *spec_deltaR_bb_min, *spec_max_deltaR_topprod;
	TH2F *spec_top_pt_vs_max_dR;
};
