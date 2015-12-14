#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/TauHists.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/LuminosityHists.h"


#include "UHH2/VLQSemiLepPreSel/include/CustomizableGenHists.h"

#include "TH1F.h"


using namespace uhh2;
using namespace std;

class ExtendedEventHists : public EventHists {
public:
    ExtendedEventHists(uhh2::Context & ctx, const std::string & dirname,
        const std::string & h_btags_loose = "n_btags",
        const std::string & h_btags_medium = "",
        const std::string & h_btags_tight = ""
        ) :
        EventHists(ctx, dirname),
        h_btags_loose_(ctx.get_handle<int>(h_btags_loose)),
        h_btags_medium_(ctx.get_handle<int>(h_btags_medium)),
        h_btags_tight_(ctx.get_handle<int>(h_btags_tight))
        // h_toptags_(ctx.get_handle<int>("n_toptags")),
        // h_higgstags_(ctx.get_handle<int>("n_higgstags"))
        {
            Weights = book<TH1F>("Weights_own", "weights", 2000,0,200);
            MET = book<TH1F>("MET_own", "missing E_{T}", 50,0,1000);
            HTLep = book<TH1F>("HTLep_own", "H_{T} Lep", 50, 0, 1000);
            h_n_btags_loose = book<TH1F>("jets_Nbs_loose", "N_{b-tags loose}", 11, -.5, 10.5);
            h_n_btags_medium = book<TH1F>("jets_Nbs_medium", "N_{b-tags medium}", 11, -.5, 10.5);
            h_n_btags_tight = book<TH1F>("jets_Nbs_tight", "N_{b-tags tight}", 11, -.5, 10.5);
            // h_n_toptags = book<TH1F>("jets_Ntops", "N_{cms top tags}", 20, 0, 20);
            // h_n_higgstags = book<TH1F>("jets_Nhiggs", "N_{higgs tags}", 20, 0, 20);

            h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
        }

    virtual void fill(const uhh2::Event & event) override {
        EventHists::fill(event);
        double w = event.weight;
        int n_btags_loose = event.is_valid(h_btags_loose_) ? event.get(h_btags_loose_) : 0;
        int n_btags_medium = event.is_valid(h_btags_medium_) ? event.get(h_btags_medium_) : 0;
        int n_btags_tight = event.is_valid(h_btags_tight_) ? event.get(h_btags_tight_) : 0;
        // int n_toptags = event.is_valid(h_toptags_) ? event.get(h_toptags_) : 0;
        // int n_higgstags = event.is_valid(h_higgstags_) ? event.get(h_higgstags_) : 0;        

        h_n_btags_loose->Fill(n_btags_loose, w);
        h_n_btags_medium->Fill(n_btags_medium, w);
        h_n_btags_tight->Fill(n_btags_tight, w);
        // h_n_toptags->Fill(n_toptags, w);
        // h_n_higgstags->Fill(n_higgstags, w);

    }

private:
    TH1F *h_n_btags_loose, *h_n_btags_medium, *h_n_btags_tight; // , *h_n_toptags, *h_n_higgstags
    uhh2::Event::Handle<int> h_btags_loose_, h_btags_medium_, h_btags_tight_; // , h_toptags_, h_higgstags_
};

// new el hists class based on the ElectronHists class in the common package but changing the histogram ranges for some of the histograms

class ExtendedElectronHists : public ElectronHists {
public:
    ExtendedElectronHists(uhh2::Context & ctx, const std::string & dirname, bool gen_plots = true) :
        ElectronHists(ctx, dirname, gen_plots)
        {
            isolation   = book<TH1F>("isolation_own",   "relIso electron",          20,0,4);
            isolation_1 = book<TH1F>("isolation_1_own", "relIso electron 1",        20,0,4);
            isolation_2 = book<TH1F>("isolation_2_own", "relIso electron 2",        20,0,4);
        }

};

// new muon hists class based on the MuonHists class in the common package but changing the histogram ranges for some of the histograms

class ExtendedMuonHists : public MuonHists {
public:
    ExtendedMuonHists(uhh2::Context & ctx, const std::string & dirname, bool gen_plots = true) :
        MuonHists(ctx, dirname, gen_plots)
        {
            isolation   = book<TH1F>("isolation_own",   "relIso electron",          20,0,4);
            isolation_1 = book<TH1F>("isolation_1_own", "relIso electron 1",        20,0,4);
            isolation_2 = book<TH1F>("isolation_2_own", "relIso electron 2",        20,0,4);
            // deltaR_ak8 = book<TH1F>("deltaR_ak8_cleaned", "dR closest Ak8 jet", 40, 0., 2.);
            // deltaR_ak8_uncleaned = book<TH1F>("deltaR_ak8_uncleaned", "dR closest Ak8 jet", 40, 0., 2.);
            // h_ak8 = ctx.get_handle<vector<TopJet>>("topjets");
            // h_ak8_uncleaned = ctx.get_handle<vector<TopJet>>("ak8jets_uncleaned");
        }
    // virtual void fill(const uhh2::Event & event) override {
    //     MuonHists::fill(event);
    //     if (event.muons->size()) {
    //         const Muon & muon = (*event.muons)[0];
    //         if (event.is_valid(h_ak8)) {
    //             vector<TopJet> const & ak8_jets = event.get(h_ak8);
    //             auto nj = closestParticle(muon, ak8_jets);
    //             auto drmin_val = nj ? deltaR(muon, *nj) : numeric_limits<float>::infinity();
    //             deltaR_ak8->Fill(drmin_val, event.weight);
    //         }
    //         if (event.is_valid(h_ak8_uncleaned)) {
    //             vector<TopJet> const & ak8_jets = event.get(h_ak8_uncleaned);
    //             auto nj = closestParticle(muon, ak8_jets);
    //             auto drmin_val = nj ? deltaR(muon, *nj) : numeric_limits<float>::infinity();
    //             deltaR_ak8_uncleaned->Fill(drmin_val, event.weight);
    //         }
    //     }
    // }

protected:
    TH1F *deltaR_ak8, *deltaR_ak8_uncleaned;
    Event::Handle<std::vector<TopJet>> h_ak8, h_ak8_uncleaned;

};

// new jet hists class based on the JetHists class in the common package; add jet histos w/o JEC

class ExtendedJetHists : public JetHists {
public:
    ExtendedJetHists(uhh2::Context & ctx, const std::string & dirname, const unsigned int NumberOfPlottedJets=4, const std::string & collection = "") :
        JetHists(ctx, dirname, NumberOfPlottedJets, collection)
    {
        alljets_nojec = book_jetHist("jet_nojec","_jet_nojec",20,1500);

        std::vector<double> minPt {20,20,20,20};
        std::vector<double> maxPt {1500,1000,500,350};
        std::vector<std::string> axis_suffix {"first jet no jec","second jet no jec","third jet no jec","fourth jet no jec"};

        for(unsigned int i =0; i<NumberOfPlottedJets; i++){
            if(i<4){
                single_jetHists_nojec.push_back(book_jetHist(axis_suffix[i],"_"+to_string(i+1)+"_nojec",minPt[i],maxPt[i]));
            }
            else {
                single_jetHists_nojec.push_back(book_jetHist(to_string(i+1)+"-th jet","_"+to_string(i+1)+"_nojec",20,500));
            }
        }

    }

    virtual void fill(const uhh2::Event & event) override
    {
        JetHists::fill(event);

        auto w = event.weight;
        const auto jets = collection.empty() ? event.jets : &event.get(h_jets);
        assert(jets);
        for(unsigned int i = 0; i <jets->size(); i++){
            auto jet = (*jets)[i];
            jet.set_v4(jet.v4() * jet.JEC_factor_raw());
            fill_jetHist(jet,alljets_nojec,w);
            if(i < single_jetHists_nojec.size()){
                fill_jetHist(jet, single_jetHists_nojec[i], w);
            }
        }

    }

protected:
    jetHist alljets_nojec;
    std::vector<jetHist> single_jetHists_nojec;

};

// new jet hists class based on the MuonHists class in the common package; add jet histos w/o JEC

class ExtendedTopJetHists : public TopJetHists {
public:

    struct tag_variables_hists
    {
        // for higgs-tagging
        TH1F *higgs_n_subjet_btags, *higgs_mass_two_leading_subjets;
        // for cmstoptagging
        TH1F *cmstoptag_mass_allsubjets, *cmstoptag_min_mass_disubj;
    };
    ExtendedTopJetHists(uhh2::Context & ctx, const std::string & dirname, const JetId & b_tag = CSVBTag(CSVBTag::WP_MEDIUM), const unsigned int NumberOfPlottedJets=4, const std::string & collection = "") :
        TopJetHists(ctx, dirname, NumberOfPlottedJets, collection), b_tag_(b_tag)
    {
        alljets_tagvars = book_tagvarHist("all topjets","");

        string axis_suffix = "topjet";
        vector<string> axis_singleSubjetSuffix {"first ","second ","third ","fourth "};

        for(unsigned int i =0; i<NumberOfPlottedJets; i++){
            if(i<4){
                single_tagvars.push_back(book_tagvarHist(axis_singleSubjetSuffix[i]+axis_suffix,string("_")+to_string(i+1)));
            }
            else{
                single_tagvars.push_back(book_tagvarHist(to_string(i+1)+string("-th jet"),string("_")+to_string(i+1)));
            }
        }

    }

    tag_variables_hists book_tagvarHist(const std::string & axisSuffix, const std::string & histSuffix)
    {
        tag_variables_hists tag_vars;
        tag_vars.higgs_n_subjet_btags = book<TH1F>("higgs_n_subjet_btags"+histSuffix, "number b-tagged subjets "+axisSuffix, 10, -.5, 9.5);
        tag_vars.higgs_mass_two_leading_subjets = book<TH1F>("higgs_mass_two_leading_subjets"+histSuffix, "mass of two leading b-tagged subjets "+axisSuffix, 100, 0., 500.);
        tag_vars.cmstoptag_mass_allsubjets = book<TH1F>("cmstoptag_mass_allsubjets"+histSuffix, "mass of all subjets "+axisSuffix, 100, 0., 500.);
        tag_vars.cmstoptag_min_mass_disubj = book<TH1F>("cmstoptag_min_mass_disubj"+histSuffix, "min mass of two subjets "+axisSuffix, 100, 0., 500.);
        return tag_vars;

    }

    void fill_tagvarHist(const uhh2::Event & event, const TopJet & topjet, tag_variables_hists & tag_vars, double weight)
    {
        auto subjets = topjet.subjets();
        if(subjets.size() >= 2)
        {
            clean_collection(subjets, event, b_tag_);
            tag_vars.higgs_n_subjet_btags->Fill(subjets.size(), weight);
            if (subjets.size() >= 2)
            {
                sort_by_pt(subjets);

                LorentzVector firsttwosubjets = subjets[0].v4() + subjets[1].v4();
                if(firsttwosubjets.isTimelike()) {
                    auto mjet = firsttwosubjets.M();
                    tag_vars.higgs_mass_two_leading_subjets->Fill(mjet, weight);
                }
            }
        }

        if(topjet.subjets().size() >= 3)
        {
            LorentzVector allsubjets;
            for(const auto & subjet : topjet.subjets()) {
                allsubjets += subjet.v4();
            }
            if(allsubjets.isTimelike())
            {
                auto mjet = allsubjets.M();
                tag_vars.cmstoptag_mass_allsubjets->Fill(mjet, weight);
                auto mmin = m_disubjet_min(topjet);
                tag_vars.cmstoptag_min_mass_disubj->Fill(mmin, weight);
            }
        }
    }

    virtual void fill(const uhh2::Event & event) override
    {
        if (!collection.empty() && !event.is_valid(h_topjets)) {
            return;
        }

        TopJetHists::fill(event);

        auto w = event.weight;
        const auto jets = collection.empty() ? event.topjets : &event.get(h_topjets);
        assert(jets);
        for(unsigned int i = 0; i <jets->size(); i++){
            const auto & jet = (*jets)[i];
            fill_tagvarHist(event,jet,alljets_tagvars,w);
            if(i < single_tagvars.size()){
                fill_tagvarHist(event,jet,single_tagvars[i],w);
            }
        }

    }

protected:
    JetId b_tag_;
    tag_variables_hists alljets_tagvars;
    std::vector<tag_variables_hists> single_tagvars;
};


class HistCollector : public uhh2::Hists {
public:
    HistCollector(uhh2::Context & ctx, const std::string & dirname, bool gen_plots = true, JetId const & btag_id = CSVBTag(CSVBTag::WP_MEDIUM));
    virtual ~HistCollector();

    virtual void fill(const uhh2::Event & ev) override;

private:
    ExtendedElectronHists * el_hists;
    ExtendedMuonHists * mu_hists;
    TauHists * tau_hists;
    ExtendedEventHists * ev_hists;
    ExtendedJetHists * jet_hists;
    ExtendedTopJetHists * cmstopjet_hists;
    ExtendedTopJetHists * heptopjet_hists;
    ExtendedTopJetHists * ca8prunedtopjet_hists;
    ExtendedTopJetHists * ca15filteredtopjet_hists;
    CustomizableGenHists * gen_hists;
    LuminosityHists * lumi_hist;
};