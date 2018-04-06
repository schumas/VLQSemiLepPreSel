#pragma once

#include "UHH2/core/include/Hists.h"
#include "TH1F.h"

class VLQSemiLepPreSelHists: public uhh2::Hists {
public:
    VLQSemiLepPreSelHists(uhh2::Context & ctx, const std::string & dir):
        Hists(ctx, dir),
        h_primlep(ctx.get_handle<FlavorParticle>("PrimaryLepton")),
        h_st(ctx.get_handle<double>("ST")),
        lepPt(book<TH1F>(
            "PrimaryLeptonPt",
            ";primary lepton p_{T};events",
            100, 0, 1000
        )),
        muoPt(book<TH1F>(
            "PrimaryMuonPt",
            ";primary muon p_{T};events",
            100, 0, 1000
        )),
        elePt(book<TH1F>(
            "PrimaryElePt",
            ";primary electron p_{T};events",
            100, 0, 1000
        )),
        muoEta(book<TH1F>(
            "PrimaryMuonEta",
            ";primary muon #eta;events",
	    100, -6., 6.
        )),
        eleEta(book<TH1F>(
            "PrimaryEleEta",
            ";primary electron #eta;events",
	    100, -6., 6.
        )),
        leadingJetPt(book<TH1F>(
            "LeadingJetPt",
            ";leading jet p_{T};events",
            100, 0, 2000
        )),
        st(book<TH1F>(
            "ST",
            ";ST;events",
            100, 0, 4000
        )) {}

    virtual void fill(const uhh2::Event & event) override {
        if (!event.jets->size()) {
            return;
        }
        double w = event.weight;
        float lep_pt = event.get(h_primlep).pt();
        lepPt->Fill(lep_pt, w);
        leadingJetPt->Fill(event.jets->at(0).pt(), w);
        st->Fill(event.get(h_st), w);
        if (event.electrons->size()
            && fabs(event.electrons->at(0).pt() - lep_pt) < 1e-40) {
            elePt->Fill(event.electrons->at(0).pt(), w);
	    eleEta->Fill(event.electrons->at(0).eta(), w);
        }
        if (event.muons->size()
            && fabs(event.muons->at(0).pt() - lep_pt) < 1e-40) {
            muoPt->Fill(event.muons->at(0).pt(), w);
            muoEta->Fill(event.muons->at(0).eta(), w);
        }
    }

    Event::Handle<FlavorParticle> h_primlep;
    Event::Handle<double> h_st;
    TH1F * lepPt;
    TH1F * muoPt;
    TH1F * elePt;
    TH1F * muoEta;
    TH1F * eleEta;
    TH1F * leadingJetPt;
    TH1F * st;
};
