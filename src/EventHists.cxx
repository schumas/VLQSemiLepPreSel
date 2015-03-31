#include "UHH2/VLQSemiLepPreSel/include/EventHists.h"

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"

#include "UHH2/VLQSemiLepPreSel/include/VLQCommonModules.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

HistCollector::HistCollector(Context & ctx, const string & dirname, bool gen_plots, JetId const & btag_id) :
    Hists(ctx, dirname),
    el_hists(new ExtendedElectronHists(ctx, dirname+"/ElectronHists", gen_plots)),
    mu_hists(new ExtendedMuonHists(ctx, dirname+"/MuonHists", gen_plots)),
    tau_hists(new TauHists(ctx, dirname+"/TauHists")),
    ev_hists(new ExtendedEventHists(ctx, dirname+"/EventHists")),
    jet_hists(new ExtendedJetHists(ctx, dirname+"/JetHists")),
    cmstopjet_hists(new ExtendedTopJetHists(ctx, dirname+"/CMSTopJetHists", btag_id, 4)),
    heptopjet_hists(new ExtendedTopJetHists(ctx, dirname+"/HEPTopJetHists", btag_id, 4, "patJetsHepTopTagCHSPacked")),
    ca8prunedtopjet_hists(new ExtendedTopJetHists(ctx, dirname+"/CA8PrunedTopJetHists", btag_id, 4, "patJetsCa8CHSJetsPrunedPacked")),
    ca15filteredtopjet_hists(new ExtendedTopJetHists(ctx, dirname+"/CA15FilteredTopJetHists", btag_id, 4, "patJetsCa15CHSJetsFilteredPacked")),
    gen_hists(gen_plots ? new CustomizableGenHists(ctx, dirname+"/GenHists", "parton_ht") : NULL)
    {
        if (gen_hists)
        {
            gen_hists->add_genhistcoll(8, 1, true, true, true, true);
            gen_hists->add_genhistcoll(8, 2, true, true, true, true);
            gen_hists->add_genhistcoll(6, 1, true, true, true, true);
            gen_hists->add_genhistcoll(6, 2, true, true, true, true);
            gen_hists->add_genhistcoll(25, 1, true, true, true, true);
            gen_hists->add_genhistcoll(25, 2, true, true, true, true);
            gen_hists->add_genhistcoll(0, 0);
            gen_hists->add_genhistcoll(11, 1, false, true, false, false, GenParticleId(GenParticleMotherId(6)), "_from_top");
            gen_hists->add_genhistcoll(11, 2, false, true, false, false, GenParticleId(GenParticleMotherId(6)), "_from_top");
            gen_hists->add_genhistcoll(13, 1, false, true, false, false, GenParticleId(GenParticleMotherId(6)), "_from_top");
            gen_hists->add_genhistcoll(13, 2, false, true, false, false, GenParticleId(GenParticleMotherId(6)), "_from_top");
            gen_hists->add_genhistcoll(11, 1, false, true, false, false,
                GenParticleId(
                    AndId<GenParticle>(
                        GenParticleMotherId(24,6),
                        GenParticleMotherId(24,25)
                    )
                ), "_from_tpW");
            gen_hists->add_genhistcoll(11, 2, false, true, false, false, GenParticleId(AndId<GenParticle>(GenParticleMotherId(24,6),GenParticleMotherId(24,25))), "_from_tpW");
            gen_hists->add_genhistcoll(13, 1, false, true, false, false, GenParticleId(AndId<GenParticle>(GenParticleMotherId(24,6),GenParticleMotherId(24,25))), "_from_tpW");
            gen_hists->add_genhistcoll(13, 2, false, true, false, false, GenParticleId(AndId<GenParticle>(GenParticleMotherId(24,6),GenParticleMotherId(24,25))), "_from_tpW");

        }
    } 

void HistCollector::fill(const Event & event) {
    el_hists->fill(event);
    mu_hists->fill(event);
    tau_hists->fill(event);
    ev_hists->fill(event);
    jet_hists->fill(event);
    cmstopjet_hists->fill(event);
    heptopjet_hists->fill(event);
    ca8prunedtopjet_hists->fill(event);
    ca15filteredtopjet_hists->fill(event);
    if (gen_hists) gen_hists->fill(event);
}

HistCollector::~HistCollector(){
    delete el_hists;
    delete mu_hists;
    delete tau_hists;
    delete ev_hists;
    delete jet_hists;
    delete cmstopjet_hists;
    delete heptopjet_hists;
    delete ca8prunedtopjet_hists;
    delete ca15filteredtopjet_hists;
    delete gen_hists;
}