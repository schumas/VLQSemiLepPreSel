#!/usr/bin/env python

# YOU NEED VARIAL FOR THIS! https://github.com/HeinAtCERN/Varial

import glob
import time
from os.path import basename

import cutflow_tables
import vlq_settings
import common

import varial.generators as gen
import varial.tools
import varial.plotter
import varial.history


# varial.settings.debug_mode = True
input_pat = './Hadd/*.root'
# input_pat = '*.root'


# don't plot all signal samples
input_pat = glob.glob(input_pat)
input_pat = filter(
    lambda s: ('_TH_' not in s and 'TpTp_' not in s) or '_lepDecay' in s,
    input_pat
)


def fix_mc_norm(wrps):

    @varial.history.track_history
    def fix_norm(w):
        w.lumi *= .5
        return w

    for w in wrps:
        if not w.is_data:
            w = fix_norm(w)
        yield w


def loader_hook(wrps, rebin_max_bins=40):
    if rebin_max_bins:
        wrps = varial.gen.gen_noex_rebin_nbins_max(wrps, rebin_max_bins)

    # wrps = common.yield_n_objs(wrps, 20)
    wrps = common.add_wrp_info(wrps)
    wrps = sorted(wrps, key=lambda w: w.in_file_path + '__' + w.sample)

    #wrps = gen.imap_conditional(
    #    wrps, lambda w: w.name == 'ST', 
    #    gen.op.rebin, bin_bounds=range(0, 1550, 50)
    #)
    #wrps = gen.imap_conditional(
    #    wrps, lambda w: w.name == 'leading_jet_pt', 
    #    gen.op.rebin, bin_bounds=range(0, 760, 20)
    #)
    #wrps = gen.imap_conditional(
    #    wrps, lambda w: w.name == 'primary_lepton_pt', 
    #    gen.op.rebin, bin_bounds=range(0, 520, 20)
    #)

    wrps = common.label_axes(wrps)
    wrps = gen.gen_make_th2_projections(wrps)
    #wrps = gen.gen_make_eff_graphs(wrps)
    return wrps


def loader_hook_norm(wrps):
    wrps = loader_hook(wrps)
    wrps = gen.switch(
        wrps,
        lambda w: 'TH1' in w.type,
        gen.gen_noex_norm_to_integral
    )
    return wrps


def plotter_factory(**kws):
    kws['hook_loaded_histos'] = loader_hook
    # kws['canvas_decorators'] += [rnd.TitleBox(
    #       text='CMS Simulation 20fb^{-1} @ 13TeV')]
    return varial.tools.Plotter(**kws)


def plotter_factory_norm(**kws):
    kws['hook_loaded_histos'] = loader_hook_norm
    return varial.tools.Plotter(**kws)


def plotter_factory_stack(**kws):
    kws['hook_loaded_histos'] = loader_hook
    kws['plot_setup'] = gen.mc_stack_n_data_sum
    return varial.tools.Plotter(**kws)


if __name__ == '__main__':
    tc = varial.tools.ToolChainParallel(
        'VLQ_presel', [
            varial.tools.mk_rootfile_plotter(
                pattern=input_pat,
                name='VLQ_presel_stack',
                plotter_factory=plotter_factory_stack,
                combine_files=True,
                auto_legend=False,
            ),
            varial.tools.mk_rootfile_plotter(
                pattern=input_pat,
                name='VLQ_presel_norm',
                plotter_factory=plotter_factory_norm,
                combine_files=True,
                auto_legend=False,
            ),
            varial.tools.mk_rootfile_plotter(
                pattern=input_pat,
                name='VLQ_presel_norm_no_signal',
                plotter_factory=plotter_factory_norm,
                combine_files=True,
                auto_legend=False,
                filter_keyfunc=lambda w: not common.is_signal(w.file_path)
            ),
            cutflow_tables.mk_cutflow_chain(input_pat, loader_hook),
        ]
    )

    #import cProfile
    #varial.settings.use_parallel_chains = False
    #cProfile.runctx('p1.run()',globals(),locals(),'prof_varial_plotting.txt')
    #print 'done profiling'

    varial.tools.Runner(tc)
    varial.tools.WebCreator().run()
    print 'done.'






