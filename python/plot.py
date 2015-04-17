#!/usr/bin/env python

# YOU NEED VARIAL FOR THIS! https://github.com/HeinAtCERN/Varial

import common

import time
import varial.tools
import varial.generators as gen

varial.settings.max_open_root_files = 100
varial.settings.defaults_Legend.update({
    'x_pos': 0.85,
    'y_pos': 0.5,
    'label_width': 0.28,
    'label_height': 0.04,
    'opt': 'f',
    'opt_data': 'p',
    'reverse': True
})
varial.settings.canvas_size_x = 600
varial.settings.canvas_size_y = 500
varial.settings.root_style.SetPadRightMargin(0.3)
varial.settings.rootfile_postfixes = ['.root', '.png']

varial.settings.max_num_processes = 20

dir_input = '/nfs/dust/cms/user/tholenhe/VLQSemiLepPreSel/PHYS14-ntuple2-v2/'


def only_leading_n_obj(wrps, n=50):
    i = 0
    for w in wrps:
        i += 1
        if i < n:
            yield w
        else:
            break


def apply_match_eff(wrps):
    factors = {
        # QCD
        '_HT250to500': 0.1685,
        '_HT500to1000': 0.2103,
        '_HT1000ToInf': 0.2358,

        # WJets
        '_LNu_HT100to200_20x25': 0.096,
        '_LNu_HT200to400_20x25': 0.084,
        '_LNu_HT400to600_20x25': 0.075,
        '_LNu_HT600toInf_20x25': 0.063,

        # ZJets
        '_LL_HT100to200_20x25': 0.099,
        '_LL_HT200to400_20x25': 0.081,
        '_LL_HT400to600_20x25': 0.067,
        '_LL_HT600toInf_20x25': 0.062,
    }
    for w in wrps:
        for k in factors:
            if k in w.file_path:
                w.lumi /= factors[k]
        w = varial.op.norm_to_lumi(w)
        yield w


def merge_samples(wrps):
    wrps = common.merge_decay_channels(wrps, (
        '_LNu_HT100to200_20x25',
        '_LNu_HT200to400_20x25',
        '_LNu_HT400to600_20x25',
        '_LNu_HT600toInf_20x25',
    ))
    wrps = common.merge_decay_channels(wrps, (
        '_LL_HT100to200_20x25',
        '_LL_HT200to400_20x25',
        '_LL_HT400to600_20x25',
        '_LL_HT600toInf_20x25',
    ))
    wrps = common.merge_decay_channels(wrps, (
        '_HT250to500',
        '_HT500to1000',
        '_HT1000ToInf',
    ))
    return wrps


def loader_hook(wrps):
    #wrps = only_leading_n_obj(wrps, 20)
    #wrps = apply_match_eff(wrps)
    wrps = common.add_wrp_info(wrps)
    wrps = merge_samples(wrps)
    wrps = (w for w in wrps if w.histo.Integral() > 1e-5)
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
    kws['save_lin_log_scale'] = True
    return varial.tools.Plotter(**kws)


def plotter_factory_norm(**kws):
    kws['hook_loaded_histos'] = loader_hook_norm
    kws['save_lin_log_scale'] = True
    return varial.tools.Plotter(**kws)


def plotter_factory_stack(**kws):
    kws['hook_loaded_histos'] = loader_hook
    kws['save_lin_log_scale'] = True
    kws['plot_setup'] = gen.mc_stack_n_data_sum
    return varial.tools.Plotter(**kws)

p1 = varial.tools.mk_rootfile_plotter(
    pattern=dir_input + '*.root',
    name='VLQ_presel_stack',
    plotter_factory=plotter_factory_stack,
    combine_files=True,
)

p2 = varial.tools.mk_rootfile_plotter(
    pattern=dir_input + '*.root',
    name='VLQ_presel_norm',
    plotter_factory=plotter_factory_norm,
    combine_files=True,
)

p3 = varial.tools.mk_rootfile_plotter(
    pattern=dir_input + '*.root',
    name='VLQ_presel_norm_no_signal',
    plotter_factory=plotter_factory,
    combine_files=True,
    filter_keyfunc=lambda w: common.is_signal(w.file_path)
)

tc_inner = varial.tools.ToolChainParallel(
    'VLQ_presel', [
        p1.tool_chain[0],
        p2.tool_chain[0],
        p3.tool_chain[0]
    ]
)
tc = varial.tools.ToolChain(
    'host_toolchain', [tc_inner]
)


if __name__ == '__main__':
    time.sleep(1)

    #import cProfile
    #varial.settings.use_parallel_chains = False
    #cProfile.runctx('p1.run()',globals(),locals(),'prof_varial_plotting.txt')
    #print 'done profiling'

    tc.run()
    varial.tools.WebCreator().run()
    print 'done.'






