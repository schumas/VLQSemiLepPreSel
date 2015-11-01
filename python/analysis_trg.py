#!/usr/bin/env python

from varial.extensions.sframe import SFrame
import varial.extensions.make
import varial.tools
import os

import plot_trg


dir_name = 'VLQSLPS'
uhh_base = os.getenv('CMSSW_BASE') + '/src/UHH2/'


tc = varial.tools.ToolChain(
    dir_name,
    [
        varial.extensions.make.Make([
            uhh_base + 'core',
            uhh_base + 'common',
            uhh_base + 'VLQSemiLepPreSel',
        ]),
        SFrame(uhh_base + 'VLQSemiLepPreSel/config/VLQTrigStudy.xml'),
        varial.tools.ToolChainParallel(
            'Plots', 
            lazy_eval_tools_func=lambda: 
                plot_trg.mk_plttrs('%s/../SFrame/*.root' % varial.analysis.cwd)
        ),
        varial.tools.WebCreator(),
    ]
)


# varial.settings.max_num_processes = 1
# varial.settings.rootfile_postfixes += ['.pdf']
# varial.tools.Runner(tc, True)
import varial.main
varial.main.main(toolchain=tc, try_reuse_results=True)
