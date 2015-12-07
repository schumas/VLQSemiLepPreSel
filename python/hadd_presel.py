#!/usr/bin/env python

from varial.extensions.hadd import Hadd
import varial.tools
import sys

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Please give me exactly one input pattern!'
        exit(-1)

    input_pat = sys.argv[1]
    hadd = Hadd(
        input_pat, 
        [
            'uhh2.AnalysisModuleRunner.DATA.Run2015D', 
            'uhh2.AnalysisModuleRunner.MC.TTbar',
            'uhh2.AnalysisModuleRunner.MC.WJets',
            'uhh2.AnalysisModuleRunner.MC.DYJets',
            'uhh2.AnalysisModuleRunner.MC.QCD',
            'uhh2.AnalysisModuleRunner.MC.SingleT',
            'uhh2.AnalysisModuleRunner.MC.TpB_TH_LH_M0700_lepDecay',
            'uhh2.AnalysisModuleRunner.MC.TpB_TH_LH_M1700_lepDecay',
        ], 
        add_aliases_to_analysis=False,
    )
    varial.tools.Runner(hadd)


