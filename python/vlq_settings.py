import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = kError;')

from varial import settings


settings.box_text_size = 0.03
settings.colors = {
    'TTbar': 632 - 7,
    'WJets': 901,
    'ZJets': 856,
    'DYJets': 856,
    'DYJetsToLL': 856,
    'SingleT': 433,
    'SingleTop': 433,
    'QCD': 851,

#     'TTbar': 800 + 8,
#     'WJets': 856,
#     'ZJets': 851,
#     'DYJets': 851,
#     'DYJetsToLL': 851,
#     'SingleT': 800 - 3,
#     'QCD': 433,

    # Heiners colors
    'TpB_TH_0700': ROOT.kGreen + 1,
    'TpB_TH_1200': ROOT.kBlue,
    'TpB_TH_1700': ROOT.kRed,
    'Signal_TpB_TH_LH_M0700': ROOT.kGreen + 1,
    'Signal_TpB_TH_LH_M1200': ROOT.kBlue,
    'Signal_TpB_TH_LH_M1700': ROOT.kRed,

    'Bkg': ROOT.kYellow - 6,
    'BkgPostFit': ROOT.kSpring + 6,
}
# 596, , 797, 800, 891, 401, 800,
# 838, 420, 403, 893, 881, 804, 599, 615, 831, 403, 593, 872

settings.defaults_Legend.update({
    'x_pos': 0.81,
    'y_pos': 0.5,
    'label_width': 0.2,
    'label_height': 0.07,
    'opt': 'f',
    'opt_data': 'p',
    'reverse': True
})

settings.stacking_order = [
    'TTbar',
    'WJets',
    'SingleT',
    'DYJets',
    'QCD',
]

settings.box_text_size = 0.05
settings.canvas_size_x = 550
settings.canvas_size_y = 400
settings.root_style.SetPadRightMargin(0.3)
settings.rootfile_postfixes = ['.root', '.png']

#settings.max_open_root_files = 100
#settings.max_num_processes = 20
