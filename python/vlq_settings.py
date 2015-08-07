import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = kError;')

from varial import settings


settings.box_text_size = 0.03
settings.colors = {
    'MC_TTbar': 632,
    'MC_WJets': 596,
    'MC_ZJets': 870,
    'MC_T': 434,
}
# 840, 902, 797, 800, 891, 401, 800,
# 838, 420, 403, 893, 881, 804, 599, 615, 831, 403, 593, 872

settings.defaults_Legend.update({
    'x_pos': 0.85,
    'y_pos': 0.5,
    'label_width': 0.28,
    'label_height': 0.04,
    'opt': 'f',
    'opt_data': 'p',
    'reverse': True
})

settings.canvas_size_x = 550
settings.canvas_size_y = 400
settings.root_style.SetPadRightMargin(0.3)
settings.rootfile_postfixes = ['.root', '.png']

#settings.max_open_root_files = 100
#settings.max_num_processes = 20