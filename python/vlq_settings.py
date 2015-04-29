from varial import settings


settings.box_text_size = 0.03
#settings.colors = {
#    'TTJets': 632,
#    'WJets': 878,
#    'ZJets': 596,
#    'TpJ_TH_M800_Tlep': 870,
#    'TpJ_TH_M800_NonTlep': 434,
#}


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

settings.max_open_root_files = 100
settings.max_num_processes = 20