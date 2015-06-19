import varial.tools

# subfolders
p_nm1sel = 'VLQ_presel_stack/Nm1Selection/'
p_cutflow = 'CutflowTools/'

ext = '.png'
target_ext = '.png'

images = {
    'Nm1Sel_vars.tex': (
        p_nm1sel + 'primary_lepton_pt_log' + ext,
        p_nm1sel + 'leading_jet_pt_log' + ext,
        p_nm1sel + 'ST_log' + ext,
        p_nm1sel + 'n_btags_log' + ext,
    ),
    'Nm1Sel_2dcut.tex': (
        p_nm1sel + 'TwoDCut_px_log' + ext,
        p_nm1sel + 'TwoDCut_py_log' + ext,
    ),
}

plain_files = {
    'cutflow_tabular.tex':
        p_cutflow + 'CutflowTableTex/cutflow_tabular.tex',
    'cutflow_stack'+target_ext:
        p_cutflow + 'CutflowStack/Cutflow_cutflow'+ext,
}


dirname = 'AutoContentVLQPreSel'
tex_content = varial.tools.TexContent(
    images,
    plain_files,
    r"\includegraphics[width=0.49\textwidth]{" + dirname + "/%s}",
    dest_dir=None,
    name=dirname,
)


if __name__ == '__main__':
    varial.tools.Runner(tex_content)