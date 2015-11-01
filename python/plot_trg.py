#!/usr/bin/env python

import ROOT
import varial
import varial.tools
import varial.util

# varial.settings.max_num_processes = 1
varial.raise_root_error_level()
# varial.settings.rootfile_postfixes += ['.pdf', '.eps', '.C']
varial.settings.canvas_size_x = 1138
varial.settings.canvas_size_y = 744
#varial.settings.root_style.SetPadLeftMargin(0.1)
#varial.settings.root_style.SetPadRightMargin(0.1)
#varial.settings.root_style.SetPadTopMargin(0.1)
#varial.settings.root_style.SetPadBottomMargin(0.1)
y_axis_max_val = 1.4
varial.settings.defaults_Legend.update({
    'x_pos': 0.5,
    'y_pos': 0.8,
    'label_width': 0.2,
    'label_height': 0.04,
    'opt': 'f',
    'opt_data': 'p',
    'reverse': True
})

bins_lepton_pt = (
      range(0, 60, 4)
    + range(60, 100, 10)
    + range(100, 200, 20)
    + range(200, 350, 50)
    + [400]
)
bins_lead_jet_pt = (
      range(50, 100, 20)
    + range(100, 300, 10)
    + range(300, 400, 20)
    + range(400, 650, 50)
    + [700]
)
bins_sublead_jet_pt = (
      range(0, 250, 10)
    + range(250, 350, 20)
    + range(350, 450, 50)
    + [500]
)
bins_st = (
      [100, 200, 250, 300, 350]
    + range(400, 1000, 20)
    + range(1000, 1400, 50)
    + range(1400, 2200, 200)
)


def rebin_plots(wrps):
    for w in wrps:
        if '1leptonPt' in w.in_file_path:
            w = varial.op.rebin(w, bins_lepton_pt, True)
        elif '2leadJetPt' in w.in_file_path:
            w = varial.op.rebin(w, bins_lead_jet_pt, True)
        elif '3subleadJetPt' in w.in_file_path:
            w = varial.op.rebin(w, bins_sublead_jet_pt, True)
        elif '4ST' in w.in_file_path:
            w = varial.op.rebin(w, bins_st, True)
        yield w


def format_y_axis(wrps):
    for w in wrps:
        w.obj.GetYaxis().SetTitle('Efficiency')
        w.val_y_max = y_axis_max_val
        yield w


def set_legend_name(wrps):
    for w in wrps:
        w.legend = w.in_file_path.split('/')[0]

        if w.name.endswith('Eff'):
            if 'COMBO' in w.in_file_path:
                if w.in_file_path[:2] == 'El':
                    w.legend = 'Trigger combination for electrons'
                else:
                    w.legend = 'Trigger combination for muons'
            else:
                if w.in_file_path[:2] == 'El':
                    w.legend = 'Non-iso. ele. p_{T} > 45 GeV, cleaned AK4 jets: p_{T} [> 200, > 50] GeV'
                else:
                    w.legend = 'Non-iso. muon p_{T} > 40 GeV, cleaned AK4 jets: p_{T} [> 200, > 50] GeV'

        else:
            if w.in_file_path[:2] == 'El':
                w.legend = 'Raw distribution (only electrons, MC: T\'(M=800) b > t H b)'
                #if 'ST' in w.in_file_path and 'COMBO' in w.in_file_path:
                if 'PFJet' in w.in_file_path:
                    continue
            else:
                continue
                # w.draw_option_legend = ''  # do not make legend entry

        yield w

def set_colors(wrps):
    for w in wrps:
        if w.name.endswith('Eff'):
            if 'IsoMu' in w.in_file_path or ('_WP75_' in w.in_file_path):
                col = ROOT.kBlue + 3
            elif 'PFJet' in w.in_file_path:
                col = ROOT.kGreen + 1
            else:
                col = ROOT.kRed

            if '_OR_' in w.in_file_path:
                col -= 3

            if 'VVVL' in w.in_file_path:
                col = ROOT.kMagenta

            if any(w.in_file_path.startswith(t) for t in ('PFHT', 'trg_PFHT')):
                col = ROOT.kOrange

            w.obj.SetMarkerColor(col)
            w.obj.SetLineColor(col)

        else:
            # these are the raw distributions
            w.obj.SetFillColor(17)
            w.obj.SetLineColor(15)
            #w.obj.SetLineWidth(0)

        yield w

def set_plot_style(wrps):
    for w in wrps:
        if w.name.endswith('Eff') and 'TH1' not in w.type:
            # here are the TGraphAsymmErrors (only for error, not line)
            w.obj.SetMarkerStyle(1)
            w.obj.SetLineWidth(2)
            w.draw_option_legend = 'LP'
            w.draw_option = 'ZP'

        elif w.name.endswith('Eff'):
            # here are the TH1F (only for line, not errors)
            w.obj.SetMarkerStyle(1)
            w.obj.SetLineWidth(2)
            w.draw_option_legend = None
            w.draw_option = 'hist'

        elif w.name.endswith('_tot'):
            w.draw_option_legend = None

        yield w

def plot_maker(wrps):
    # wrps = rebin_plots(wrps)
    wrps = varial.gen.gen_add_wrp_info(
        wrps,
        name=lambda w: w.name + (
            '_sub' if ('/passing/' in w.in_file_path) else '_tot'),
        in_file_path=lambda w: w.in_file_path + (
            '_sub' if ('/passing/' in w.in_file_path) else '_tot'),
    )
    wrps = varial.gen.gen_add_wrp_info(
        wrps,
        in_file_path=lambda w: w.in_file_path.split('/')[0] + '/' + w.in_file_path.split('/')[2]
    )
    wrps = varial.gen.gen_make_eff_graphs(
        wrps, new_postfix='Eff', yield_everything=True)
    wrps = filter(
        lambda w: w.name.endswith('Eff') or w.name.endswith('_tot'), wrps)
    # wrps = varial.gen.gen_make_eff_graphs(
    #    wrps, new_postfix='Eff', eff_func=varial.op.div)
    wrps = varial.gen.imap_conditional(
        wrps,
        lambda w: 'Eff' not in w.name,
        varial.op.norm_to_integral,
    )
    wrps = format_y_axis(wrps)
    # wrps = set_legend_name(wrps)
    wrps = varial.gen.gen_add_wrp_info(wrps, legend=lambda w: w.in_file_path.split('/')[0])
    wrps = set_colors(wrps)
    wrps = set_plot_style(wrps)
    return wrps


def plot_grouper(wrps):
    #group_key = lambda w: str(w.in_file_path.split('/')[0][-3:]) + w.name.replace('Eff', '')
    group_key = lambda w: w.name.replace('Eff', '').replace('_tot', '')
    wrps = sorted(wrps, key=lambda w: w.name[::-1])  # All Eff stuff to back
    # for w in wrps: print w.name
    wrps = sorted(wrps, key=group_key)
    # for w in wrps: print w.name
    wrps = varial.gen.group(wrps, group_key)

    # efficiencies to be plotted on top of the raw distributions
    wrps = (sorted(w, key=lambda w: w.name.endswith('Eff')) for w in wrps)
    return wrps


def mk_txtbx_cms():
    txtbx = ROOT.TLatex(
        0.6137566,0.3793003,"CMS")  # default y: 0.6793003
    txtbx.SetNDC()
    txtbx.SetTextFont(61)
    txtbx.SetTextSize(0.06122449)
    txtbx.SetLineWidth(2)
    return txtbx


def mk_txtbx_sim():
    txtbx = ROOT.TLatex(
        0.617284,0.3355685,"Simulation Preliminary")  # default y: 0.6355685
    txtbx.SetNDC()
    txtbx.SetTextFont(52)
    txtbx.SetTextSize(0.0451895)
    txtbx.SetLineWidth(2)
    return txtbx


def mk_txtbx_2015():
    txtbx = ROOT.TLatex(0.9,0.92," 2015, 13 TeV")
    txtbx.SetNDC()
    txtbx.SetTextAlign(31)
    txtbx.SetTextFont(42)
    txtbx.SetTextSize(0.06)
    txtbx.SetLineWidth(2)
    return txtbx


class MyPlotStyler(varial.util.Decorator):
    def do_final_cosmetics(self):
        self.decoratee.do_final_cosmetics()
        h = self.first_drawn
        h.GetXaxis().SetLabelFont(42)
        h.GetXaxis().SetLabelSize(0.035)
        h.GetXaxis().SetTitleSize(0.035)
        h.GetXaxis().SetTitleOffset(1.4)
        h.GetXaxis().SetTitleFont(42)
        h.GetYaxis().SetLabelFont(42)
        h.GetYaxis().SetLabelSize(0.035)
        h.GetYaxis().SetTitleSize(0.035)
        h.GetYaxis().SetTitleOffset(1.0)
        h.GetYaxis().SetTitleFont(42)


def plotter_factory(**kws):
    kws['hook_loaded_histos'] = plot_maker
    kws['plot_setup'] = None
    kws['plot_grouper'] = plot_grouper
    kws['save_name_func'] = lambda c: c.name
    #kws['canvas_decorators'] = [
    #    varial.rnd.Legend(
    #        # xy_coords=(0.3294533,0.1676385,0.574162,0.2959184),
    #        x_pos=(0.574162 - 0.3294533)/2. + 0.3294533,
    #        label_width=(0.574162 - 0.3294533),
    #        y_pos=(0.2959184 - 0.1676385)/2. + 0.1676385,
    #        label_height=0.035,
    #        text_size=0.028,
    #        text_font=42,
    #    ),
    #    varial.rnd.TextBox(textbox=mk_txtbx_cms()),
    #    varial.rnd.TextBox(textbox=mk_txtbx_sim()),
    #    varial.rnd.TextBox(textbox=mk_txtbx_2015()),
    #    MyPlotStyler(),
    #]
    return varial.tools.Plotter(**kws)


def input_filter_keyfunc_general(w):
    return (
        w.in_file_path.count('/') == 2  # only trg_name/type/variable allowed
        and 'leading' not in w.name  # don't look at jet right now
    )


def input_filter_keyfunc_ele(w):
    ifp = w.in_file_path
    res = input_filter_keyfunc_general(w) and (
        any(w.in_file_path.startswith(t) for t in ('PFHT', 'trg_PFHT')) or
        ifp.startswith('trg_Ele') or ifp.startswith('Ele')
    )
    return res
 

def input_filter_keyfunc_mu(w):
    ifp = w.in_file_path
    return input_filter_keyfunc_general(w) and (
        any(w.in_file_path.startswith(t) for t in ('PFHT', 'trg_PFHT')) or
        not (ifp.startswith('trg_Ele') or ifp.startswith('Ele')) 
    )


def mk_plttrs(path='*.root'):
    plot_ele = varial.tools.mk_rootfile_plotter(
        pattern=path,
        filter_keyfunc=input_filter_keyfunc_mu,
        plotter_factory=plotter_factory,
        flat=True,
        name='VLQTrigMu',
        auto_legend=False,
    )
    plot_mu = varial.tools.mk_rootfile_plotter(
        pattern=path,
        filter_keyfunc=input_filter_keyfunc_ele,
        plotter_factory=plotter_factory,
        flat=True,
        name='VLQTrigEle',
        auto_legend=False,
    )
    return [plot_ele, plot_mu]

if __name__ == '__main__':
    varial.tools.Runner(varial.tools.ToolChainParallel('VLQTrig', mk_plttrs()))
    varial.tools.WebCreator().run()


# HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50 : Non-isolated electron p_{T} > 45 GeV, cleaned AK4 jet p_{T} > 200 GeV, cleaned AK4 jet p_{T} > 50 GeV
# HLT_Mu40_eta2p1_PFJet200_PFJet50 : Non-isolated muon p_{T} > 40 GeV, cleaned AK4 jet p_{T} > 200 GeV, cleaned AK4 jet p_{T} > 50 GeV
# HLT_PFHT800 : H_{T} > 800 GeV
