import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = kError;')


import varial
import varial.history
import varial.tools
import glob
import os


def label_axes(wrps):
    for w in wrps:
        if 'TH1' in w.type and w.histo.GetXaxis().GetTitle() == '':
            w.histo.GetXaxis().SetTitle(w.histo.GetTitle())
            w.histo.GetYaxis().SetTitle('events')
            w.histo.SetTitle('')
        yield w


signal_indicators = ['_TH_', 'TpTp_',]


def get_samplename(fname):
    fname = os.path.basename(fname)
    if fname.startswith('uhh2'):
        return fname.split('.')[-2]
    else:
        return os.path.splitext(fname)[0]


def add_wrp_info(wrps, sig_ind=signal_indicators):
    return varial.generators.gen_add_wrp_info(
        wrps,
        sample=lambda w: get_samplename(w.file_path),
        legend=lambda w: w.sample,
        is_signal=lambda w: any(s in w.sample for s in sig_ind),
        is_data=lambda w: 'Run20' in w.sample,
        variable=lambda w: w.in_file_path.split('/')[-1]
    )


def merge_decay_channels(wrps, postfixes=('_Tlep', '_NonTlep'), suffix='', print_warning=True):
    """histos must be sorted!!"""

    @varial.history.track_history
    def merge_decay_channel(w):
        return w

    def do_merging(buf):
        res = varial.operations.merge(buf)
        res.sample = next(res.sample[:-len(p)]
                          for p in postfixes
                          if res.sample.endswith(p))+suffix
        res.legend = next(res.legend[:-len(p)]
                          for p in postfixes
                          if res.legend.endswith(p))+suffix
        res.file_path = ''
        del buf[:]
        return merge_decay_channel(res)

    buf = []
    for w in wrps:
        if any(w.sample.endswith(p) for p in postfixes):
            buf.append(w)
            if len(buf) == len(postfixes):
                yield do_merging(buf)
        else:
            if buf:
                if print_warning:
                    print 'WARNING In merge_decay_channels: buffer not empty.\n' \
                          'postfixes:\n' + str(postfixes) + '\n' \
                          'Flushing remaining items:\n' + '\n'.join(
                        '%s, %s' % (w.sample, w.in_file_path) for w in buf
                    )
                yield do_merging(buf)
            yield w
    if buf:
        yield do_merging(buf)

def is_signal(name):
    return any(t in name for t in signal_indicators)


def yield_n_objs(wrps, n=50):
    i = 0
    for w in wrps:
        i += 1
        if i < n:
            yield w
        else:
            break


######################################################### limit calculation ###
try:
    from varial.extensions.limits import ThetaLimits
except ImportError:
    from varial import monitor
    monitor.message(
        'UHH2.VLQSemiLepPreSel.common',
        'WARNING Theta is not working'
    )
    ThetaLimits = object


class TpTpThetaLimits(ThetaLimits):
    def __init__(self, brs=None, *args ,**kws):
        super(TpTpThetaLimits, self).__init__(*args, **kws)
        self.brs = brs

    def run(self):
        super(TpTpThetaLimits, self).run()
        self.result = varial.wrappers.Wrapper(
            name=self.result.name,
            res_exp_x=self.result.res_exp_x,  # TODO only TObjects or native python objects (list, dict, int, str ...) allowed
            res_exp_y=self.result.res_exp_y,
            res_exp_xerrors=self.result.res_exp_xerrors,
            res_exp_yerrors=self.result.res_exp_yerrors,
            brs=self.brs,
            masses=list(int(x) for x in self.result.res_exp_x)
        )


class TriangleLimitPlots(varial.tools.Tool):
    def __init__(self,
        name=None,
        limit_rel_path=''
    ):
        super(TriangleLimitPlots, self).__init__(name)
        self.limit_rel_path = limit_rel_path

    def make_tri_hist(self, wrps, mass_ind):
        tri_hist = ROOT.TH2F("triangular_limits", ";br to th;br to tz", 11, -0.05, 1.05, 11, -0.05, 1.05)
        for w in wrps:
            br_th = float(w.brs['th'])
            br_tz = float(w.brs['tz'])
            tri_hist.Fill(br_th, br_tz, w.res_exp_y[mass_ind])
        return varial.wrappers.HistoWrapper(tri_hist,
                legend='M-'+str(w.masses[mass_ind]),
                mass=w.masses[mass_ind])


    def run(self):
        # parents = os.listdir(self.cwd+'/..')
        theta_tools = glob.glob(os.path.join(self.cwd+'..', self.limit_rel_path))
        wrps = list(self.lookup_result(k) for k in theta_tools)
        filename = os.path.join(varial.analysis.cwd, self.name + ".root")
        # f = ROOT.TFile.Open(filename, "RECREATE")
        # f.cd()
        list_hists=[]
        for i, m in enumerate(wrps[0].masses):
            list_hists.append(self.make_tri_hist(wrps, i))
        # tri_hist.Write()
        self.result = list_hists
        # f.Close()
