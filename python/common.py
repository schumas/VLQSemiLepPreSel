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


def get_samplename(wrp):
    if hasattr(wrp, 'sample') and wrp.sample:
        return wrp.sample
    fname = os.path.basename(wrp.file_path)
    if fname.startswith('uhh2'):
        return fname.split('.')[-2]
    else:
        return os.path.splitext(fname)[0]


def get_legend(wrp, sig_ind):
    smpl = get_samplename(wrp)
    if 'Run20' in smpl:
        return 'Data'
    elif any(s in smpl for s in sig_ind):
        mass = smpl.split('_')[-1]
        if mass.startswith('M'):
            mass = mass[1:]
        hnd = 'rh' if '_RH_' in smpl else 'lh'
        return 'T_{%s}(%d)#rightarrowtH' % (hnd, int(mass))
    else:
        return smpl


def get_sys_info(w):
    def get_info(tok):
        if tok in w.file_path:
            return next(s
                        for s in w.file_path.split('/') 
                        if s.endswith(tok))
        else:
            return ''
    return get_info('__minus') or get_info('__plus')


def add_wrp_info(wrps, sig_ind=None):
    sig_ind = sig_ind or signal_indicators
    return varial.generators.gen_add_wrp_info(
        wrps,
        sample=get_samplename,
        legend=lambda w: get_legend(w, sig_ind),
        is_signal=lambda w: any(s in w.sample for s in sig_ind),
        is_data=lambda w: 'Run20' in w.sample,
        variable=lambda w: w.in_file_path.split('/')[-1],
        sys_info=get_sys_info,
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
        self.result.__dict__.update({
            'brs' : self.brs,
            # 'masses' : list(int(x) for x in self.result.res_exp_x)
            }
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
