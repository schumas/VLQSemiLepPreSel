import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = kError;')


import varial
import varial.history
import varial.tools
from varial.extensions.limits import ThetaLimits


def label_axes(wrps):
    for w in wrps:
        if 'TH1' in w.type and w.histo.GetXaxis().GetTitle() == '':
            w.histo.GetXaxis().SetTitle(w.histo.GetTitle())
            w.histo.GetYaxis().SetTitle('events')
            w.histo.SetTitle('')
        yield w


signal_indicators = ['_TH_', 'TpTp_',]


def add_wrp_info(wrps):
    return varial.generators.gen_add_wrp_info(
        wrps,
        sample=lambda w: w.file_path.split('.')[-2],
        legend=lambda w: w.sample,
        is_signal=lambda w: any(s in w.sample for s in signal_indicators),
        is_data=lambda w: 'Run20' in w.sample,
    )


def merge_decay_channels(wrps, postfixes=('_Tlep', '_NonTlep')):
    """histos must be sorted!!"""

    @varial.history.track_history
    def merge_decay_channel(w):
        return w

    def do_merging(buf):
        res = varial.operations.merge(buf)
        res.sample = next(res.sample[:-len(p)]
                          for p in postfixes
                          if res.sample.endswith(p))
        res.legend = next(res.legend[:-len(p)]
                          for p in postfixes
                          if res.legend.endswith(p))
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
                print 'WARNING In merge_decay_channels: buffer not empty.\n' \
                      'postfixes:\n' + str(postfixes) + '\n' \
                      'Flushing remaining items:\n' + '\n'.join(
                    '%s, %s' % (w.sample, w.in_file_path) for w in buf
                )
                yield do_merging(buf)
            yield w


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





class TpTpThetaLimits(ThetaLimits):
    def __init__(self, brs=None, *args ,**kws):
        super(TpTpThetaLimits, self).__init__(*args, **kws)
        self.brs = brs

    def run(self):
        super(TpTpThetaLimits, self).run()
        self.result = varial.wrappers.Wrapper(
            name=self.result.name,
            _res_exp=self.result._res_exp,
            _res_obs=self.result._res_obs,
            brs=self.brs
        )


class TriangleLimitPlots(varial.tools.Tool):
    def __init__(self, name=None):
        super(TriangleLimitPlots, self).__init__(name)


    def run(self):
        # parent = varial.analysis.lookup_tool('../.')
        # varial.analysis.print_tool_tree()
        parents = os.listdir(self.cwd+'/..')
        # print parents
        theta_tools = list(k for k in parents if k.startswith("ThetaLimit"))
        # print theta_tools
        wrps = list(self.lookup_result('../' + k) for k in theta_tools)
        filename = os.path.join(varial.analysis.cwd, self.name + ".root")
        f = ROOT.TFile.Open(filename, "RECREATE")
        f.cd()
        tri_hist = ROOT.TH2F("triangular_limits", ";br to th;br to tz", 10, 0., 1., 10, 0., 1.)
        for w in wrps:
            br_th = float(w.brs['th'])
            br_tz = float(w.brs['tz'])
            # limit_f = float(w.res_exp.y[0])
            tri_hist.Fill(br_th, br_tz, w.res_exp.y[0])
        tri_hist.Write()
        f.Close()


