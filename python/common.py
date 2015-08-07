import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = kError;')


import varial
import varial.history


def label_axes(wrps):
    for w in wrps:
        if 'TH1' in w.type and w.histo.GetXaxis().GetTitle() == '':
            w.histo.GetXaxis().SetTitle(w.histo.GetTitle())
            w.histo.GetYaxis().SetTitle('events')
            w.histo.SetTitle('')
        yield w


signal_indicators = ['_M800_', '_M1000_', '_M1200_']


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
                print 'WARNING In merge_decay_channels: buffer not empty. ' \
                      'Flushing remaining items:' + str(buf)
                yield do_merging(buf)
            yield w


def is_signal(name):
    return any(t in name for t in ['_M800_', '_M1000_', '_M1200_'])


def yield_n_objs(wrps, n=50):
    i = 0
    for w in wrps:
        i += 1
        if i < n:
            yield w
        else:
            break
