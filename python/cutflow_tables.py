from math import sqrt
import subprocess
import itertools
import shutil

import varial.tools
import varial.util
import varial.plotter
import varial.generators as gen


############################################ fetch histos and table content ###
class CutflowTableContent(varial.tools.Tool):
    """Generates cutflow table data."""
    can_reuse = False

    def __init__(self, name=None):
        super(CutflowTableContent, self).__init__(name)
        self.head_line      = []

        self._input_data    = []
        self.titles_data    = []
        self.table_data     = []

        self._input_mc      = []
        self.titles_mc      = []
        self.table_mc       = []
        self.table_mc_err   = []

        self._input_sig     = []
        self.titles_sig     = []
        self.table_sig      = []
        self.table_sig_err  = []

    def get_input_histos(self):
        wrps = self.lookup_result('../CutflowHistos')
        assert wrps

        mcee = itertools.ifilter(lambda w: w.is_background, wrps)
        mcee = varial.gen.gen_norm_to_data_lumi(mcee)
        self._input_mc = varial.gen.sort(mcee)

        sig = itertools.ifilter(lambda w: w.is_signal, wrps)
        sig = varial.gen.gen_norm_to_data_lumi(sig)
        self._input_sig = varial.gen.sort(sig)

        data = itertools.ifilter(lambda w: w.is_data, wrps)
        self._input_data = varial.gen.sort(data)

    def fill_head_line(self):
        bin_label = self._input_mc[0].histo.GetXaxis().GetBinLabel
        self.head_line = list(
            bin_label(i + 1)
            for i in xrange(self._input_mc[0].histo.GetNbinsX())
            if bin_label(i + 1)
        )

    @staticmethod
    def _get_value_list(wrp):
        bin_cont = wrp.histo.GetBinContent
        bin_label = wrp.histo.GetXaxis().GetBinLabel
        return list(
            bin_cont(i + 1)
            for i in xrange(wrp.histo.GetNbinsX())
            if bin_label(i + 1)
        )

    @staticmethod
    def _get_error_list(wrp):
        err_cont = wrp.histo.GetBinError
        bin_label = wrp.histo.GetXaxis().GetBinLabel
        return list(
            err_cont(i + 1)
            for i in xrange(wrp.histo.GetNbinsX())
            if bin_label(i + 1)
        )

    def fill_tables(self):
        for wrp in self._input_data:
            self.titles_data.append(wrp.sample)
            self.table_data.append(self._get_value_list(wrp))

        for wrp in self._input_mc:
            self.titles_mc.append(wrp.sample)
            self.table_mc.append(self._get_value_list(wrp))
        for wrp in self._input_mc:
            self.table_mc_err.append(self._get_error_list(wrp))

        for wrp in self._input_sig:
            self.titles_sig.append(wrp.sample)
            self.table_sig.append(self._get_value_list(wrp))
        for wrp in self._input_sig:
            self.table_sig_err.append(self._get_error_list(wrp))

    @staticmethod
    def _make_column_sum(table, squared = False):
        def gen_sum(A, B):
            return list(a + b for a, b in itertools.izip(A, B))
        def gen_sum_sq(A, B):
            return map(sqrt, ((a*a) + (b*b) for a, b in itertools.izip(A, B)))
        row = []
        for ro in table:
            if not row:
                row = ro[:]
            elif squared:
                row = gen_sum_sq(row, ro)
            else:
                row = gen_sum(row, ro)
        return row

    def fill_sum_rows(self):
        if len(self._input_data) > 1:
            self.titles_data.append("Data Sum")
            self.table_data.append(self._make_column_sum(self.table_data))

        self.titles_mc.append("MC Sum")
        self.table_mc.append(self._make_column_sum(self.table_mc))
        self.table_mc_err.append(self._make_column_sum(self.table_mc_err, True))

    def calc_permillage(self):
        self.head_line.append("1000 X output/input")
        for row in self.table_data:
            if not row:
                continue
            row.append(1000. * row[-1] / row[0])

        for row, row_err in itertools.izip(self.table_mc, self.table_mc_err):
            if not (row and row_err):
                continue
            row.append(1000. * row[-1] / row[0])
            row_err.append(0.)  # TODO uncert for binomial

        for row, row_err in itertools.izip(self.table_sig, self.table_sig_err):
            if not (row and row_err):
                continue
            row.append(1000. * row[-1] / row[0])
            row_err.append(0.)  # TODO uncert for binomial

    def run(self):
        self.get_input_histos()
        self.fill_head_line()
        self.fill_tables()
        self.fill_sum_rows()
        self.calc_permillage()
        self.result = self

    def mc_title_val_err_iterator(self):
        for title, vals, errs in itertools.izip(
            self.titles_mc,
            self.table_mc,
            self.table_mc_err,
        ):
            yield title, vals, errs

    def sig_title_val_err_iterator(self):
        for title, vals, errs in itertools.izip(
            self.titles_sig,
            self.table_sig,
            self.table_sig_err,
        ):
            yield title, vals, errs

    def data_title_value_iterator(self):
        for title, vals in itertools.izip(
            self.titles_data,
            self.table_data
        ):
            yield title, vals


################################################################# txt table ###
class CutflowTableTxt(varial.tools.Tool):
    """Reads cutflow table and creates a txt table."""

    def __init__(self, name=None):
        super(CutflowTableTxt, self).__init__(name)
        self.cont           = None
        self.table_lines    = []
        self.sep            = ", "

    def configure(self):
        self.cont = self.lookup_result('../CutflowTableContent')

    def make_header(self):
        line = self.sep.join(itertools.imap(lambda s: "%24s"%s,
                                            self.cont.head_line))
        line = 32*" " + self.sep + line + " \n"
        self.table_lines.append(line)

    def _add_mc_stuff(self, iterator):
        for title, vals, errs in iterator:
            zipped = ((a, b) for a, b in itertools.izip(vals, errs))
            self.table_lines.append(
                "%32s" % title
                + self.sep
                + self.sep.join("%13.2f +- %7.2f" % p for p in zipped)
                + " \n"
            )

    def make_center(self):
        self.table_lines.append("\n")
        self._add_mc_stuff(self.cont.sig_title_val_err_iterator())
        self.table_lines.append("\n")
        self._add_mc_stuff(self.cont.mc_title_val_err_iterator())
        self.table_lines.append("\n")
        for title, vals in self.cont.data_title_value_iterator():
            self.table_lines.append(
                "%32s" % title
                + self.sep
                + self.sep.join("%24d" % v for v in vals)
            )
        self.table_lines.append("\n")

    def write_out(self):
        with open(self.cwd + "cutflow_table.txt", "w") as f:
            f.writelines(self.table_lines)
        wrp = varial.wrappers.Wrapper(name="CutflowTableTxt")
        for i, line in enumerate(self.table_lines):
            setattr(wrp, "line_%2d"%i, line)
        self.result = wrp

    def run(self):
        self.configure()
        self.make_header()
        self.make_center()
        self.write_out()


################################################################# tex table ###
tex_template = r"""
\documentclass[10pt,fullpage]{article}
\pagestyle{empty}
\usepackage[landscape]{geometry}  % [margin=5mm, ]
\usepackage[usenames]{color} % Farbunterstuetzung
\usepackage{amssymb}	% Mathe
\usepackage{amsmath} % Mathe
\usepackage[utf8]{inputenc} % Direkte Eingabe von Umlauten
\begin{document}
\begin{table}
\input{cutflow_tabular.tex}
\end{table}
\end{document}
"""


class CutflowTableTex(varial.tools.Tool):
    """Reads cutflow histos from pool and creates latex table code."""

    def __init__(self, compile_latex=False, name=None):
        super(CutflowTableTex, self).__init__(name)
        self.cont           = None
        self.table_lines    = []
        self.target_dir     = ''
        self.sep            = ' $&$ '
        self.compile_latex  = compile_latex

    def configure(self):
        self.cont = self.lookup_result('../CutflowTableContent')

    def make_header(self):
        self.table_lines += (
            r"\begin{tabular}{l | "
                + len(self.cont.head_line)*"r "
                + "}",
            30*" "
                + " & "
                + " & ".join(itertools.imap(
                    lambda s: "%30s" % varial.analysis.get_pretty_name(
                        s + "_tex").replace("_", r"\_"),
                    self.cont.head_line
                ))
                + r" \\",
            r"\hline",
            r"\hline",
        )

    def _add_mc_stuff(self, iterator):
        for title, vals, errs in iterator:
            self.table_lines += (
                "%30s" % title.replace("_", r"\_")
                    + " &$ "
                    + self.sep.join("%17.1f" % v for v in vals)
                    + r" $ \\",
                30*" "
                    + " &$ "
                    + self.sep.join("\\pm%17.1f" % e for e in errs)
                    + r" $ \\",
            )

    def make_center(self):
        self._add_mc_stuff(self.cont.sig_title_val_err_iterator())
        self.table_lines.append(r"\hline")
        self._add_mc_stuff(self.cont.mc_title_val_err_iterator())
        self.table_lines.insert(-2, r"\hline")
        self.table_lines.append(r"\hline")
        for title, vals in self.cont.data_title_value_iterator():
            self.table_lines.append(
                "%17s" % title.replace("_", r"\_")
                + " &$ "
                + " $&$ ".join("%17d" % v for v in vals)
                + r" $ \\"
            )
        self.table_lines.insert(-1, r"\hline")
        self.table_lines.append(r"\hline")

    def make_footer(self):
        self.table_lines += (
            r"\end{tabular}",
        )

    def write_out(self):
        with open(self.cwd + "cutflow_tabular.tex", "w") as f:
            f.writelines(map(lambda l: l + "\n", self.table_lines))

        if self.compile_latex:
            with open(self.cwd + "cutflow.tex", "w") as f:
                f.write(tex_template)

            ret = subprocess.call(
                ["pdflatex", "cutflow.tex"],
                cwd=self.cwd
            )
            if ret == 0:
                # convert to png for webcreator
                subprocess.call(
                    ["convert", "cutflow.pdf", "cutflow.png"],
                    cwd=self.cwd
                )
                # make cutflow.info, so the webcreator picks it up
                self.io.write(varial.wrappers.Wrapper(name='cutflow'))

        wrp = varial.wrp.Wrapper(name="CutflowTableTex")
        for i, line in enumerate(self.table_lines):
            setattr(wrp, "line_%2d"%i, line)
        self.result = wrp

    def copy_to_target_dir(self):
        if not self.target_dir:
            return
        self.message("INFO Copying cutflow_tabular.tex to " + self.target_dir)
        shutil.copy2(
            self.cwd + "cutflow_tabular.tex",
            self.target_dir
        )

    def run(self):
        self.configure()
        self.make_header()
        self.make_center()
        self.make_footer()
        self.write_out()
        self.copy_to_target_dir()


def rebin_cutflow(wrps):
    wrps = list(wrps)
    last_bin = 2
    for w in wrps:
        maybe_last_bin = 0
        for i in xrange(w.histo.GetNbinsX(), last_bin, -1):
            if w.histo.GetBinContent(i) == w.histo.GetBinContent(i-1):
                continue
            maybe_last_bin = i
            break

        if last_bin < maybe_last_bin:
            last_bin = maybe_last_bin

    return gen.gen_trim(wrps, left=False, right=last_bin)


def mk_cutflow_chain(input_pat, loader_hook, filter_keyfunc=None):

    cutflow_fk = lambda w: (w.in_file_path.split('/')[-1] == 'cutflow' and
                            not w.sample.startswith('Signal_'))

    if filter_keyfunc:
        final_fk = lambda w: filter_keyfunc(w) and cutflow_fk(w)
    else:
        final_fk = cutflow_fk

    cutflow_histos = varial.tools.HistoLoader(
        name='CutflowHistos',
        pattern=input_pat,
        filter_keyfunc=final_fk,
        hook_loaded_histos=lambda w: rebin_cutflow(loader_hook(w, 0)),
    )

    #class AxisTitles(varial.util.Decorator):
    #    def do_final_cosmetics(self):
    #        self.decoratee.do_final_cosmetics()
    #        self.first_drawn.GetXaxis().SetTitle("cutflow")
    #        if hasattr(self, "bottom_hist"):
    #            self.bottom_hist.GetXaxis().SetTitle("cutflow")
    #        self.first_drawn.GetYaxis().SetTitle("selected events / step")

    cutflow_stack_plots = varial.tools.Plotter(
        'CutflowStack',
        stack=True,
        input_result_path='../CutflowHistos',
    )

    cutflow_normed_plots = varial.tools.Plotter(
        'CutflowNormed',
        stack=False,
        plot_grouper=varial.plotter.plot_grouper_by_in_file_path,
        hook_loaded_histos=gen.gen_norm_to_max_val,
        input_result_path='../CutflowHistos',
        canvas_decorators=[varial.rendering.Legend]
    )

    return varial.tools.ToolChain("CutflowTools", [
        cutflow_histos,
        cutflow_stack_plots,
        # cutflow_normed_plots,
        CutflowTableContent(),
        CutflowTableTxt(),
        CutflowTableTex(),
    ])

