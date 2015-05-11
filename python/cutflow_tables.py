
import varial.tools
import varial.util

import itertools
import re
import shutil
from math import sqrt


############################################ fetch histos and table content ###
cutflow_histos = varial.tools.HistoLoader(
    name='CutflowHistos',
    pattern='uhh2.*.root',
    filter_keyfunc=lambda w: 'cutflow' in w.in_file_path,
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
    input_result_path='CutflowHistos',
    save_log_scale=True,
)


class CutflowTableContent(varial.tools.Tool):
    """Generates cutflow table data."""
    can_reuse = False

    def __init__(self, name=None):
        super(CutflowTableContent, self).__init__()
        self.input_mc       = []
        self.input_data     = []
        self.head_line      = []
        self.table_data     = []
        self.table_mc_err   = []
        self.table_mc       = []
        self.titles_data    = []
        self.titles_mc      = []

    def get_input_histos(self):
        wrps = self.lookup_result('CutflowHistos')
        assert wrps
        mcee = itertools.ifilter(lambda w: not w.is_data, wrps)
        mcee = varial.gen.gen_norm_to_data_lumi(mcee)
        data = itertools.ifilter(lambda w: w.is_data, wrps)
        self.input_mc    = varial.gen.sort(mcee)
        self.input_data  = varial.gen.sort(data)

    def fill_head_line(self):
        bin_label = self.input_mc[0].histo.GetXaxis().GetBinLabel
        self.head_line = list(
            bin_label(i + 1)
            for i in xrange(self.input_mc[0].histo.GetNbinsX())
        )

    @staticmethod
    def _get_value_list(wrp):
        bin_cont = wrp.histo.GetBinContent
        return list(bin_cont(i + 1) for i in xrange(wrp.histo.GetNbinsX()))

    @staticmethod
    def _get_error_list(wrp):
        err_cont = wrp.histo.GetBinError
        return list(err_cont(i + 1) for i in xrange(wrp.histo.GetNbinsX()))

    def fill_tables(self):
        for wrp in self.input_data:
            self.titles_data.append(wrp.sample)
            self.table_data.append(self._get_value_list(wrp))
        for wrp in self.input_mc:
            self.titles_mc.append(wrp.sample)
            self.table_mc.append(self._get_value_list(wrp))
        for wrp in self.input_mc:
            self.table_mc_err.append(self._get_error_list(wrp))

    def _make_column_sum(self, table, squared = False):
        def gen_sum(A, B):
            return list(a + b for a,b in itertools.izip(A, B))
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
        self.table_data.append(self._make_column_sum(self.table_data))
        self.table_mc.append(self._make_column_sum(self.table_mc))
        self.table_mc_err.append(self._make_column_sum(self.table_mc_err, True))
        self.titles_data.append("Data Sum")
        self.titles_mc.append("MC Sum")

    def run(self):
        self.get_input_histos()
        self.fill_head_line()
        self.fill_tables()
        self.fill_sum_rows()
        self.result = self

    def mc_title_val_err_iterator(self):
        for title, vals, errs in itertools.izip(
            self.titles_mc,
            self.table_mc,
            self.table_mc_err,
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
        self.cont = self.lookup_result('CutflowTableContent')

    def make_header(self):
        line = self.sep.join(itertools.imap(lambda s: "%17s"%s,
                                            self.cont.head_line))
        line = 17*" " + self.sep + line + " \n"
        self.table_lines.append(line)

    def make_center(self):
        self.table_lines.append("\n")
        for title, vals, errs in self.cont.mc_title_val_err_iterator():
            zipped = ((a, b) for a, b in itertools.izip(vals, errs))
            self.table_lines.append(
                "%17s" % title
                + self.sep
                + self.sep.join("%8.1f +- %5.1f" % p for p in zipped)
                + " \n"
            )
        self.table_lines.append("\n")
        for title, vals in self.cont.data_title_value_iterator():
            self.table_lines.append(
                "%17s" % title
                + self.sep
                + self.sep.join("%17d" % v for v in vals)
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
tex_template = [
    r"\documentclass[10pt,fullpage]{article}",
    r"\pagestyle{empty}",
    r"\usepackage[landscape]{geometry}  % [margin=5mm, ]",
    r"\usepackage[usenames]{color} % Farbunterstuetzung",
    r"\usepackage{amssymb}	% Mathe",
    r"\usepackage{amsmath} % Mathe",
    r"\usepackage[utf8]{inputenc} % Direkte Eingabe von Umlauten",
    r"\begin{document}",
    r"\begin{table}",
    r"\input{cutflow_tabular.tex}",
    r"\end{table}",
    r"\end{document}",
]


class CutflowTableTex(varial.tools.Tool):
    """Reads cutflow histos from pool and creates latex table code."""

    def __init__(self, name=None):
        super(CutflowTableTex, self).__init__(name)
        self.cont           = None
        self.table_lines    = []
        self.target_dir     = ''
        self.sep            = ' $&$ '

    def configure(self):
        self.cont = self.lookup_result('CutflowTableContent')

    def make_header(self):
        self.table_lines += (
            r"\begin{tabular}{l | "
                + len(self.cont.head_line)*"r "
                + "}",
            17*" "
                + " & "
                + " & ".join(itertools.imap(
                    lambda s: "%17s" %
                              varial.analysis.get_pretty_name(s + "_tex"),
                    self.cont.head_line
                ))
                + r" \\",
            r"\hline",
            r"\hline",
        )

    def make_center(self):
        for title, vals, errs in self.cont.mc_title_val_err_iterator():
            self.table_lines += (
                "%17s" % title.replace("_", r"\_")
                    + " &$ "
                    + self.sep.join("%17.1f" % v for v in vals)
                    + r" $ \\",
                17*" "
                    + " &$ "
                    + self.sep.join("\\pm%17.1f" % e for e in errs)
                    + r" $ \\",
            )
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

        with open(self.cwd + "cutflow.tex", "w") as f:
            f.writelines(map(lambda l: l + "\n", tex_template))

#        subprocess.call(
#            ["pdflatex", "cutflow.tex"],
#            cwd=self.cwd
#        )

        wrp = wrappers.Wrapper(name="CutflowTableTex")
        for i, line in enumerate(self.table_lines):
            setattr(wrp, "line_%2d"%i, line)
        diskio.write(wrp, self.cwd + "cutflow_table.info")

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


cutflow_chain = varial.tools.ToolChain("CutflowTools", [
    cutflow_histos,
    cutflow_stack_plots,
    CutflowTableContent(),
    CutflowTableTxt(),
    CutflowTableTex(),
])

