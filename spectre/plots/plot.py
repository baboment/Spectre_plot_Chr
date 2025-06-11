import logging
from spectre.util import logger
import matplotlib.pyplot as plot_engine
from matplotlib import gridspec
import numpy as np
from typing import Optional


class CoveragePlot:
    def __init__(self, as_dev=False):
        logging.getLogger('matplotlib.font_manager').disabled = True
        self.logger = logger.setup_log(__name__, as_dev)
        # the plot
        self.figure = plot_engine.figure(figsize=(8, 4))
        gs = gridspec.GridSpec(1, 1)
        self.main_plot = plot_engine.subplot(gs[0])        # colors
        self.coverage_color = "#67a9cf"
        # legends
        self.main = ""
        self.x_axis = ""
        self.y_axis = ""
        self.file_prefix = "test"
        self.output_directory = "./"

    def plot_coverage(self, current_chromosome="", coverage=None):
        self.logger.debug("plotting coverage")
        # main plot
        self.main_plot.plot(np.array(coverage["pos"]), np.array(coverage["cov"]), color=self.coverage_color,
                            linewidth='0.5')
        # save and close
        self.figure.savefig(f'{self.output_directory}/{self.file_prefix}-{current_chromosome}.png', dpi=300)
        self.logger.info(f'Plot saved: {self.file_prefix}-{current_chromosome}.png')
        self.figure.clf()


class CNVPlot:
    def __init__(self, as_dev=False):
        self.logger = logger.setup_log(__name__, as_dev)
        # the plot
        self.figure = plot_engine.figure(figsize=(8, 4))
        gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
        self.main_plot = plot_engine.subplot(gs[0])
        self.candidates_plot = plot_engine.subplot(gs[1])
        self.candidates_plot.axes.get_yaxis().set_visible(False)
        # colors
        self.coverage_color = "#67a9cf"
        self.cnv_color = {"DUP": "#d73027", "DEL": "#1a9850"}
        # legends
        self.main = ""
        self.x_axis = ""
        self.y_axis = ""
        self.axis_ylim = {"bottom": 0, "top": 6}  # not showing over 6x coverage, min can not be lower than 0
        self.file_prefix = "test"
        self.output_directory = "./"

    def plot_coverage_cnv(self, current_chromosome="", stats=None, coverage=None, cnv_cand_list=None, bounds=None):
        if stats is None or coverage is None:
            self.logger.error("bot parameters are needed")
        self.logger.debug("plotting coverage + CNV")
        # main plot
        self.main_plot.plot(np.array(coverage["pos"]), np.array(coverage["cov"]), color=self.coverage_color,
                            linewidth='0.5')
        self.main_plot.axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])
        start = coverage["pos"][0]
        end = coverage["pos"][-1]
        self.candidates_plot.plot(np.array([start, end]), np.array([0, 0]), linewidth='0', color="#ffffff")
        # add CNV candidates
        self.logger.info(f"CNVs in chromosome: {current_chromosome}")
        if cnv_cand_list is not None:
            for cnv in cnv_cand_list:
                start = cnv.start
                end = cnv.end
                cnv_color = self.cnv_color[cnv.type]
                self.candidates_plot.plot(np.array([start, end]), np.array([0, 0]), linewidth='5', color=cnv_color)
        # save and close
        self.main_plot.plot(np.array([1, stats.chromosome_len]), np.array([stats.average, stats.average]),
                            linewidth='1', color="#000000")
        if bounds is not None:
            [upperb, lowerb] = bounds if len(bounds) == 2 else [np.NaN, np.NaN]
            self.main_plot.plot(np.array([1, stats.chromosome_len]), np.array([lowerb, lowerb]),
                                linewidth='1', color="#dd3497")
            self.main_plot.plot(np.array([1, stats.chromosome_len]), np.array([upperb, upperb]),
                                linewidth='1', color="#dd3497")
        self.figure.suptitle(f'{self.file_prefix} chromosome: {current_chromosome}')
        self.figure.savefig(f'{self.output_directory}/img/{self.file_prefix}_plot_cnv_{current_chromosome}.png', dpi=300)
        self.logger.info(f'Plot saved: img/{self.file_prefix}_plot_cnv_{current_chromosome}.png')
        self.figure.clf()


class GenomeCNVPlot:
    """Plot coverage and CNV calls for the whole genome in one figure."""

    def __init__(self, as_dev: bool = False):
        logging.getLogger('matplotlib.font_manager').disabled = True
        self.logger = logger.setup_log(__name__, as_dev)
        # wider figure for genome wide plots
        self.figure = plot_engine.figure(figsize=(32, 6))
        gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
        self.main_plot = plot_engine.subplot(gs[0])
        self.candidates_plot = plot_engine.subplot(gs[1])
        self.candidates_plot.axes.get_yaxis().set_visible(False)
        self.coverage_color = "#67a9cf"
        self.cnv_color = {"DUP": "#d73027", "DEL": "#1a9850"}
        self.axis_ylim = {"bottom": 0, "top": 6}
        self.file_prefix = "test"
        self.output_directory = "./"

    def plot_genome(self, coverage_per_chr: dict, cnv_per_chr: dict, chr_lengths: dict,
                    baseline: Optional[float] = None, bounds: Optional[list] = None):
        """Create genome wide CNV plot.

        Parameters
        ----------
        coverage_per_chr : dict
            Dictionary with chromosome -> {"pos": [...], "cov": [...]}.
        cnv_per_chr : dict
            Dictionary with chromosome -> list[CNVCandidate].
        chr_lengths : dict
            Dictionary with chromosome lengths.
        baseline : float, optional
            Genome wide baseline ploidy. If ``None`` the mean ploidy of all
            windows is used and shown as a green line.
        bounds : list, optional
            Lower and upper coverage bounds, used for plotting thresholds.
            Both thresholds are displayed in red.
        
        Notes
        -----
        A black horizontal line is drawn for each chromosome representing its
        average ploidy. The x-axis is scaled per chromosome starting at ``0``
        with tick marks every 20 Mbp and major labels every 100 Mbp.
        """

        if len(coverage_per_chr) == 0:
            self.logger.error("No coverage provided for genome plot")
            return

        chromosomes = list(coverage_per_chr.keys())

        all_pos = []
        all_cov = []
        window_cov = []
        xticks = []
        labels = []
        boundaries = []
        scale_ticks = []
        scale_labels = []

        # compute average ploidy for each chromosome for later plotting
        chr_means = {}

        offset = 0
        gap_size = 0  # no space between chromosomes
        for chrom in chromosomes:
            raw_pos = np.array(coverage_per_chr[chrom]["pos"])
            cov = np.array(coverage_per_chr[chrom]["cov"])
            step = np.median(np.diff(raw_pos)) if len(raw_pos) > 1 else 1
            win_green = max(1, int(round(100000 / step)))
            win_blue = max(1, int(round(20000 / step)))

            sm_cov = np.convolve(cov, np.ones(win_blue) / win_blue, mode="same")
            green_cov = np.convolve(cov, np.ones(win_green) / win_green, mode="same")

            pos = raw_pos + offset
            chr_means[chrom] = np.nanmean(cov)
            all_pos.append(pos)
            all_cov.append(sm_cov)
            window_cov.append(green_cov)
            length = chr_lengths.get(chrom, pos[-1] if len(pos) > 0 else 0)
            xticks.append(offset + length / 2)
            labels.append(chrom)
            start_off = offset
            offset += length
            boundaries.append(offset)
            tick = 0
            while tick <= length:
                scale_ticks.append(start_off + tick)
                if tick == 0:
                    scale_labels.append("0")
                elif tick % 100000000 == 0:
                    scale_labels.append(f"{int(tick/1000000)}m")
                else:
                    scale_labels.append("")
                tick += 20000000

        all_pos = np.concatenate(all_pos)
        all_cov = np.concatenate(all_cov)
        window_cov = np.concatenate(window_cov)

        # plot coverage and 100 kb window average
        self.main_plot.plot(all_pos, all_cov, color=self.coverage_color, linewidth='0.5')
        self.main_plot.scatter(all_pos, window_cov, color="#1a9850", s=3)
        self.main_plot.axes.set_ylim(bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"])

        # optional baseline across the genome and bounds
        genome_end = boundaries[-1]
        if baseline is not None:
            self.main_plot.plot(np.array([0, genome_end]),
                                np.array([baseline, baseline]),
                                linewidth='1', color="#0000ff")
        if bounds is not None and len(bounds) == 2:
            upperb, lowerb = bounds[1], bounds[0]
            self.main_plot.plot(np.array([0, genome_end]), np.array([lowerb, lowerb]),
                                linewidth='1', color="#d73027")
            self.main_plot.plot(np.array([0, genome_end]), np.array([upperb, upperb]),
                                linewidth='1', color="#d73027")

        self.main_plot.set_xlim(left=0, right=genome_end)
        self.candidates_plot.set_xlim(left=0, right=genome_end)
        self.main_plot.margins(x=0)
        self.candidates_plot.margins(x=0)

        # plot CNV segments and chromosome averages
        offset = 0
        for chrom in chromosomes:
            length = chr_lengths.get(chrom, 0)
            chr_mean = chr_means.get(chrom, np.nan)
            # chromosome baseline represented as a black line
            self.main_plot.plot(np.array([offset, offset + length]),
                                np.array([chr_mean, chr_mean]),
                                linewidth='1', color="#000000")
            if chrom in cnv_per_chr:
                for cnv in cnv_per_chr[chrom]:
                    start = cnv.start + offset
                    end = cnv.end + offset
                    cnv_color = self.cnv_color.get(cnv.type, "#000000")
                    self.candidates_plot.plot(np.array([start, end]), np.array([0, 0]),
                                              linewidth='5', color=cnv_color)
            offset += length

        # draw chromosome boundaries
        for boundary in boundaries[:-1]:
            self.main_plot.axvline(boundary, color="#303030", linewidth=1)
            self.candidates_plot.axvline(boundary, color="#303030", linewidth=1)

        self.main_plot.set_xticks(xticks)
        self.main_plot.set_xticklabels(labels, rotation=90, fontsize=6)
        self.candidates_plot.set_xticks(scale_ticks)
        self.candidates_plot.set_xticklabels(scale_labels, rotation=45, fontsize=6)

        self.figure.tight_layout()
        output_path = f'{self.output_directory}/img/{self.file_prefix}_plot_cnv_genome.png'
        self.figure.savefig(output_path, dpi=300)
        self.logger.info(f'Plot saved: {output_path}')
        self.figure.clf()

