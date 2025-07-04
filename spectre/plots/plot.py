import logging
from spectre.util import logger
import matplotlib.pyplot as plot_engine
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from typing import Optional


class CoveragePlot:
    def __init__(self, as_dev: bool = False):
        logging.getLogger('matplotlib.font_manager').disabled = True
        self.logger = logger.setup_log(__name__, as_dev)
        # the plot
        self.figure = plot_engine.figure(figsize=(8, 4))
        gs = gridspec.GridSpec(1, 1)
        self.main_plot = plot_engine.subplot(gs[0])
        # colors
        self.coverage_color = "#67a9cf"
        # legends
        self.main = ""
        self.x_axis = ""
        self.y_axis = ""
        self.file_prefix = "test"
        self.output_directory = "./"

    def plot_coverage(self, current_chromosome: str = "", coverage=None) -> None:
        self.logger.debug("plotting coverage")
        # main plot
        self.main_plot.plot(
            np.array(coverage["pos"]),
            np.array(coverage["cov"]),
            color=self.coverage_color,
            linewidth="0.5",
        )
        # save and close
        self.figure.savefig(
            f"{self.output_directory}/{self.file_prefix}-{current_chromosome}.png",
            dpi=300,
        )
        self.logger.info(
            f"Plot saved: {self.file_prefix}-{current_chromosome}.png"
        )
        self.figure.clf()


class CNVPlot:
    def __init__(self, as_dev: bool = False):
        self.logger = logger.setup_log(__name__, as_dev)
        # the plot
        self.figure = plot_engine.figure(figsize=(8, 4))
        gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
        self.main_plot = plot_engine.subplot(gs[0])
        self.candidates_plot = plot_engine.subplot(gs[1])
        self.candidates_plot.axes.get_yaxis().set_visible(False)
        # colors
        self.coverage_color = "#67a9cf"
        self.cnv_color = {"DUP": "#313695", "DEL": "#b2182b"}
        # legends
        self.main = ""
        self.x_axis = ""
        self.y_axis = ""
        self.axis_ylim = {"bottom": 0, "top": 6}
        self.file_prefix = "test"
        self.output_directory = "./"

    def plot_coverage_cnv(
        self,
        current_chromosome: str = "",
        stats=None,
        coverage=None,
        cnv_cand_list=None,
        bounds=None,
    ) -> None:
        if stats is None or coverage is None:
            self.logger.error("bot parameters are needed")
        self.logger.debug("plotting coverage + CNV")
        # main plot
        self.main_plot.plot(
            np.array(coverage["pos"]),
            np.array(coverage["cov"]),
            color=self.coverage_color,
            linewidth="0.5",
        )
        self.main_plot.axes.set_ylim(
            bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"]
        )
        start = coverage["pos"][0]
        end = coverage["pos"][-1]
        self.candidates_plot.plot(
            np.array([start, end]), np.array([0, 0]), linewidth="0", color="#ffffff"
        )
        # add CNV candidates
        self.logger.info(f"CNVs in chromosome: {current_chromosome}")
        if cnv_cand_list is not None:
            for cnv in cnv_cand_list:
                start = cnv.start
                end = cnv.end
                cnv_color = self.cnv_color[cnv.type]
                self.candidates_plot.plot(
                    np.array([start, end]),
                    np.array([0, 0]),
                    linewidth="5",
                    color=cnv_color,
                )
        # save and close
        self.main_plot.plot(
            np.array([1, stats.chromosome_len]),
            np.array([stats.average, stats.average]),
            linewidth="1",
            color="#000000",
        )
        if bounds is not None:
            [upperb, lowerb] = bounds if len(bounds) == 2 else [np.NaN, np.NaN]
            self.main_plot.plot(
                np.array([1, stats.chromosome_len]),
                np.array([lowerb, lowerb]),
                linewidth="1",
                color="#dd3497",
            )
            self.main_plot.plot(
                np.array([1, stats.chromosome_len]),
                np.array([upperb, upperb]),
                linewidth="1",
                color="#dd3497",
            )
        self.figure.suptitle(f"{self.file_prefix} chromosome: {current_chromosome}")
        self.figure.savefig(
            f"{self.output_directory}/img/{self.file_prefix}_plot_cnv_{current_chromosome}.png",
            dpi=300,
        )
        self.logger.info(
            f"Plot saved: img/{self.file_prefix}_plot_cnv_{current_chromosome}.png"
        )
        self.figure.clf()


class GenomeCNVPlot:
    """Plot coverage and CNV calls for the whole genome in one figure."""

    def __init__(self, as_dev: bool = False):
        logging.getLogger('matplotlib.font_manager').disabled = True
        self.logger = logger.setup_log(__name__, as_dev)
        # wider figure for genome wide plots
        self.figure = plot_engine.figure(figsize=(32, 6))
        # two panels, one for coverage and one for CNV candidates
        gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1])
        self.main_plot = plot_engine.subplot(gs[0])
        self.candidates_plot = plot_engine.subplot(gs[1])
        # hide y-axis on the CNV candidate plot
        self.candidates_plot.axes.get_yaxis().set_visible(False)
        self.coverage_color = "#67a9cf"
        self.cnv_color = {"DUP": "#313695", "DEL": "#b2182b"}
        self.ploidy_cmap = LinearSegmentedColormap.from_list(
            "ploidy",
            [
                "#b2182b",
                "#f46d43",
                "#fbea6cf9",
                "#2166ac",
                "#313695",
            ],
        )
        self.axis_ylim = {"bottom": 0, "top": 4}
        self.file_prefix = "test"
        self.output_directory = "./"

    def plot_genome(
        self,
        coverage_per_chr: dict,
        cnv_per_chr: dict,
        chr_lengths: dict,
        baseline: Optional[float] = None,
        bounds: Optional[list] = None,
    ) -> None:
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
        The coverage input is expected at approximately 1 kb intervals as
        produced by Mosdepth. For the genome wide plot these values are
        smoothed by averaging in a window of roughly one megabase derived from
        the median spacing of the coverage data.
        A black horizontal line is drawn for each chromosome representing its
        average ploidy. The x-axis is scaled per chromosome starting at ``0``
        with tick marks every 20 Mbp and major labels every 100 Mbp. Scatter
        points are coloured with a nine step palette ranging from ``#d73027`` at
        zero (deletions) through ``#ffffbf`` at two to ``#4575b4`` at four
        (duplications).
        """

        if len(coverage_per_chr) == 0:
            self.logger.error("No coverage provided for genome plot")
            return

        chromosomes = list(coverage_per_chr.keys())

        all_pos = []
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
            win_green = max(1, int(round(1000000 / step)))

            # smooth coverage using a simple moving average of approx 1 Mb
            green_cov = np.convolve(cov, np.ones(win_green) / win_green, mode="same")

            pos = raw_pos + offset
            chr_means[chrom] = np.nanmean(cov)
            all_pos.append(pos)
            window_cov.append(green_cov)
            length = chr_lengths.get(chrom, pos[-1] if len(pos) > 0 else 0)
            xticks.append(offset + length / 2)
            labels.append(chrom)
            start_off = offset
            if cnv_per_chr is not None and chrom in cnv_per_chr:
                for cnv in cnv_per_chr[chrom]:
                    cnv_start = start_off + cnv.start
                    cnv_end = start_off + cnv.end
                    self.candidates_plot.plot(
                        [cnv_start, cnv_end],
                        [0, 0],
                        linewidth=5,
                        color=self.cnv_color[cnv.type],
                    )
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
        window_cov = np.concatenate(window_cov)

        # plot genome coverage averaged in 1 Mb windows using colour coded scatter
        scatter = self.main_plot.scatter(
            all_pos,
            window_cov,
            c=window_cov,
            cmap=self.ploidy_cmap,
            vmin=0,
            vmax=4,
            s=0.1,
            zorder=2,
        )
        self.main_plot.axes.set_ylim(
            bottom=self.axis_ylim["bottom"], top=self.axis_ylim["top"]
        )

        # optional baseline across the genome and bounds
        genome_end = boundaries[-1]
        if baseline is not None:
            self.main_plot.plot(
                np.array([0, genome_end]),
                np.array([baseline, baseline]),
                linewidth="1",
                color="#3f51b5",
            )
        if bounds is not None and len(bounds) == 2:
            upperb, lowerb = bounds[1], bounds[0]
            self.main_plot.plot(
                np.array([0, genome_end]),
                np.array([lowerb, lowerb]),
                linewidth="1",
                color="#f44336",
            )
            self.main_plot.plot(
                np.array([0, genome_end]),
                np.array([upperb, upperb]),
                linewidth="1",
                color="#f44336",
            )

        self.main_plot.set_xlim(left=0, right=genome_end)
        self.main_plot.margins(x=0)

        # plot chromosome averages as black lines
        offset = 0
        for chrom in chromosomes:
            length = chr_lengths.get(chrom, 0)
            chr_mean = chr_means.get(chrom, np.nan)
            self.main_plot.plot(
                np.array([offset, offset + length]),
                np.array([chr_mean, chr_mean]),
                linewidth=1,
                color="#000000",
            )
            offset += length

        # draw chromosome boundaries
        for boundary in boundaries[:-1]:
            self.main_plot.axvline(boundary, color="#303030", linewidth=1)

        self.main_plot.set_xticks(scale_ticks)
        self.main_plot.set_xticklabels(scale_labels, rotation=45, fontsize=8)

        chr_axis = self.main_plot.secondary_xaxis("top")
        chr_axis.set_xticks(xticks)
        chr_axis.set_xticklabels(labels, rotation=90, fontsize=10)

        cbar = self.figure.colorbar(
            scatter,
            ax=self.main_plot,
            label="Ploidy",
            ticks=[0, 1, 2, 3, 4],
            fraction=0.025,
            pad=0.01,
        )

        # synchronize candidate subplot with coverage plot
        self.candidates_plot.set_xlim(self.main_plot.get_xlim())
        self.candidates_plot.set_xticks(scale_ticks)
        self.candidates_plot.set_xticklabels(scale_labels, rotation=45, fontsize=8)

        self.figure.suptitle(self.file_prefix)
        self.figure.tight_layout(rect=[0, 0, 0.99, 0.95])
        output_path = f"{self.output_directory}/img/{self.file_prefix}_plot_cnv_genome.png"
        self.figure.savefig(output_path, dpi=350)
        self.logger.info(f"Plot saved: {output_path}")
        self.figure.clf()
