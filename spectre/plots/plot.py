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


import logging
from spectre.util import logger
import matplotlib.pyplot as plot_engine
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from typing import Optional
import os

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


import os

class GenomeCNVPlot:
    """Plot coverage and CNV calls for the whole genome in one figure."""

    def __init__(self, as_dev: bool = False):
        logging.getLogger('matplotlib.font_manager').disabled = True
        self.logger = logger.setup_log(__name__, as_dev)

        # Solid white background; constrained layout
        plot_engine.rcParams['savefig.facecolor'] = 'white'
        plot_engine.rcParams['figure.facecolor']  = 'white'
        self.figure = plot_engine.figure(figsize=(40, 5), constrained_layout=True, facecolor="white")

        gs = gridspec.GridSpec(1, 1, figure=self.figure)
        self.main_plot = plot_engine.subplot(gs[0])

        # Colors (your palette)
        self.coverage_color = "#67a9cf"
        self.cnv_color = {"DUP": "#313695", "DEL": "#b2182b"}
        self.ploidy_cmap = LinearSegmentedColormap.from_list(
            "Aneuploidy",
            ["#980f20", "#bb1f1f", "#e4b832", "#3e84de", "#064883"]
        )

        self.axis_ylim = {"bottom": 0, "top": 4}
        self.file_prefix = "test"
        self.output_directory = "./"

    def _smooth_same_len_safe(self, cov: np.ndarray, step: float) -> np.ndarray:
        """~500 kb moving average, robust, same length as input."""
        n = cov.size
        if n <= 3:
            return cov
        win = max(1, int(round(1_000_000 / max(step, 1e-9))))
        win = min(win, n)
        if win % 2 == 0:
            win = win - 1 if win > 1 else 1
        if win <= 1:
            return cov
        kernel = np.ones(win, dtype=float) / float(win)
        return np.convolve(cov, kernel, mode="same")

    def plot_genome(
        self,
        coverage_per_chr: dict,
        cnv_per_chr: dict,
        chr_lengths: dict,
        baseline: Optional[float] = None,
        bounds: Optional[list] = None,
    ) -> None:

        if not coverage_per_chr:
            self.logger.error("No coverage provided for genome plot")
            return

        # Ensure output folder exists
        out_dir = os.path.join(self.output_directory, "img")
        os.makedirs(out_dir, exist_ok=True)

        chromosomes = list(coverage_per_chr.keys())

        all_pos = []
        window_cov = []
        xticks, labels = [], []
        boundaries = []
        chr_means = {}

        offset = 0
        for chrom in chromosomes:
            raw_pos = np.asarray(coverage_per_chr[chrom].get("pos", []))
            cov     = np.asarray(coverage_per_chr[chrom].get("cov", []))

            if raw_pos.size == 0 or cov.size == 0:
                self.logger.warning(f"[{chrom}] empty coverage; skipping")
                length = int(chr_lengths.get(chrom, 0))
                xticks.append(offset + length / 2)
                # label even if empty, so spacing stays consistent
                labels.append(chrom)
                offset += length
                boundaries.append(offset)
                continue

            step = np.median(np.diff(raw_pos)) if raw_pos.size > 1 else 1.0
            green_cov = self._smooth_same_len_safe(cov, step)

            pos = raw_pos + offset
            chr_means[chrom] = float(np.nanmean(cov))
            all_pos.append(pos)
            window_cov.append(green_cov)

            length = int(chr_lengths.get(chrom, int(pos[-1]) if pos.size else 0))
            xticks.append(offset + length / 2)
            labels.append(chrom)

            offset += length
            boundaries.append(offset)

        if not all_pos:
            self.logger.error("All chromosomes were empty; nothing to plot.")
            return

        all_pos    = np.concatenate(all_pos)
        window_cov = np.concatenate(window_cov)

        # Scatter (rasterized so points don't vanish at high DPI)
        scatter = self.main_plot.scatter(
            all_pos, window_cov,
            c=window_cov, cmap=self.ploidy_cmap, vmin=0, vmax=4,
            s=4.0, marker='o', edgecolors='none', linewidths=0,
            rasterized=True, zorder=2
        )

        self.main_plot.set_ylim(self.axis_ylim["bottom"], self.axis_ylim["top"])

        self.main_plot.spines['bottom'].set_visible(False)
        self.main_plot.spines['top'].set_visible(False)
        
        self.main_plot.set_yticks([0, 1, 2, 3, 4])
        self.main_plot.tick_params(axis='y', which='both', left=True, right=False)
        genome_end = boundaries[-1]
        self.main_plot.set_xlim(0, genome_end)
        self.main_plot.margins(x=0)
        self.main_plot.set_axisbelow(True)
        # Optional baseline & bounds
        if baseline is not None:
            self.main_plot.plot([0, genome_end], [baseline, baseline],
                                linewidth=2.0, color="#3659FF", zorder=3, solid_capstyle="butt")
        if bounds is not None and len(bounds) == 2:
            lowerb, upperb = bounds[0], bounds[1]
            self.main_plot.plot([0, genome_end], [lowerb, lowerb],
                                linewidth=1.4, color="#EB4444", zorder=3)
            self.main_plot.plot([0, genome_end], [upperb, upperb],
                                linewidth=1.4, color="#EB4444", zorder=3)

        # Per-chromosome mean lines
        offset = 0
        for chrom in chromosomes:   
            length   = int(chr_lengths.get(chrom, 0))
            chr_mean = chr_means.get(chrom, np.nan)
            if length > 0 and np.isfinite(chr_mean):
                self.main_plot.plot([offset, offset + length], [chr_mean, chr_mean],
                                    linewidth=4.8, color="#585858", zorder=4, solid_capstyle="butt")
            offset += length

        # Chromosome boundaries
        for boundary in boundaries[:-1]:
            self.main_plot.axvline(boundary, color="#303030", linewidth=1.5, zorder=5)

        # === X axis labels at bottom ONLY (no numeric Mb ticks) ===
        # Convert labels like 'chr1' -> '1', 'chrX' -> 'X', etc.
        simple_labels = []
        for name in labels:
            lab = str(name)
            if lab.lower().startswith("chr"):
                lab = lab[3:]
            simple_labels.append(lab)

        # Clear any existing tick settings, then set our bottom chromosome labels
        self.main_plot.set_xticks(xticks)
        self.main_plot.set_xticklabels(simple_labels, rotation=0, ha="center", fontsize=23)
        # Remove the top axis (we're not using secondary_xaxis anymore)
        self.main_plot.tick_params(axis='x', which='both', top=False, labeltop=False)

        # Colorbar
        cbar = self.figure.colorbar(
            scatter, ax=self.main_plot, label="Aneuploidy",
            ticks=[0, 1, 2, 3, 4], fraction=0.03, pad=0.01
        )
        cbar.ax.tick_params(labelsize=15)
        cbar.set_label("Aneuploidy", fontsize=25)

        self.main_plot.set_ylabel("Aneuploidy", fontsize=25)
        self.main_plot.tick_params(axis="y", labelsize=20)
        self.figure.suptitle(self.file_prefix, fontsize=20)

        output_path = os.path.join(self.output_directory, "img", f"{self.file_prefix}_plot_cnv_genome.png")
        self.logger.info(f"Saving genome plot to: {output_path}")
        self.figure.savefig(output_path, dpi=800, bbox_inches="tight", transparent=False, facecolor="white")
        self.logger.info(f"Plot saved: {output_path}")
        self.figure.clf()
