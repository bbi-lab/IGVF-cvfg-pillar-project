import matplotlib as mpl

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = [
    'Arial',            # preferred
    'Helvetica',
    'Nimbus Sans',
    'DejaVu Sans'       # guaranteed fallback
]

import logging
logging.getLogger('matplotlib').setLevel(logging.ERROR)
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch
from matplotlib.colors import to_rgba, LinearSegmentedColormap, BoundaryNorm
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import pandas as pd
import json
# from src.assay_calibration.plot_utils.utils import sample_density, compute_classification_metrics
# from src.assay_calibration.data_utils.dataset import Scoreset
from data_utils import Scoreset
from scipy.stats import beta, truncnorm, cauchy, skewnorm

# Global styling constants - optimized for Nature print
FONTSIZE_PANEL_LETTER = 24
FONTSIZE_TITLE = 18
FONTSIZE_SUBTITLE = 12
FONTSIZE_AXIS_LABEL = 11
FONTSIZE_TICK = 9
FONTSIZE_LEGEND = 9
FONTSIZE_ANNOTATION = 10

TITLE_SPACER = 0.10
PANEL_BE_HSPACE = 0.22
PANEL_BE_HRATIOS = [1.5, 0.3, 0.3]#, 0.2]
PANEL_BE_TITLE_SPACER_FACTOR = 3

PANEL_LETTER_X = [0.01, 0.52]
PANEL_LETTER_Y = [0.9825, 0.795, 0.2675]

CMAP = ["#f0e5d0", "#f2f2f2", "#4a3d5f"]
WHITE_TEXT_THRESHOLD_4C = 0.7


SAMPLE_NAMES = ["ClinVar PLP", "ClinVar BLB", "gnomAD", "Synonymous"]
SAMPLE_COLORS = ['#CA7682', '#1D7AAB', '#A0A0A0', '#6BAA75']
SAMPLE_ALPHAS = [0.5, 0.5, 0.15, 0.4]

STRENGTH_COLOR = {
    -8: '#4b91a6', -7: '#5DA3BD', -6: '#6FAACE', -5: '#74ABCE',
    -4: '#7ab5d1', -3: '#99c8dc', -2: '#d0e8f0', -1: '#e4f1f6',
    0: '#e0e0e0',
    1: '#e6b1b8', 2: '#d68f99', 3: '#ca7682', 4: '#b85c6b',
    5: '#B1535F', 6: '#AA4E58', 7: '#A2484F', 8: '#943744'
}

BENIGN_THRESHOLD_COLOR = '#2166AC'
PATHOGENIC_THRESHOLD_COLOR = '#B2182B'
CARTOON_FIT_COLORS = ['#8B3A47', '#0D4A6B']

def joint_densities(x, params, weights):
    """
    weighted pdfs of a mixture of skew normal distributions
    """
    return np.array(
        [
            w * skewnorm.pdf(x, a, loc, scale)
            for (a, loc, scale), w in zip(params, weights)
        ]
    )


def sample_density(x, fits, sampleNum):
    _density = np.stack([joint_densities(x, _fit['fit']['component_params'],_fit['fit']['weights'][sampleNum])
                        for _fit in fits])
    density = np.full(_density.shape,np.nan)
    for fitIdx,fit in enumerate(fits):
        fit_xmin,fit_xmax = fit['fit']['xlims']
        mask = (x >= fit_xmin) & (x <= fit_xmax)
        density[fitIdx,:,mask] = _density[fitIdx,:,mask]
    return density
    

def compute_classification_metrics(confusion_matrix):
    """
    Compute classification performance metrics from a confusion matrix.
    
    Parameters:
    -----------
    confusion_matrix : DataFrame
        Confusion matrix with rows = ClinVar (B/LB, P/LP)
                         and columns = Evidence categories (Benign, Indeterminate, Pathogenic)
                         OR Functional categories (Normal, Indeterminate, Abnormal)
    
    Returns:
    --------
    dict with metrics: accuracy, sensitivity, specificity, MCC, DOR (determinate only),
                       LR+ and DOR for dichotomized predictions (including uncertain)
    """
    
    # Get counts
    if confusion_matrix.shape[0] != 2 or confusion_matrix.shape[1] != 3:
        raise ValueError(f"Expected 2×3 confusion matrix, got {confusion_matrix.shape}")
    
    # Row 0 = B/LB (negatives), Row 1 = P/LP (positives)
    # Col 0 = Benign/Normal, Col 1 = Indeterminate, Col 2 = Pathogenic/Abnormal
    
    TN = confusion_matrix.iloc[0, 0]  # B/LB → Benign/Normal
    FP = confusion_matrix.iloc[0, 2]  # B/LB → Pathogenic/Abnormal
    FN = confusion_matrix.iloc[1, 0]  # P/LP → Benign/Normal
    TP = confusion_matrix.iloc[1, 2]  # P/LP → Pathogenic/Abnormal
    
    uncertain_benign = confusion_matrix.iloc[0, 1]  # B/LB → Indeterminate
    uncertain_path = confusion_matrix.iloc[1, 1]    # P/LP → Indeterminate
    
    total_variants = confusion_matrix.sum().sum()
    determinate_calls = TP + TN + FP + FN
    
    # Compute metrics
    coverage = determinate_calls / total_variants if total_variants > 0 else 0.0
    
    # ========================================================================
    # STANDARD METRICS (determinate calls only, excluding uncertain)
    # ========================================================================
    if determinate_calls > 0:
        accuracy = (TP + TN) / determinate_calls
        sensitivity = TP / (TP + FN) if (TP + FN) > 0 else 0.0
        specificity = TN / (TN + FP) if (TN + FP) > 0 else 0.0
        
        # Matthews Correlation Coefficient
        numerator = (TP * TN) - (FP * FN)
        denominator = np.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
        mcc = numerator / denominator if denominator > 0 else 0.0
        
        # Standard DOR (on determinate calls only)
        dor_standard = (TP * TN) / (FP * FN) if (FP > 0 and FN > 0) else float('inf')
        
        # Standard LR+ (on determinate calls only)
        lr_plus_standard = sensitivity / (1 - specificity) if specificity < 1.0 else float('inf')
        
    else:
        accuracy = sensitivity = specificity = mcc = 0.0
        dor_standard = lr_plus_standard = 0.0
    
    # ========================================================================
    # DICHOTOMIZED METRICS (including uncertain variants)
    # ========================================================================
    
    # For PATHOGENIC prediction: Pathogenic vs. non-Pathogenic (Benign + Uncertain)
    TP_path_dich = TP                                    # P/LP → Pathogenic
    FN_path_dich = FN + uncertain_path                   # P/LP → (Benign or Uncertain)
    FP_path_dich = FP                                    # B/LB → Pathogenic
    TN_path_dich = TN + uncertain_benign                 # B/LB → (Benign or Uncertain)
    
    sens_path_dich = TP_path_dich / (TP_path_dich + FN_path_dich) if (TP_path_dich + FN_path_dich) > 0 else 0.0
    spec_path_dich = TN_path_dich / (TN_path_dich + FP_path_dich) if (TN_path_dich + FP_path_dich) > 0 else 0.0
    
    lr_plus_pathogenic = sens_path_dich / (1 - spec_path_dich) if spec_path_dich < 1.0 else float('inf')
    dor_pathogenic = (TP_path_dich * TN_path_dich) / (FP_path_dich * FN_path_dich) if (FP_path_dich > 0 and FN_path_dich > 0) else float('inf')
    
    # For BENIGN prediction: Benign vs. non-Benign (Pathogenic + Uncertain)
    TP_ben_dich = TN                                     # B/LB → Benign
    FN_ben_dich = FP + uncertain_benign                  # B/LB → (Pathogenic or Uncertain)
    FP_ben_dich = FN                                     # P/LP → Benign
    TN_ben_dich = TP + uncertain_path                    # P/LP → (Pathogenic or Uncertain)
    
    sens_ben_dich = TP_ben_dich / (TP_ben_dich + FN_ben_dich) if (TP_ben_dich + FN_ben_dich) > 0 else 0.0
    spec_ben_dich = TN_ben_dich / (TN_ben_dich + FP_ben_dich) if (TN_ben_dich + FP_ben_dich) > 0 else 0.0
    
    lr_plus_benign = sens_ben_dich / (1 - spec_ben_dich) if spec_ben_dich < 1.0 else float('inf')
    dor_benign = (TP_ben_dich * TN_ben_dich) / (FP_ben_dich * FN_ben_dich) if (FP_ben_dich > 0 and FN_ben_dich > 0) else float('inf')
    
    return {
        'TP': int(TP),
        'TN': int(TN),
        'FP': int(FP),
        'FN': int(FN),
        'uncertain': int(uncertain_benign + uncertain_path),
        'total': int(total_variants),
        'determinate': int(determinate_calls),
        'coverage': float(coverage),
        'accuracy': float(accuracy),
        'sensitivity': float(sensitivity),
        'specificity': float(specificity),
        'mcc': float(mcc),
        'dor_standard': float(dor_standard),
        'lr_plus_standard': float(lr_plus_standard),
        'lr_plus_pathogenic': float(lr_plus_pathogenic),
        'lr_plus_benign': float(lr_plus_benign),
        'dor_pathogenic': float(dor_pathogenic),
        'dor_benign': float(dor_benign),
    }

def plot_panel_letters(fig, letters):
    """Plots panel letters on a figure using fig.text."""
    
    pos_idx = [(i, j) for i in range(len(PANEL_LETTER_X)) 
                      for j in range(len(PANEL_LETTER_Y))]
    
    for letter, (ix, iy) in zip(letters, pos_idx):
        x = PANEL_LETTER_X[ix]
        y = PANEL_LETTER_Y[iy]
        fig.text(x, y, letter, fontsize=FONTSIZE_PANEL_LETTER,
                 fontweight='bold', va='top', ha='left')


def plot_panel_a(gs_spec, scoreset_2018, indv_summary, fits, score_range, flipped, fig):
    """Panel A: Multi-sample mixture model calibration"""
    # Add title row
    gs_outer = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_spec, 
                                                height_ratios=[0.08, 1], hspace=TITLE_SPACER)
    gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs_outer[1], wspace=0.10)
    
    x_min, x_max = score_range[0], score_range[-1]
    bin_width = (x_max - x_min) / 50
    point_ranges = indv_summary['point_ranges']
    
    # Pre-compute thresholds
    linestyles = ['dotted', 'dashed', 'dashdot', (5, (10, 3)), (0, (3, 5, 1, 5))]
    threshold_info = []
    
    for idx, point_val in enumerate([1, 2, 3, 4, 8]):
        for pv, score_ranges_pr in point_ranges.items():
            if int(pv) == -point_val and score_ranges_pr:
                threshold_score = score_ranges_pr[0][0] if not flipped else score_ranges_pr[0][1]
                threshold_info.append((int(pv), threshold_score, BENIGN_THRESHOLD_COLOR, linestyles[idx], 1.0))
                break
        for pv, score_ranges_pr in point_ranges.items():
            if int(pv) == point_val and score_ranges_pr:
                threshold_score = score_ranges_pr[0][1] if not flipped else score_ranges_pr[0][0]
                threshold_info.append((int(pv), threshold_score, PATHOGENIC_THRESHOLD_COLOR, linestyles[idx], 1.0))
                break
    
    num_skipped = 0
    
    for ax_idx,sample_num in enumerate([1,0,2]):
        if scoreset_2018.sample_counts[sample_num] == 0:
            num_skipped += 1
            continue
        
        sample_idx = sample_num - num_skipped
        ax = plt.subplot(gs[ax_idx])
        
        sample_mask = scoreset_2018.sample_assignments[:, sample_idx]
        color = SAMPLE_COLORS[sample_num]
        hist_data = scoreset_2018.scores[sample_mask]
        n_count = sample_mask.sum()
        
        sns.histplot(hist_data, binwidth=bin_width, stat='density', ax=ax,
                   alpha=0.5, color=color)
        
        density_sample = sample_density(score_range, fits, sample_idx)
        d_total = np.nansum(density_sample, axis=1)
        d_total_perc = np.percentile(d_total, [5, 50, 95], axis=0)
        
        ax.fill_between(score_range, d_total_perc[0], d_total_perc[2], color='gray', alpha=0.3)
        ax.plot(score_range, d_total_perc[1], color='black', alpha=0.65, linewidth=2.5)
        
        ax.set_xlim(x_min, x_max)
        max_hist_density = max([patch.get_height() for patch in ax.patches]) if ax.patches else 1.0
        ax.set_ylim(0, 1.1*max_hist_density)
        ax.set_xlabel('Assay score', fontsize=FONTSIZE_AXIS_LABEL)
        ax.set_ylabel('')#'Density' if ax_idx == 0 else '', fontsize=FONTSIZE_AXIS_LABEL)
        ax.tick_params(axis='both', labelsize=FONTSIZE_TICK, left=False, labelleft=False)
        
        for pv, thresh_score, thresh_color, thresh_ls, thresh_lw in threshold_info:
            if abs(pv) in [1, 2, 4, 8]:
                ax.axvline(thresh_score, color=thresh_color, linestyle=thresh_ls, 
                          linewidth=1.5, alpha=0.8)
        
        face_rgba = to_rgba(color, 0.5)
        hist_patch = Patch(facecolor=face_rgba, edgecolor='black')
        
        if sample_num == 2:
            legend_label = f'{SAMPLE_NAMES[sample_idx]}\nprior: {indv_summary["prior"]:.3f}\n(n={n_count:,})'
        else:
            legend_label = f'{SAMPLE_NAMES[sample_idx]}\n(n={n_count:,})'
        
        loc = 'upper left' if (sample_num == 0 and flipped) or (sample_num != 0 and not flipped) else 'upper right'
        ax.legend([hist_patch], [legend_label], loc=loc, fontsize=FONTSIZE_LEGEND, framealpha=0.9)
        
        # if sample_idx == 1:
        #     ax.text(-0.25, 1.15, "a", transform=ax.transAxes,
        #            fontsize=FONTSIZE_PANEL_LETTER, fontweight='bold', va='top', ha='left')
    
    # Add centered title in title row
    bbox = plt.subplot(gs_outer[0])
    bbox.axis('off')
    bbox.text(0.5, 2.5, 'ExCALIBR experimental data calibration',
             transform=bbox.transAxes, ha='center', va='center',
             fontsize=FONTSIZE_TITLE, fontweight='bold')


def plot_panel_b(gs_spec, scoreset, all_scores, point_ranges, score_range, flipped, fig):
    """Panel B: Experimental score calibration comparison"""
    # Add title row
    gs_outer = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_spec,
                                                height_ratios=[0.06, 1], hspace=TITLE_SPACER-PANEL_BE_HSPACE/PANEL_BE_TITLE_SPACER_FACTOR)
    gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_outer[1],
                                         height_ratios=PANEL_BE_HRATIOS, hspace=PANEL_BE_HSPACE)
    
    x_min, x_max = score_range[0], score_range[-1]
    bin_width = (x_max - x_min) / 50
    
    # Histogram
    ax_hist = plt.subplot(gs[0])
    ax_twin = ax_hist.twinx()  # Create twin axis for SNV counts
    
    sample_handles = []
    
    num_skipped = 0
    for sample_num in [1,0,2]:
        if scoreset.sample_counts[sample_num] == 0:
            num_skipped += 1
            continue
        
        sample_idx = sample_num - num_skipped
        sample_mask = scoreset.sample_assignments[:, sample_idx]
        color = SAMPLE_COLORS[sample_num]
        alpha = SAMPLE_ALPHAS[sample_num]
        
        if sample_num == 2:
            hist_data = all_scores
            display_name = 'All SNVs'
            n_count = len(all_scores)
            # Plot All SNVs on twin axis (right)
            sns.histplot(hist_data, binwidth=bin_width, stat='count', ax=ax_twin,
                       alpha=alpha, color=color)
        else:
            hist_data = scoreset.scores[sample_mask]
            display_name = SAMPLE_NAMES[sample_num]
            n_count = sample_mask.sum()
            # Plot PLP/BLB on main axis (left)
            sns.histplot(hist_data, binwidth=bin_width, stat='count', ax=ax_hist,
                       alpha=alpha, color=color)
        
        face_rgba = to_rgba(color, alpha)
        hist_patch = Patch(facecolor=face_rgba, edgecolor='black')
        sample_handles.append((hist_patch, f'{display_name}\n(n={n_count:,})'))
    
    ax_hist.set_xlim(x_min, x_max)
    ax_twin.set_xlim(x_min, x_max)
    ax_hist.set_ylim(0,  1.18*max([patch.get_height() for patch in ax_hist.patches]) if ax_hist.patches else 1.0)
    ax_twin.set_ylim(0,  1.18*max([patch.get_height() for patch in ax_twin.patches]) if ax_twin.patches else 1.0)
    
    ax_hist.set_xlabel('')
    ax_hist.set_ylabel('Control variant count', fontsize=FONTSIZE_AXIS_LABEL)
    ax_twin.set_ylabel('SNV count', fontsize=FONTSIZE_AXIS_LABEL)
    
    ax_hist.tick_params(axis='both', labelsize=FONTSIZE_TICK)
    ax_twin.tick_params(axis='both', labelsize=FONTSIZE_TICK)
    
    # Combine legends from both axes
    ax_hist.legend([h[0] for h in sample_handles], [h[1] for h in sample_handles],
                  loc='upper right', fontsize=FONTSIZE_LEGEND)

    ax_hist.text(0.5, 0.98, rf'$\mathbfit{{MSH2}}$'+' experimental scores',
            transform=ax_hist.transAxes,
            fontsize=FONTSIZE_SUBTITLE, 
            fontweight='bold',
            va='top', 
            ha='center')
    
    # Add centered title in title row
    bbox = plt.subplot(gs_outer[0])
    bbox.axis('off')
    
    # Scott et al. bar
    ax_scott = plt.subplot(gs[1])
    ax_scott.axvspan(x_min, 0, color=STRENGTH_COLOR[-4], alpha=1.0)
    ax_scott.axvspan(0, 0.4, color=STRENGTH_COLOR[0], alpha=1.0)
    ax_scott.axvspan(0.4, x_max, color=STRENGTH_COLOR[4], alpha=1.0)
    
    count_below_0 = (all_scores < 0).sum()
    count_0_to_04 = ((all_scores >= 0) & (all_scores < 0.4)).sum()
    count_above_04 = (all_scores >= 0.4).sum()
    ax_scott.text((x_min + 0) / 2, 0.5, f'{count_below_0:,}', ha='center', va='center', fontsize=FONTSIZE_ANNOTATION)
    ax_scott.text((0 + 0.4) / 2, 0.5, f'{count_0_to_04:,}', ha='center', va='center', fontsize=FONTSIZE_ANNOTATION)
    ax_scott.text((0.4 + x_max) / 2, 0.5, f'{count_above_04:,}', ha='center', va='center', fontsize=FONTSIZE_ANNOTATION)
    
    ax_scott.set_xlim(x_min, x_max)
    ax_scott.set_ylim(0, 1)
    ax_scott.set_yticks([])
    ax_scott.set_xticks([])
    ax_scott.set_title('Scott et al. (2022)', loc='left', pad=3, fontsize=FONTSIZE_SUBTITLE, style='italic')
    
    # ExCALIBR bar
    ax_excalibr = plt.subplot(gs[2])
    
    intervals = []
    for pv in sorted([int(p) for p in point_ranges.keys() if p != 0]):
        if point_ranges[pv]:
            sr = point_ranges[pv][0]
            intervals.append((pv, sr[0], sr[1]))
    
    neg_int = [(pv, s, e) for pv, s, e in intervals if pv < 0]
    pos_int = [(pv, s, e) for pv, s, e in intervals if pv > 0]
    
    if neg_int and pos_int:
        neg_sorted = sorted(neg_int, key=lambda x: x[2])
        pos_sorted = sorted(pos_int, key=lambda x: x[1])
        ir_start = neg_sorted[-1][2] if flipped else pos_sorted[-1][2]
        ir_end = pos_sorted[0][1] if flipped else neg_sorted[0][1]
        intervals.append((0, ir_start, ir_end))
    
    intervals_sorted = sorted(intervals, key=lambda x: x[1])
    
    for point_val, start, end in intervals_sorted:
        ax_excalibr.axvspan(start, end, color=STRENGTH_COLOR[point_val], alpha=1.0)
        count = ((all_scores >= start) & (all_scores < end)).sum()
        if (end - start) > 0.3:
            text_color = 'white' if abs(point_val) >= 7 else 'black'
            ax_excalibr.text((start + end) / 2, 0.5, f'{count:,}',
                           ha='center', va='center', fontsize=FONTSIZE_ANNOTATION, color=text_color)
    
    ax_excalibr.set_xlim(x_min, x_max)
    ax_excalibr.set_ylim(0, 1)
    ax_excalibr.set_yticks([])
    ax_excalibr.set_xlabel(rf'$\mathit{{MSH2}}$ experimental score', fontsize=FONTSIZE_AXIS_LABEL)
    ax_excalibr.tick_params(axis='x', labelsize=FONTSIZE_TICK)
    ax_excalibr.set_title('ExCALIBR', loc='left', pad=3, fontsize=FONTSIZE_SUBTITLE, style='italic')
    
    legend_order = [-8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8] if flipped else [8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8]
    
    point_labels = {
        -8: "-8 (very strong)", -4: "-4 (strong)", -3: "-3", -2: "-2 (moderate)", -1: "-1 (supporting)",
        0: "0 (indeterminate)",
        1: "+1 (supporting)", 2: "+2 (moderate)", 3: "+3", 4: "+4 (strong)", 8: "+8 (very strong)"
    }

    scott_legend_info = [(-4,None,None),(0,None,None),(4,None,None)]
    legend_handles = [Patch(facecolor=STRENGTH_COLOR[pv], label=point_labels.get(pv, f"{pv:+d}"), edgecolor='none')
                     for pv in legend_order if any(p == pv for p, _, _ in intervals_sorted+scott_legend_info)]

    return legend_handles


def plot_panel_c(gs_spec, danzs_oob, auths_oob, fig):
    """Panel C: Confusion matrices with purple gradient"""
    # Add title row
    gs_outer = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_spec,
                                                height_ratios=[0.00, 1], hspace=0)
    gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs_outer[1], wspace=0.15)
    
    danz_agg = sum([d for d in danzs_oob if d is not None])
    auth_agg = sum([a for a in auths_oob if a is not None])
    danz_metrics = compute_classification_metrics(danz_agg)
    auth_metrics = compute_classification_metrics(auth_agg)
    
    # Purple gradient colormap
    purple_cmap = LinearSegmentedColormap.from_list("purple", CMAP[1:])
    max_val = max(danz_agg.values.max(), auth_agg.values.max())
    
    def plot_confusion(df, ax, title, metrics, show_cbar=False, cbar_ax=None):
        def get_text_color(value, max_value):
            return 'white' if value / max_value > WHITE_TEXT_THRESHOLD_4C else 'black'
        
        annot = df.copy().astype(str)
        for row in range(len(df)):
            for col in range(len(df.columns)):
                annot.iloc[row, col] = f"{df.iloc[row, col]:,}"
        
        sns.heatmap(df, annot=annot, fmt='', cmap=purple_cmap, vmin=0, vmax=max_val,
                   ax=ax, cbar=show_cbar, cbar_ax=cbar_ax,
                   cbar_kws=None,#{'label': 'Count'} if show_cbar else None,
                   linewidths=3.0, linecolor='white',
                   annot_kws={'fontsize': FONTSIZE_ANNOTATION + 1, 'ha': 'center', 'va': 'center'}, square=True)
        
        for text_obj in ax.texts:
            x, y = text_obj.get_position()
            row, col = int(y), int(x)
            if row < len(df) and col < len(df.columns):
                value = df.iloc[row, col]
                text_obj.set_color(get_text_color(value, max_val))
        
        ax.set_facecolor('#F9F9F9')
        
        xlabels = ["Benign", "Indeterminate", "Pathogenic"] if not show_cbar else ["Normal", "Indeterminate", "Abnormal"]
        ylabels = ["BLB", "PLP"]
        
        xlabel_text = 'Evidence direction' if not show_cbar else 'Functional class'
        ax.set_xlabel(xlabel_text, fontsize=FONTSIZE_AXIS_LABEL)
        ax.set_ylabel('ClinVar Classification' if title.startswith("ExCALIBR") else '', fontsize=FONTSIZE_AXIS_LABEL)
        ax.set_xticklabels(xlabels, rotation=0, fontsize=FONTSIZE_TICK)
        ax.set_yticklabels(ylabels, rotation=0, fontsize=FONTSIZE_TICK)
        ax.set_title(title, fontsize=FONTSIZE_SUBTITLE, fontweight='bold', pad=6)
        ax.tick_params(length=0, labelsize=FONTSIZE_TICK)
        
        coverage_text = f"DOR: {metrics['dor_standard']:.1f}\nDeterminate: Controls {100*metrics['coverage']:.1f}%, VUS "
        coverage_text += "79.7%" if "ExCALIBR" in title else "93.2%"
        
        ax.text(0.5, -0.26, coverage_text, transform=ax.transAxes,
               fontsize=FONTSIZE_LEGEND, ha='center', va='top', color='#555555')
        
        for spine in ax.spines.values():
            spine.set_visible(False)
    
    ax_danz = plt.subplot(gs[0])
    ax_auth = plt.subplot(gs[1])
    ax_danz.set_aspect('equal', adjustable='box')
    ax_auth.set_aspect('equal', adjustable='box')
    
    # Create position for colorbar manually (after axes created)
    fig.canvas.draw()
    bbox_auth = ax_auth.get_position()
    cbar_ax = fig.add_axes([bbox_auth.x1 + 0.01, bbox_auth.y0, 0.015, bbox_auth.height])
    
    plot_confusion(danz_agg, ax_danz, "ExCALIBR Evidence", danz_metrics, show_cbar=False)
    plot_confusion(auth_agg, ax_auth, "Author Annotations", auth_metrics, show_cbar=True, cbar_ax=cbar_ax)
    
    # ax_danz.text(-0.25, 1.15, "c", transform=ax_danz.transAxes,
    #             fontsize=FONTSIZE_PANEL_LETTER, fontweight='bold', va='top', ha='left')
    
    # Add centered title in title row
    bbox_title = plt.subplot(gs_outer[0])
    bbox_title.axis('off')
    # bbox_title.text(0.5, 0.5, 'ExCALIBR calibration vs. author functional class',
    #                transform=bbox_title.transAxes, ha='center', va='center',
    #                fontsize=FONTSIZE_TITLE, fontweight='bold')


def plot_panel_d(gs_spec, prior, Post_p, Post_b, p_data_sim, b_data_sim, fig):
    """Panel D: Simulation-based optimization"""
    # Add title row
    gs_outer = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_spec,
                                                height_ratios=[0.08, 1], hspace=TITLE_SPACER)
    gs = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs_outer[1], 
                                         width_ratios=[1, 0.05, 0.5, 0.5], wspace=0.4)
    
    # Add centered title in title row
    bbox = plt.subplot(gs_outer[0])
    bbox.axis('off')
    bbox.text(0.5, 2.5, 'Single-gene predictor calibration',
             transform=bbox.transAxes, ha='center', va='center',
             fontsize=FONTSIZE_TITLE, fontweight='bold')
    
    # Simulation fitting
    ax_sim = plt.subplot(gs[0])
    
    def add_fits(ax, data, color):
        """Add fits using threshold line colors"""
        x_fit = np.linspace(0, 1, 500)
        a_beta, b_beta, loc, scale = beta.fit(data, floc=0, fscale=1)
        ax.plot(x_fit, beta.pdf(x_fit, a_beta, b_beta, loc=loc, scale=scale),
               linestyle='-', color=color, lw=2.5, alpha=0.8)
        mu, sigma = np.mean(data), np.std(data)
        a_trunc, b_trunc = (0 - mu) / sigma, (1 - mu) / sigma
        ax.plot(x_fit, truncnorm.pdf(x_fit, a_trunc, b_trunc, loc=mu, scale=sigma),
               linestyle='--', color=color, lw=2.5, alpha=0.8)
        alpha_skew, loc_skew, scale_skew = skewnorm.fit(data)
        ax.plot(x_fit, skewnorm.pdf(x_fit, alpha_skew, loc=loc_skew, scale=scale_skew),
               linestyle='-.', color=color, lw=2.5, alpha=0.8)
        loc_c, scale_c = np.median(data), np.std(data)
        ax.plot(x_fit, cauchy.pdf(x_fit, loc=loc_c, scale=scale_c),
               linestyle=':', color=color, lw=2.5, alpha=0.8)
    
    # Use same colors and alphas as other histograms
    ax_sim.hist(p_data_sim, bins=30, density=True, alpha=SAMPLE_ALPHAS[0], 
               color=SAMPLE_COLORS[0], edgecolor="black")
    ax_sim.hist(b_data_sim, bins=30, density=True, alpha=SAMPLE_ALPHAS[1], 
               color=SAMPLE_COLORS[1], edgecolor="black")
    
    # Fit lines in threshold colors (blue for benign, red for pathogenic)
    add_fits(ax_sim, p_data_sim, CARTOON_FIT_COLORS[0])  # Red for pathogenic
    add_fits(ax_sim, b_data_sim, CARTOON_FIT_COLORS[1])  # Blue for benign
    
    ax_sim.set_xlabel("Prediction score", fontsize=FONTSIZE_AXIS_LABEL)
    ax_sim.set_ylabel('')#"Density", fontsize=FONTSIZE_AXIS_LABEL)
    ax_sim.tick_params(labelsize=FONTSIZE_TICK, left=False, labelleft=False)
    
    # Legend with same style as other panels
    face_p = to_rgba(SAMPLE_COLORS[0], SAMPLE_ALPHAS[0])
    face_b = to_rgba(SAMPLE_COLORS[1], SAMPLE_ALPHAS[1])
    custom_patches = [
        Patch(facecolor=face_b, edgecolor='black', label="ClinVar BLB"),
        Patch(facecolor=face_p, edgecolor='black', label="ClinVar PLP")
    ]
    ax_sim.legend(handles=custom_patches, fontsize=FONTSIZE_LEGEND, frameon=True, 
                 framealpha=0.9, loc='upper left')

    max_hist_density = max([patch.get_height() for patch in ax_sim.patches]) if ax_sim.patches else 1.0
    ax_sim.set_ylim(0, 1.1*max_hist_density)
    
    # Posterior plots
    x = np.linspace(0, 1, 100)
    y_para = 1 - x**2
    y_exp = np.exp(3 * x) / np.exp(3)
    
    ax_benign = plt.subplot(gs[2])
    ax_benign.plot(x, y_para, linewidth=3.5, color='#13506b')
    
    linestyles_post = ['dotted', 'dashed', 'dashdot', (5, (10, 3)), 'solid']
    for i in range(5):
        ax_benign.axhline(Post_b[(4-i)], linestyle=linestyles_post[i], color=BENIGN_THRESHOLD_COLOR, linewidth=2.5)
        # if i < 3:
        ax_benign.text(0.95, Post_b[(4-i)] - 0.0003, f'-{i+1}' if i < 4 else "-8", ha='center', va='top',
                          color=BENIGN_THRESHOLD_COLOR, fontsize=FONTSIZE_ANNOTATION, fontweight='bold')
    
    ax_benign.set_ylim(0.8, 1.0)
    ax_benign.set_xlabel('Score', fontsize=FONTSIZE_AXIS_LABEL)
    ax_benign.set_ylabel('Posterior', fontsize=FONTSIZE_AXIS_LABEL)
    # ax_benign.set_title('Benign', color='#1D7AAB', fontweight='bold', fontsize=FONTSIZE_SUBTITLE)
    ax_benign.tick_params(labelsize=FONTSIZE_TICK)
    
    ax_path = plt.subplot(gs[3])
    ax_path.plot(x, y_exp, linewidth=3.5, color='#7a2f37')
    
    for i in range(5):
        ax_path.axhline(Post_p[(4-i)], linestyle=linestyles_post[i], color=PATHOGENIC_THRESHOLD_COLOR, linewidth=2.5)
        label = f'+{i+1}' if i < 3 else ('+4' if i == 3 else '+8')
        ax_path.text(0.05, Post_p[(4-i)] - 0.005, label, ha='center', va='top',
                    color=PATHOGENIC_THRESHOLD_COLOR, fontsize=FONTSIZE_ANNOTATION, fontweight='bold')
    
    ax_path.set_ylim(0, 1.001)
    ax_path.set_xlabel('Score', fontsize=FONTSIZE_AXIS_LABEL)
    # ax_path.set_title('Pathogenic', color='#943744', fontweight='bold', fontsize=FONTSIZE_SUBTITLE)
    ax_path.tick_params(labelsize=FONTSIZE_TICK)
    
    # ax_sim.text(-0.25, 1.15, "d", transform=ax_sim.transAxes,
    #            fontsize=FONTSIZE_PANEL_LETTER, fontweight='bold', va='top', ha='left')


def plot_panel_e(gs_spec, gene, dist, labdat, snvdf, sorted_thresholds, oldsorted_thresholds, fig):
    """Panel E: REVEL comparison"""
    strenth_to_point = {
        "BP4_Very Strong": -8, "BP4_Strong": -4, "BP4_Moderate+": -3,
        "BP4_Moderate": -2, "BP4_Supporting": -1, "IR": 0,
        "PP3_Supporting": 1, "PP3_Moderate": 2, "PP3_Moderate+": 3,
        "PP3_Strong": 4, "PP3_Very Strong": 8
    }
    
    # Add title row
    gs_outer = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_spec,
                                                height_ratios=[0.06, 1], hspace=TITLE_SPACER-PANEL_BE_HSPACE/PANEL_BE_TITLE_SPACER_FACTOR)
    gs = gridspec.GridSpecFromSubplotSpec(3, 1, subplot_spec=gs_outer[1],
                                         height_ratios=PANEL_BE_HRATIOS, hspace=PANEL_BE_HSPACE)
    
    ax_hist = plt.subplot(gs[0])
    ax_twin = ax_hist.twinx()
    
    bin_width = (labdat[0].max() - labdat[0].min()) / 50
    sns.histplot(labdat[labdat[1] == 0][0], binwidth=bin_width, color=SAMPLE_COLORS[1],
                alpha=SAMPLE_ALPHAS[1], ax=ax_hist, label=f'ClinVar BLB\n(n={len(labdat[labdat[1] == 0])})')
    sns.histplot(labdat[labdat[1] == 1][0], binwidth=bin_width, color=SAMPLE_COLORS[0],
                alpha=SAMPLE_ALPHAS[0], ax=ax_hist, label=f'ClinVar PLP\n(n={len(labdat[labdat[1] == 1])})')
    sns.histplot(snvdf[dist], binwidth=bin_width, color=SAMPLE_COLORS[2], alpha=SAMPLE_ALPHAS[2],
                ax=ax_twin, label=f'All SNVs\n(n={len(snvdf):,})')
    
    ax_hist.set_xlim(0, 1)
    ax_hist.set_ylim(0,  1.18*max([patch.get_height() for patch in ax_hist.patches]) if ax_hist.patches else 1.0)
    ax_twin.set_ylim(0,  1.18*max([patch.get_height() for patch in ax_twin.patches]) if ax_twin.patches else 1.0)
    ax_hist.set_xlabel('')
    ax_hist.set_ylabel('Control variant count', fontsize=FONTSIZE_AXIS_LABEL)
    ax_twin.set_ylabel('SNV count', fontsize=FONTSIZE_AXIS_LABEL)
    ax_hist.tick_params(labelsize=FONTSIZE_TICK)
    ax_twin.tick_params(labelsize=FONTSIZE_TICK)

    ax_hist.text(0.5, 0.98, rf'$\mathbfit{{MSH2}}$'+' REVEL scores',
            transform=ax_hist.transAxes,
            fontsize=FONTSIZE_SUBTITLE, 
            fontweight='bold',
            va='top', 
            ha='center')  # Centered at top
    
    lines1, labels1 = ax_hist.get_legend_handles_labels()
    lines2, labels2 = ax_twin.get_legend_handles_labels()
    ax_hist.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=FONTSIZE_LEGEND)
    ax_twin.get_legend().remove() if ax_twin.get_legend() else None
    
    # ax_hist.text(-0.05, 1.15, "e", transform=ax_hist.transAxes,
    #             fontsize=FONTSIZE_PANEL_LETTER, fontweight='bold', va='top', ha='left')
    
    # Add centered title in title row
    bbox = plt.subplot(gs_outer[0])
    bbox.axis('off')
    # bbox.text(0.5, 0.5, rf'Comparison of $\mathbfit{{{gene}}}$ {dist} score calibration',
    #          transform=bbox.transAxes, ha='center', va='center',
    #          fontsize=FONTSIZE_TITLE, fontweight='bold')
    
    # Stacked bars
    def create_stacked_bar(ax, thresholds, data, title):
        TEXT_THRESH = 0.05
        trimmed_thresh = thresholds[1:]
        trimmed_thresh_keys = list(trimmed_thresh.keys())
        for i, (index, thresh) in enumerate(trimmed_thresh.items()):
            if 'BP4' in index and (i == 0 or np.isnan(thresholds.iloc[i])):
                start, end = 0, thresh
                xpos = thresh / 2
                threshdiff = thresh
            elif index == 'PP3_Very Strong' or (i == len(trimmed_thresh)-1 and 'PP3' in index):
                start, end = thresholds.iloc[i], 1
                xpos = (thresholds.iloc[i] + 1) / 2
                threshdiff = 1 - thresholds.iloc[i]
            elif index == 'PP3_Supporting':
                start, end = thresholds.iloc[i], thresh
                ax.axvspan(start, end, color=STRENGTH_COLOR[strenth_to_point['IR']], alpha=1)
                xpos = (start + end) / 2
                threshdiff = end - start
                mask = (data >= start) & (data < end)
                count = mask.sum()
                if count > 3 and threshdiff > TEXT_THRESH:
                    ax.text(xpos, 0.5, f'{count:,}', ha='center', va='center', fontsize=FONTSIZE_ANNOTATION)
                continue
            else:
                start, end = thresholds.iloc[i], thresh
                xpos = (start + end) / 2
                threshdiff = end - start

            if index.startswith('PP3'):
                ax.axvspan(start, end, color=STRENGTH_COLOR[strenth_to_point[trimmed_thresh_keys[i-1]]], alpha=1)
            else:
                ax.axvspan(start, end, color=STRENGTH_COLOR[strenth_to_point[index]], alpha=1)
            
            mask = (data >= start) & (data < end)
            count = mask.sum()
            if count > 3 and threshdiff > TEXT_THRESH:
                ax.text(xpos, 0.5, f'{count:,}', ha='center', va='center', fontsize=FONTSIZE_ANNOTATION)
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        if "single" in title.lower():
            ax.tick_params(axis='x', labelsize=FONTSIZE_TICK)
        else:
            ax.set_xticks([])
        ax.set_title(title, loc='left', pad=3, fontsize=FONTSIZE_SUBTITLE, style='italic')
    
    ax_agg = plt.subplot(gs[1])
    ax_single = plt.subplot(gs[2])
    
    create_stacked_bar(ax_agg, oldsorted_thresholds, snvdf[dist].values, 'Gene aggregation')
    create_stacked_bar(ax_single, sorted_thresholds, snvdf[dist].values, 'Single gene')
    ax_single.set_xlabel(rf'$\mathit{{MSH2}}$ {dist} score', fontsize=FONTSIZE_AXIS_LABEL)
    
    # Legend - two rows
    # ax_leg = plt.subplot(gs[3])
    # ax_leg.axis('off')
    strens = ['-4 (strong)', '-3', '-2 (moderate)', '-1 (supporting)', '0 (indeterminate)',
              '+1 (supporting)', '+2 (moderate)', '+3', '+4 (strong)']
    fill_colors = ['#7ab5d1', '#99c8dc', '#d0e8f0', '#e4f1f6', '#E0E0E0',
                  '#e6b1b8', '#d68f99', '#ca7682', '#b85c6b']
    legend_elements = [Patch(facecolor=color, label=label) for color, label in zip(fill_colors, strens)]
    ncol = (len(legend_elements) + 1) // 2
    # ax_leg.legend(handles=legend_elements, loc='center', ncol=ncol, frameon=False,
    #              fontsize=FONTSIZE_LEGEND-1, columnspacing=0.6, handletextpad=0.4, handlelength=0.8)
    return legend_elements


def plot_panel_f(gs_spec, dist, finalout, fig):
    """Panel F: Calibration differences heatmap"""
    # Add title row
    gs_outer = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs_spec,
                                                height_ratios=[0.00, 1], hspace=0)
    
    ax = plt.subplot(gs_outer[1])
    
    classifications = ['BLB', 'PLP', 'VUS', 'gnomAD', 'allSNVs']
    score_levels = [-4, -3, -2, -1, 0, 1, 2, 3, 4]
    
    def compute_percentage(df, score_col):
        table = df.pivot_table(index='merg_clinvar_sig', columns=score_col,
                               aggfunc='size', fill_value=0)
        table = table.div(table.sum(axis=1), axis=0) * 100
        return table.reindex(index=classifications, columns=score_levels, fill_value=0)
    
    old_data = compute_percentage(finalout, 'Old_scr')
    new_data = compute_percentage(finalout, 'New_scr')
    diff_data = new_data - old_data
    
    annot_data = diff_data.map(lambda x: f"{x:+.1f}%" if abs(x) >= 0.1 else f"{x:+.2f}%")

    diff_data = diff_data.rename(index={'allSNVs': 'SNV'})
    annot_data = annot_data.rename(index={'allSNVs': 'SNV'})
    
    # cmap_custom = LinearSegmentedColormap.from_list("PurpleYellow", ["gold", "whitesmoke", "purple"])
    cmap_custom = LinearSegmentedColormap.from_list(
        "TealBrown",
        CMAP
    )
    max_abs = np.abs(diff_data.values).max()
    max_abs_rounded = np.ceil(max_abs / 5) * 5
    # bounds = np.arange(-max_abs_rounded, max_abs_rounded + 5, 5)
    # norm = BoundaryNorm(bounds, cmap_custom.N)
    # norm = plt.Normalize(vmin=-max_abs_rounded, vmax=max_abs_rounded)

    from matplotlib.colors import TwoSlopeNorm

    norm = TwoSlopeNorm(
        vmin=-max_abs_rounded,
        vcenter=0,
        vmax=max_abs_rounded
    )

    
    # Plot without colorbar first
    sns.heatmap(diff_data, annot=annot_data, fmt='', cmap=cmap_custom, norm=norm,
               linewidths=0.5, linecolor='lightgray',
               cbar=False,
               annot_kws={'size': FONTSIZE_ANNOTATION - 1},
               ax=ax)
    
    # Create manual colorbar axis outside grid (like 4c)
    fig.canvas.draw()
    bbox_heatmap = ax.get_position()
    cbar_ax = fig.add_axes([bbox_heatmap.x1 + 0.01, bbox_heatmap.y0, 0.015, bbox_heatmap.height])
    
    # Add colorbar manually
    import matplotlib.colorbar as mcolorbar
    import matplotlib.cm as cm
    sm = cm.ScalarMappable(cmap=cmap_custom, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.outline.set_visible(False)      # removes outer box
    # cbar.ax.tick_params(size=0)          # removes tick marks

    cbar.set_label('Percentage Point Difference', fontsize=FONTSIZE_AXIS_LABEL)
    cbar.ax.tick_params(labelsize=FONTSIZE_TICK)
    
    # Square cells
    ax.set_aspect('auto', adjustable='box')
    
    ax.set_ylabel('')
    ax.set_xlabel('Calibration points', fontsize=FONTSIZE_AXIS_LABEL)
    ax.set_yticklabels(diff_data.index, rotation=90, fontsize=FONTSIZE_TICK)
    ax.set_xticklabels(score_levels, rotation=0, fontsize=FONTSIZE_TICK)
    
    # ax.text(-0.12, 1.15, "f", transform=ax.transAxes,
    #        fontsize=FONTSIZE_PANEL_LETTER, fontweight='bold', va='top', ha='left')
    
    # Add centered title in title row
    bbox_title = plt.subplot(gs_outer[0])
    bbox_title.axis('off')
    # bbox_title.text(0.5, 0.5, 'Single gene vs. gene aggregate calibration',
    #                transform=bbox_title.transAxes, ha='center', va='center',
    #                fontsize=FONTSIZE_TITLE, fontweight='bold')


def plot_figure4(
    scoreset_2018, scoreset, indv_summary, fits, score_range, n_c, n_samples, flipped,
    danzs_oob, auths_oob, dataset_names,
    prior, Post_p, Post_b, p_data_sim, b_data_sim,
    gene_4e, dist_4e, labdat_4e, snvdf_4e, sorted_thresholds_4e, oldsorted_thresholds_4e,
    dist_4f, finalout_4f,
    figsize=(13, 15)
):
    """
    Create unified Figure 4 with modular subfigures.
    Left column: a, b, c | Right column: d, e, f
    Optimized for Nature print publication.
    """
    
    fig = plt.figure(figsize=figsize)
    
    # Main grid: 2 columns, adjusted for shared legend row
    main_gs = gridspec.GridSpec(
        6, 2,
        figure=fig,
        height_ratios=[0.55, 0.08, 1.7, 0.35, 0.01, 0.8],  # Added 0.15 for legend row
        width_ratios=[1, 1],
        hspace=0.01,
        wspace=0.22,
        left=0.05, right=0.96,
        top=0.98, bottom=0.03
    )
    
    # Plot all panels
    plot_panel_a(main_gs[0, 0], scoreset_2018, indv_summary, fits, score_range, flipped, fig)
    
    all_scores = scoreset.snv_scores
    point_ranges = {int(k): v for k,v in indv_summary['point_ranges'].items()}
    legend_handles_b = plot_panel_b(main_gs[2, 0], scoreset, all_scores, point_ranges, score_range, flipped, fig)
    
    plot_panel_c(main_gs[5, 0], danzs_oob, auths_oob, fig)  # Updated row index
    
    plot_panel_d(main_gs[0, 1], prior, Post_p, Post_b, p_data_sim, b_data_sim, fig)
    
    legend_handles_e = plot_panel_e(main_gs[2, 1], gene_4e, dist_4e, labdat_4e, snvdf_4e, 
                sorted_thresholds_4e, oldsorted_thresholds_4e, fig)
    
    # Create shared legend spanning both columns
    ax_shared_legend = plt.subplot(main_gs[3, :])  # Spans both columns
    ax_shared_legend.axis('off')
    
    # Combine and deduplicate legend handles
    # Create a dict to maintain order and avoid duplicates by label
    combined_handles = {}
    for handle in legend_handles_b:
        label = handle.get_label()
        if label not in combined_handles:
            combined_handles[label] = handle
    
    for handle in legend_handles_e:
        label = handle.get_label()
        if label not in combined_handles:
            combined_handles[label] = handle
    
    # Convert back to list
    final_handles = sorted(list(combined_handles.values()), key=lambda x: int(x.get_label().split(' ')[0]))
    
    # Create centered legend
    ncol = len(final_handles)#(len(final_handles) + 1) // 2
    legend = ax_shared_legend.legend(
        handles=final_handles, 
        loc='center', 
        ncol=ncol,
        frameon=True, 
        fontsize=FONTSIZE_AXIS_LABEL, 
        columnspacing=0.6, 
        handletextpad=0.4, 
        handlelength=0.8,
        borderpad=0.7
    )

    # frame = legend.get_frame()
    # frame.set_edgecolor('black')    # Border color
    # frame.set_linewidth(1.0)        # Border thickness
    
    plot_panel_f(main_gs[5, 1], dist_4f, finalout_4f, fig)  # Updated row index

    plot_panel_letters(fig, ['a','b','c','d','e','f'])
    
    return fig






class fig_json_encoder(json.JSONEncoder):
    def default(self, obj):

        # Custom Scoreset class
        if obj.__class__.__name__ == "Scoreset":
            return obj.to_serializable()

        # NumPy arrays
        if isinstance(obj, np.ndarray):
            return {
                "__type__": "ndarray",
                "dtype": str(obj.dtype),
                "data": obj.tolist(),
            }

        # NumPy scalars
        if isinstance(obj, (np.integer, np.floating)):
            return {
                "__type__": "npscalar",
                "dtype": str(obj.dtype),
                "data": obj.item(),
            }

        # Pandas DataFrame
        if isinstance(obj, pd.DataFrame):
            return {
                "__type__": "dataframe",
                "data": obj.to_dict(orient="split"),
            }

        # Pandas Series
        if isinstance(obj, pd.Series):
            return {
                "__type__": "series",
                "data": obj.to_dict(),
                "index": obj.index.tolist(),
            }

        # Python set
        if isinstance(obj, set):
            return {
                "__type__": "set",
                "data": list(obj),
            }

        return super().default(obj)


def fig_json_hook(dct):
    if "__type__" not in dct:
        return dct

    t = dct["__type__"]
    
    if t == "Scoreset":
        return Scoreset.from_serializable(dct)

    if t == "ndarray":
        return np.array(dct["data"], dtype=dct["dtype"])

    if t == "npscalar":
        return np.dtype(dct["dtype"]).type(dct["data"])

    if t == "dataframe":
        return pd.DataFrame(**dct["data"])

    if t == "series":
        return pd.Series(dct["data"]).reindex(dct["index"])

    if t == "set":
        return set(dct["data"])

    return dct


