import numpy as np 
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams["pdf.fonttype"] = 42  
matplotlib.rcParams["ps.fonttype"] = 42 

from pathlib import Path
import numpy as np
import keras
import crested

import pandas as pd
import re
from typing import Optional, Tuple
from pysam import FastaFile

from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap, ScalarMappable

def reorder_predictions(pred_matrix: np.ndarray, output_names: list[str], desired_order: list[str]) -> tuple[np.ndarray, list[str]]:
    """
    Reorder columns of prediction matrix to match a desired class order,
    keeping only classes that exist in the output names.

    Parameters
    ----------
    pred_matrix : np.ndarray
        2D array of shape (N, C) with prediction scores.
    output_names : list[str]
        Current class names corresponding to columns of pred_matrix.
    desired_order : list[str]
        Desired order of class names.

    Returns
    -------
    reordered_matrix : np.ndarray
        Prediction matrix with columns reordered and filtered.
    reordered_class_names : list[str]
        Class names in the new order corresponding to reordered_matrix columns.
    """
    name_to_index = {name: i for i, name in enumerate(output_names)}
    reordered_indices = [name_to_index[name] for name in desired_order if name in name_to_index]
    reordered_names = [name for name in desired_order if name in name_to_index]
    reordered_matrix = pred_matrix[:, reordered_indices]
    return reordered_matrix, reordered_names

def load_macaque_chrom_map(filepath):
    chrom_map = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                old_name, new_name = parts[0], parts[2][:-1]
                chrom_map[old_name] = new_name
    return chrom_map

import pandas as pd
import re
from typing import Optional, Tuple
from pysam import FastaFile

def _raw_and_padded(
    coord_str: str,
    genome: FastaFile,
    target_length: int,
    chrom_map: Optional[dict],
    plasmid_id: str
) -> Tuple[Optional[str], Optional[str], Optional[int], Optional[int]]:
    if pd.isna(coord_str):
        print(f"[DEBUG][{plasmid_id}] No coordinate → skipping")
        return None, None, None, None

    coord_clean = coord_str.replace(",", "")
    m = re.match(r"(chr[\w\d_]+):(\d+)-(\d+)", coord_clean)
    if not m:
        print(f"[WARN][{plasmid_id}] Bad coord format '{coord_str}'")
        return None, None, None, None
    chrom, s_s, e_s = m.groups()
    start, end = int(s_s), int(e_s)
    print(f"[DEBUG][{plasmid_id}] Original coords = {chrom}:{start}-{end}")

    if chrom_map and chrom not in chrom_map:
        print(f"[WARN][{plasmid_id}] Chrom {chrom} not in mapping")
        return None, None, None, None
    if chrom_map:
        chrom = chrom_map[chrom]

    try:
        raw_seq = genome.fetch(chrom, start, end).upper()
    except Exception as exc:
        print(f"[ERROR][{plasmid_id}] Fetch failed {chrom}:{start}-{end} → {exc}")
        return None, None, None, None

    raw_len = len(raw_seq)
    total_pad = target_length - raw_len
    if total_pad < 0:
        raise ValueError(f"[{plasmid_id}] target {target_length} < raw length {raw_len}")

    left = total_pad // 2
    right = total_pad - left
    padded = "N"*left + raw_seq + "N"*right

    print(f"[DEBUG][{plasmid_id}] Raw length = {raw_len}")
    print(f"[DEBUG][{plasmid_id}] Original padded = {len(padded)} (L={left}, R={right})")

    return raw_seq, padded, start, end

def process_row(
    row: pd.Series,
    species_to_genome: dict,
    target_length: int,
    macaque_chrom_map: Optional[dict]
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    plasmid_id = row.get("Plasmid ID", row.name)
    coord      = row["Coordinates"]
    species    = row["Species"]
    core       = row.get("Core", None)

    if pd.isna(coord) or pd.isna(species) or "/" in coord or "/" in species:
        print(f"[DEBUG][{plasmid_id}] Missing coord/species → skipping")
        return None, None, None

    genome = species_to_genome.get(species)
    if genome is None:
        print(f"[WARN][{plasmid_id}] No genome for species '{species}'")
        return None, None, None

    chrom_map = macaque_chrom_map if species == "Macaque" else None

    # 1) Fetch raw and original padded
    raw_seq, orig_padded, gstart, gend = _raw_and_padded(
        coord, genome, target_length, chrom_map, plasmid_id
    )
    if raw_seq is None:
        return None, None, None

    # If no core or explicit “weirdo”, skip bashing
    if pd.isna(core) or str(core).lower() == "weirdo":
        print(f"[DEBUG][{plasmid_id}] core='{core}' → no bashing")
        return orig_padded, coord, orig_padded

    core_str = str(core)

    # 2a) “Nx” with no segment → repeat entire raw_seq
    m_full = re.match(r"(\d+)x$", core_str)
    if m_full:
        repeats = int(m_full.group(1))
        print(f"[DEBUG][{plasmid_id}] core='{core_str}' → repeating full sequence {repeats}×")
        bashed = raw_seq * repeats
        # record core_coords = original coords
        core_coords = coord
    else:
        # 2b) “NxM” segment logic
        m = re.match(r"(\d+)x(\d+)", core_str)
        if not m:
            print(f"[WARN][{plasmid_id}] Bad core format '{core_str}' → skipping bashing")
            return orig_padded, coord, orig_padded

        repeats, idx = map(int, m.groups())
        seg_i = idx - 1
        if not (0 <= seg_i < 3):
            print(f"[WARN][{plasmid_id}] Core index {idx} out of 1..3 → skipping")
            return orig_padded, coord, orig_padded

        raw_len = len(raw_seq)
        inner   = raw_len - 50
        third   = inner / 3.0
        start_i = int(seg_i * third)
        end_i   = int(start_i + third + 50)
        end_i   = min(end_i, raw_len)

        print(f"[DEBUG][{plasmid_id}] core='{core_str}': raw[{start_i}:{end_i}) "
              f"len={end_i-start_i}, repeat={repeats}×")

        segment = raw_seq[start_i:end_i]
        bashed  = segment * repeats

        # compute core_coords
        core_start = gstart + start_i
        core_end   = gstart + end_i
        core_coords = f"{coord.split(':')[0]}:{core_start}-{core_end}"
        print(f"[DEBUG][{plasmid_id}] Core coords = {core_coords}")

    # 3) Pad bashed
    bashed_len = len(bashed)
    pad_total  = target_length - bashed_len
    if pad_total < 0:
        raise ValueError(f"[{plasmid_id}] after bashing len {bashed_len} > target {target_length}")
    lpad = pad_total // 2
    rpad = pad_total - lpad
    bashed_padded = "N"*lpad + bashed + "N"*rpad

    print(f"[DEBUG][{plasmid_id}] After bashing unpadded len = {bashed_len}")
    print(f"[DEBUG][{plasmid_id}] Final padded len = {len(bashed_padded)} (L={lpad}, R={rpad})")

    return bashed_padded, core_coords, orig_padded

def assign_sequences(
    df: pd.DataFrame,
    species_to_genome: dict,
    target_length: int = 2114,
    macaque_chrom_map: Optional[dict] = None
) -> pd.DataFrame:
    results = df.apply(
        lambda row: process_row(row, species_to_genome, target_length, macaque_chrom_map),
        axis=1,
        result_type="expand"
    )
    results.columns = ["Sequence", "Core_Coordinates", "Original_Padded"]
    return pd.concat([df, results], axis=1)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import get_cmap, ScalarMappable

def plot_prediction_scatter_grid(prediction_strength,
                                  prediction_specificity,
                                  df_valid,
                                  output_names,
                                  size_scale=300,
                                  font_size=9,
                                  figsize=(14, 17),
                                 cmap='magma',
                                save_path=None,
    ):
    """
    Creates a labeled scatter dot plot with two colorbars:
    - Dot color encodes prediction strength
    - Dot size encodes prediction specificity
    - Side strip shows level of labeling per row

    Parameters:
    - prediction_strength: 2D array-like of shape (n_rows, n_cols)
    - prediction_specificity: 2D array-like of shape (n_rows, n_cols)
    - df_valid: DataFrame with index as row labels and a 'Level of labeling' column
    - output_names: List of column labels (length n_cols)
    - size_scale: Scalar to control dot size
    - font_size: Font size for axis and colorbar labels
    - figsize: Size of the full figure
    """
    n_rows, n_cols = prediction_strength.shape
    x, y = np.meshgrid(np.arange(n_cols), np.arange(n_rows))
    x = x.flatten()
    y = y.flatten()

    # Dot size and color
    sizes = prediction_specificity.flatten() * size_scale
    color_norm = Normalize(vmin=prediction_strength.min(), vmax=prediction_strength.max())
    colors = get_cmap(cmap)(color_norm(prediction_strength.flatten()))

    # Level of labeling
    level_series = df_valid["Level of labeling"]
    level_norm = Normalize(vmin=0, vmax=level_series.max())
    level_cmap = get_cmap("plasma")
    level_colors = level_cmap(level_norm(level_series.values))

    # Begin plotting
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(1, 4, width_ratios=[0.15, 0.02, 1, 0.05], wspace=0.05)

    # Row labels
    ax_labels = fig.add_subplot(gs[0, 0])
    ax_labels.set_xlim(0, 1)
    ax_labels.set_ylim(n_rows - 0.5, -0.5)
    for i, label in enumerate(df_valid.index):
        ax_labels.text(1, i, label, ha='right', va='center', fontsize=font_size)
    ax_labels.axis("off")

    # Color strip (flipped to match scatter)
    #level_colors_flipped = level_colors[::-1]
    #ax_colors = fig.add_subplot(gs[0, 1])
    #ax_colors.imshow(level_colors_flipped[:, np.newaxis], aspect="auto", origin="upper")
    #ax_colors.set_xticks([])
    #ax_colors.set_yticks(np.arange(n_rows))
    #ax_colors.set_yticklabels([])
    #ax_colors.invert_yaxis()

    from matplotlib.patches import Rectangle

    ax_colors = fig.add_subplot(gs[0, 1])
    ax_colors.set_xlim(0, 1)
    ax_colors.set_ylim(n_rows, 0)
    
    for i, color in enumerate(level_colors):
        rect = Rectangle((0, i), 1, 1, color=color)
        ax_colors.add_patch(rect)
    
    ax_colors.axis("off")

    # Main scatter plot
    ax_main = fig.add_subplot(gs[0, 2])
    sc = ax_main.scatter(x, y, s=sizes, c=colors)
    ax_main.set_xticks(np.arange(n_cols))
    ax_main.set_xticklabels(output_names, rotation=90)
    ax_main.set_yticks(np.arange(n_rows))
    ax_main.set_yticklabels([])  # already shown in ax_labels
    ax_main.set_ylim(n_rows - 0.5, -0.5)
    ax_main.set_xlabel("Group")
    ax_main.set_ylabel("")
    ax_main.tick_params(axis='x', labelsize=font_size)
    ax_labels.tick_params(axis='y', labelsize=font_size)

    # Colorbars
    cax_strength = fig.add_axes([0.91, 0.55, 0.015, 0.35])
    cb = plt.colorbar(ScalarMappable(norm=color_norm, cmap=cmap), cax=cax_strength)
    cb.set_label("Prediction strength (log)", fontsize=font_size)
    cb.ax.tick_params(labelsize=font_size)

    cax_level = fig.add_axes([0.91, 0.10, 0.015, 0.35])
    cb2 = plt.colorbar(ScalarMappable(norm=level_norm, cmap=level_cmap), cax=cax_level)
    cb2.set_label("Level of labeling", fontsize=font_size)
    cb2.ax.tick_params(labelsize=font_size)
    cb2.ax.yaxis.set_label_position('right')
    cb2.ax.yaxis.set_ticks_position('right')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()