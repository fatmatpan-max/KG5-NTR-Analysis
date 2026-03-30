================================================================================
SCRIPTS/ALIGNMENT_ANALYSIS.PY
================================================================================
#!/usr/bin/env python3
"""
alignment_analysis.py

Multiple sequence alignment analysis for YdgI nitroreductase genes.
Extracts sequences, performs alignment, and calculates conservation scores.

Usage:
    python alignment_analysis.py
"""

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd
import numpy as np
from collections import Counter
import os

def calculate_conservation(alignment_file):
    """
    Calculate conservation scores for each position in alignment.

    Args:
        alignment_file: Path to FASTA alignment file

    Returns:
        list: Conservation scores (0-1) for each position
    """
    alignment = AlignIO.read(alignment_file, "fasta")
    conservation = []

    for i in range(alignment.get_alignment_length()):
        column = alignment[:, i]
        # Count non-gap characters
        bases = [b for b in column if b != '-']
        if bases:
            most_common = Counter(bases).most_common(1)[0][1]
            score = most_common / len(bases)
        else:
            score = 0
        conservation.append(score)

    return conservation

def identify_conserved_motifs(conservation_scores, threshold=0.95, min_length=10):
    """
    Identify conserved motifs based on conservation scores.

    Args:
        conservation_scores: List of conservation scores
        threshold: Minimum conservation to be considered conserved
        min_length: Minimum length of conserved region

    Returns:
        list: Tuples of (start, end) for conserved motifs
    """
    motifs = []
    in_motif = False
    start = 0

    for i, score in enumerate(conservation_scores):
        if score >= threshold and not in_motif:
            start = i
            in_motif = True
        elif score < threshold and in_motif:
            if i - start >= min_length:
                motifs.append((start, i))
            in_motif = False

    # Check if alignment ends with conserved region
    if in_motif and len(conservation_scores) - start >= min_length:
        motifs.append((start, len(conservation_scores)))

    return motifs

def main():
    # File paths
    alignment_file = "data/alignment/YdgI_alignment.fasta"
    output_file = "data/conservation_scores.csv"

    # Calculate conservation
    print("Calculating conservation scores...")
    conservation = calculate_conservation(alignment_file)

    # Identify conserved motifs
    print("Identifying conserved motifs...")
    motifs = identify_conserved_motifs(conservation)

    print(f"Found {len(motifs)} conserved motifs:")
    for start, end in motifs:
        print(f"  Position {start+1}-{end}: {end-start} bp")

    # Save to CSV
    df = pd.DataFrame({
        'Position': range(1, len(conservation) + 1),
        'Conservation_Score': conservation
    })
    df.to_csv(output_file, index=False)
    print(f"Conservation scores saved to {output_file}")

if __name__ == "__main__":
    main()


================================================================================
SCRIPTS/CONSERVATION_PLOT.PY
================================================================================
#!/usr/bin/env python3
"""
conservation_plot.py

Generate schematic alignment figure with primer binding sites.
Creates publication-ready figure for Supplementary Figure S1.

Usage:
    python conservation_plot.py
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle
import pandas as pd
import numpy as np

def get_base_color(base):
    """Return color for nucleotide base."""
    colors = {'A': '#FF6B6B', 'T': '#4ECDC4', 'G': '#FFE66D', 'C': '#95E1D3'}
    return colors.get(base.upper(), '#CCCCCC')

def plot_alignment_overview(ax, conservation_scores):
    """Plot overview of alignment with conserved motifs."""

    strains = [
        ('CP126678.1', 'B. cereus'),
        ('CP121861.1', 'B. anthracis'),
        ('CP026521.1', 'B. thuringiensis'),
        ('CP025941.1', 'B. mycoides'),
        ('CP130596.1', 'B. wiedmannii'),
        ('AZQS01000004.1', 'B. paranthracis'),
        ('NP_388447.1', 'B. subtilis (Ref)')
    ]

    y_pos = 8.5
    for acc, species in strains:
        # Background
        rect = patches.Rectangle((0, y_pos-0.3), 636, 0.6, 
                                facecolor='#E8E8E8', edgecolor='black', linewidth=1)
        ax.add_patch(rect)

        # Conservation
        for i, score in enumerate(conservation_scores):
            if i % 20 == 0:
                rect = patches.Rectangle((i, y_pos-0.3), 20, 0.6,
                                       facecolor='#C41E3A', alpha=score*0.7)
                ax.add_patch(rect)

        ax.text(-5, y_pos, acc, ha='right', va='center', fontsize=9, fontweight='bold')
        ax.text(640, y_pos, species, ha='left', va='center', fontsize=9, style='italic')
        y_pos -= 1

    # Conserved motifs
    motifs = [(1, 50, '5\' UTR'), (100, 150, 'FMN-binding'),
              (300, 350, 'Catalytic'), (400, 450, 'NAD(P)H'),
              (500, 550, 'Dimerization')]

    for start, end, label in motifs:
        rect = patches.Rectangle((start, 0.5), end-start, 0.8,
                               facecolor='yellow', edgecolor='orange',
                               linewidth=2, alpha=0.6)
        ax.add_patch(rect)
        ax.text((start+end)/2, 0.9, label, ha='center', fontsize=8, fontweight='bold')

def main():
    # Load conservation scores
    df = pd.read_csv("data/conservation_scores.csv")
    conservation_scores = df['Conservation_Score'].tolist()

    # Create figure
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1.2, 1], hspace=0.3)

    # Panel A: Overview
    ax1 = fig.add_subplot(gs[0])
    ax1.set_xlim(0, 650)
    ax1.set_ylim(0, 10)
    ax1.axis('off')
    ax1.set_title('YdgI Gene Alignment Overview', fontsize=14, fontweight='bold', pad=20)
    plot_alignment_overview(ax1, conservation_scores)

    # Primer arrows
    ax1.annotate('', xy=(80, 9.5), xytext=(15, 9.5),
                arrowprops=dict(arrowstyle='->', color='red', lw=4))
    ax1.text(47, 9.8, 'Forward Primer', ha='center', fontsize=10, color='red', fontweight='bold')

    ax1.annotate('', xy=(580, 9.5), xytext=(635, 9.5),
                arrowprops=dict(arrowstyle='->', color='blue', lw=4))
    ax1.text(607, 9.8, 'Reverse Primer', ha='center', fontsize=10, color='blue', fontweight='bold')

    # Save figure
    plt.savefig("figures/Figure_S1_alignment.png", dpi=300, bbox_inches='tight')
    print("Figure saved to figures/Figure_S1_alignment.png")

if __name__ == "__main__":
    main()
