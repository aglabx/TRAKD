#!/usr/bin/env python3
"""
TRAKD Karyotype Visualization Tool
Creates a karyotype-style visualization of k-mer repeats on chromosomes
Designed for bird genomes with maternal/paternal chromosome pairs
"""

import sys
import argparse
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Rectangle, FancyBboxPatch
import numpy as np
import colorsys

def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def canonical_kmer(kmer):
    """Return the lexicographically smaller of kmer and its reverse complement"""
    rc = reverse_complement(kmer)
    return min(kmer, rc)

def generate_distinct_colors(n, saturation_range=(0.6, 1.0), value_range=(0.5, 0.9)):
    """Generate n visually distinct colors using HSV color space"""
    colors = []
    golden_ratio = (1 + 5 ** 0.5) / 2
    
    for i in range(n):
        hue = (i * golden_ratio) % 1.0
        sat_variation = (i % 5) / 5.0
        saturation = saturation_range[0] + (saturation_range[1] - saturation_range[0]) * sat_variation
        val_variation = ((i // 5) % 3) / 3.0
        value = value_range[0] + (value_range[1] - value_range[0]) * val_variation
        rgb = colorsys.hsv_to_rgb(hue, saturation, value)
        colors.append(rgb)
    
    return colors

def parse_bed_line(line):
    """Parse a line from the detailed BED file"""
    if line.startswith('track') or not line.strip():
        return None
    
    parts = line.strip().split('\t')
    if len(parts) < 4:
        return None
    
    chrom = parts[0]
    start = int(parts[1])
    end = int(parts[2])
    
    info = parts[3]
    kmer_match = re.search(r'kmer=([ACGT]+)', info)
    if not kmer_match:
        return None
    
    kmer = kmer_match.group(1)
    
    stats = {}
    for stat in ['tf', 'local_tf', 'entropy', 'locus_count', 'max_locus_size']:
        match = re.search(f'{stat}=([0-9.]+)', info)
        if match:
            stats[stat] = float(match.group(1))
    
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'kmer': kmer,
        'canonical_kmer': canonical_kmer(kmer),
        'stats': stats
    }

def parse_fai(fai_file):
    """Parse FASTA index file to get chromosome lengths"""
    chrom_sizes = {}
    with open(fai_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes

def parse_chromosome_name(chrom_name):
    """Parse chromosome name to extract number and type (maternal/paternal)"""
    # Handle special cases
    if chrom_name == 'chrMT':
        return ('MT', None, 'special')
    elif chrom_name.startswith('chrW'):
        return ('W', 'mat', 'sex')
    elif chrom_name.startswith('chrZ'):
        return ('Z', 'pat', 'sex')
    elif chrom_name.startswith('rDNA'):
        return (chrom_name, None, 'rDNA')
    
    # Parse regular chromosomes
    match = re.match(r'chr(\d+[A-Z]?)_(mat|pat)', chrom_name)
    if match:
        return (match.group(1), match.group(2), 'autosome')
    
    return (chrom_name, None, 'unknown')

def organize_chromosomes(chrom_sizes):
    """Organize chromosomes into pairs and categories"""
    pairs = defaultdict(dict)  # chr_num -> {mat: name, pat: name}
    sex_chroms = {}
    
    for chrom_name in chrom_sizes:
        # Skip mitochondrial and rDNA
        if chrom_name in ['chrMT'] or chrom_name.startswith('rDNA'):
            continue
            
        chr_num, parent, chr_type = parse_chromosome_name(chrom_name)
        
        if chr_type == 'autosome':
            if parent:
                pairs[chr_num][parent] = chrom_name
        elif chr_type == 'sex':
            sex_chroms[chr_num] = chrom_name
    
    return pairs, sex_chroms

def draw_chromosome(ax, x, y, width, height, chrom_name, chrom_size, repeats, kmer_colors, 
                   label_above=True, signal_amplification=1.0):
    """Draw a single chromosome with repeats"""
    # Draw chromosome outline with slightly rounded ends
    chrom_patch = FancyBboxPatch((x, y), width, height,
                                boxstyle="round,pad=0.005",  # Less rounded
                                facecolor='lightgray',
                                edgecolor='black',
                                linewidth=1.0,
                                alpha=0.2)
    ax.add_patch(chrom_patch)
    
    # Add chromosome label
    label = chrom_name.replace('chr', '').replace('_mat', '♀').replace('_pat', '♂')
    label_y = y + height + 0.005 if label_above else y - 0.005  # Even closer
    va = 'bottom' if label_above else 'top'
    ax.text(x + width/2, label_y, label, ha='center', va=va, 
            fontsize=7, weight='bold')
    
    # Draw repeats
    for repeat in repeats:
        if repeat['chrom'] != chrom_name:
            continue
        
        # Calculate position on chromosome
        repeat_start = (repeat['start'] / chrom_size) * width
        repeat_end = (repeat['end'] / chrom_size) * width
        repeat_width = repeat_end - repeat_start
        
        # Apply horizontal amplification - make bands wider
        amplified_width = repeat_width * signal_amplification
        width_expansion = (amplified_width - repeat_width) / 2
        
        # Adjust position to center the expanded band
        repeat_x = x + repeat_start - width_expansion
        # Make sure we don't go outside chromosome boundaries
        repeat_x = max(x, repeat_x)
        if repeat_x + amplified_width > x + width:
            amplified_width = x + width - repeat_x
        
        # Get color
        color = kmer_colors[repeat['canonical_kmer']]
        
        # Keep height proportional but not too tall
        height_scale = min(1.0, repeat['stats'].get('local_tf', 50) / 100)
        repeat_height = height * 3.5 * (0.3 + 0.7 * height_scale)  # Increase multiplier for thin chromosomes
        repeat_y = y + (height - repeat_height) / 2
        
        # Draw repeat with enhanced width
        repeat_patch = Rectangle((repeat_x, repeat_y), amplified_width, repeat_height,
                               facecolor=color, edgecolor='none', alpha=0.9)
        ax.add_patch(repeat_patch)

def create_karyotype_plot(repeats_by_chrom, chrom_sizes, kmer_colors, output_file,
                         fig_width=16, fig_height=12, signal_amplification=1.0,
                         title='K-mer Repeat Karyotype'):
    """Create karyotype-style visualization"""
    
    # Organize chromosomes
    pairs, sex_chroms = organize_chromosomes(chrom_sizes)
    
    # Sort chromosome pairs by number
    def sort_key(x):
        # Extract numeric part for sorting
        num_match = re.match(r'(\d+)', x)
        if num_match:
            return (0, int(num_match.group(1)), x)  # Regular numbered chromosomes first
        else:
            return (1, 0, x)  # Then special chromosomes
    
    sorted_pairs = sorted(pairs.keys(), key=sort_key)
    
    # Calculate layout
    chromosomes_per_row = 10  # Balance between width and height
    num_rows = int(np.ceil((len(sorted_pairs) + 2) / chromosomes_per_row))  # +2 for sex chromosomes
    
    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Set up coordinate system
    ax.set_xlim(-0.05, 1.05)  # Extend limits to show edge chromosomes
    ax.set_ylim(-0.2, 1.1)  # Extend bottom more
    ax.set_aspect('equal')
    
    # Calculate positions - much more compact
    margin = 0.02  # Very small margins
    usable_width = 1 - 2 * margin
    usable_height = 1.1  # Use more vertical space
    
    chrom_width = usable_width / (chromosomes_per_row)
    row_height = usable_height / (num_rows + 0.5)  # Much less row height
    chrom_height = row_height * 0.02  # Even thinner chromosomes
    pair_spacing = row_height * 0.08  # Minimal space between maternal/paternal
    col_spacing = chrom_width * 1.0  # Minimal horizontal spacing
    
    # Find maximum chromosome size for scaling
    max_size = max(chrom_sizes.values())
    
    # Draw autosome pairs
    current_row = 0
    current_col = 0
    
    for chr_num in sorted_pairs:
        pair = pairs[chr_num]
        
        # Calculate position with proper spacing
        x = margin + current_col * col_spacing
        y = 1 - margin - (current_row + 1) * row_height
        
        # Draw maternal chromosome (top)
        if 'mat' in pair:
            mat_name = pair['mat']
            mat_size = chrom_sizes.get(mat_name, 0)
            scaled_width = chrom_width * (mat_size / max_size)
            draw_chromosome(ax, x, y + pair_spacing/2, scaled_width, chrom_height,
                          mat_name, mat_size, repeats_by_chrom.get(mat_name, []),
                          kmer_colors, label_above=True, signal_amplification=signal_amplification)
        
        # Draw paternal chromosome (bottom)
        if 'pat' in pair:
            pat_name = pair['pat']
            pat_size = chrom_sizes.get(pat_name, 0)
            scaled_width = chrom_width * (pat_size / max_size)
            draw_chromosome(ax, x, y - pair_spacing/2, scaled_width, chrom_height,
                          pat_name, pat_size, repeats_by_chrom.get(pat_name, []),
                          kmer_colors, label_above=False, signal_amplification=signal_amplification)
        
        # Move to next position
        current_col += 1
        if current_col >= chromosomes_per_row:
            current_col = 0
            current_row += 1
    
    # Draw sex chromosomes
    if sex_chroms:
        # Skip a column for spacing
        current_col += 1
        if current_col >= chromosomes_per_row:
            current_col = 0
            current_row += 1
        
        x = margin + current_col * col_spacing
        y = 1 - margin - (current_row + 1) * row_height
        
        # Draw Z chromosome
        if 'Z' in sex_chroms:
            z_name = sex_chroms['Z']
            z_size = chrom_sizes.get(z_name, 0)
            scaled_width = chrom_width * (z_size / max_size)
            draw_chromosome(ax, x, y + pair_spacing/2, scaled_width, chrom_height,
                          z_name, z_size, repeats_by_chrom.get(z_name, []),
                          kmer_colors, label_above=True, signal_amplification=signal_amplification)
        
        # Draw W chromosome
        if 'W' in sex_chroms:
            w_name = sex_chroms['W']
            w_size = chrom_sizes.get(w_name, 0)
            scaled_width = chrom_width * (w_size / max_size)
            draw_chromosome(ax, x, y - pair_spacing/2, scaled_width, chrom_height,
                          w_name, w_size, repeats_by_chrom.get(w_name, []),
                          kmer_colors, label_above=False, signal_amplification=signal_amplification)
    
    
    # Remove axes
    ax.axis('off')
    
    # Add title
    plt.title(title, fontsize=20, weight='bold', pad=20)
    
    # Add note about centromeres
    ax.text(0.5, -0.05, 'Note: Centromere positions will be indicated in future versions',
            ha='center', va='bottom', fontsize=10, style='italic', alpha=0.7,
            transform=ax.transAxes)
    
    # Add legend
    create_karyotype_legend(fig, ax, kmer_colors, show_all=True)  # Always show all k-mers
    
    # Add "Created by TRAKD" signature
    fig.text(0.99, 0.01, 'Created by TRAKD', ha='right', va='bottom', 
             fontsize=8, alpha=0.7, style='italic', transform=fig.transFigure)
    
    # Save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def create_karyotype_legend(fig, ax, kmer_colors, show_all=True):
    """Create a compact legend for karyotype"""
    # Create legend in the right margin - adjust position
    legend_ax = fig.add_axes([0.88, 0.1, 0.12, 0.8])
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.axis('off')
    
    # Title
    legend_ax.text(0.5, 0.98, 'K-mer Types', ha='center', va='top',
                  fontsize=12, weight='bold')
    
    # Sort k-mers for consistent display
    sorted_kmers = sorted(kmer_colors.keys())
    n_items = len(sorted_kmers)
    
    if show_all:
        # Show all k-mers, even if there are many
        y_start = 0.94
        if n_items > 30:
            # Use smaller font and spacing for many items
            fontsize = max(4, min(7, 200 / n_items))
            y_step = (y_start - 0.05) / n_items
            box_height = y_step * 0.8
        else:
            fontsize = 7
            y_step = (y_start - 0.05) / n_items
            box_height = y_step * 0.6
        
        for i, kmer in enumerate(sorted_kmers):
            y_pos = y_start - (i * y_step)
            
            # Color box
            rect = Rectangle((0.1, y_pos - box_height/2), 0.15, box_height,
                            facecolor=kmer_colors[kmer], edgecolor='none')
            legend_ax.add_patch(rect)
            
            # K-mer text (truncated if too long)
            display_kmer = kmer if len(kmer) <= 10 else kmer[:8] + '..'
            legend_ax.text(0.28, y_pos, display_kmer, ha='left', va='center',
                          fontsize=fontsize, family='monospace')
    else:
        # Original behavior - show limited items
        max_items = 20
        n_items = min(len(kmer_colors), max_items)
        if len(kmer_colors) > max_items:
            legend_ax.text(0.5, 0.94, f'({max_items} of {len(kmer_colors)})',
                          ha='center', va='top', fontsize=8, style='italic')
            y_start = 0.90
        else:
            y_start = 0.94
        
        sorted_kmers = sorted_kmers[:n_items]
        y_step = (y_start - 0.05) / n_items if n_items > 0 else 0
        
        for i, kmer in enumerate(sorted_kmers):
            y_pos = y_start - (i * y_step)
            
            # Color box
            rect = Rectangle((0.1, y_pos - y_step * 0.4), 0.2, y_step * 0.6,
                            facecolor=kmer_colors[kmer], edgecolor='none')
            legend_ax.add_patch(rect)
            
            # K-mer text
            display_kmer = kmer if len(kmer) <= 10 else kmer[:8] + '..'
            legend_ax.text(0.35, y_pos, display_kmer, ha='left', va='center',
                          fontsize=7, family='monospace')

def main():
    parser = argparse.ArgumentParser(
        description='Create karyotype-style visualization of TRAKD k-mer repeats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python visualize_karyotype.py detailed_loci.bed genome.fai karyotype.png
  python visualize_karyotype.py detailed_loci.bed genome.fai karyotype.pdf --min-tf 1000
        """)
    
    parser.add_argument('bed_file', help='Detailed BED file from LocusBedGeneratorDetailed')
    parser.add_argument('fai_file', help='FASTA index file (.fai)')
    parser.add_argument('output', help='Output image file (PNG, PDF, SVG)')
    parser.add_argument('--width', type=float, default=16,
                       help='Figure width in inches (default: 16)')
    parser.add_argument('--height', type=float, default=12,
                       help='Figure height in inches (default: 12)')
    parser.add_argument('--min-tf', type=int, default=0,
                       help='Minimum k-mer frequency to display (default: 0)')
    parser.add_argument('--amplification', type=float, default=1.0,
                       help='Horizontal signal amplification factor - makes bands wider (default: 1.0)')
    parser.add_argument('--show-all-kmers', action='store_true',
                       help='Show all k-mer types in legend (default: show up to 20)')
    parser.add_argument('--title', type=str, default='K-mer Repeat Karyotype',
                       help='Title for the plot (default: K-mer Repeat Karyotype)')
    
    args = parser.parse_args()
    
    # Parse FASTA index
    print("Reading FASTA index...")
    chrom_sizes = parse_fai(args.fai_file)
    
    # Parse BED file
    print("Reading BED file...")
    repeats_by_chrom = defaultdict(list)
    unique_kmers = set()
    
    with open(args.bed_file, 'r') as f:
        for line in f:
            repeat = parse_bed_line(line)
            if repeat and repeat['stats'].get('tf', 0) >= args.min_tf:
                repeats_by_chrom[repeat['chrom']].append(repeat)
                unique_kmers.add(repeat['canonical_kmer'])
    
    print(f"Found {len(unique_kmers)} unique k-mer types (considering reverse complements)")
    print(f"Total repeats: {sum(len(r) for r in repeats_by_chrom.values())}")
    
    if len(unique_kmers) == 0:
        print("Error: No k-mer repeats found in the BED file.")
        print("Please check the BED file format and --min-tf parameter")
        return
    
    # Generate colors for k-mers
    print("Generating color palette...")
    colors = generate_distinct_colors(len(unique_kmers))
    kmer_colors = dict(zip(sorted(unique_kmers), colors))
    
    # Create visualization
    print("Creating karyotype visualization...")
    create_karyotype_plot(repeats_by_chrom, chrom_sizes, kmer_colors,
                         args.output, fig_width=args.width, fig_height=args.height,
                         signal_amplification=args.amplification, title=args.title)
    
    print(f"Karyotype visualization saved to: {args.output}")

if __name__ == '__main__':
    main()