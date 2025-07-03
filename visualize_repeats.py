#!/usr/bin/env python3
"""
TRAKD Repeat Visualization Tool
Visualizes k-mer repeats on chromosomes using signature k-mers for coloring
"""

import sys
import argparse
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
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
    
    # Golden ratio for better distribution
    golden_ratio = (1 + 5 ** 0.5) / 2
    
    for i in range(n):
        # Use golden ratio for hue distribution
        hue = (i * golden_ratio) % 1.0
        
        # Vary saturation and value for more distinct colors
        sat_variation = (i % 5) / 5.0
        saturation = saturation_range[0] + (saturation_range[1] - saturation_range[0]) * sat_variation
        
        val_variation = ((i // 5) % 3) / 3.0
        value = value_range[0] + (value_range[1] - value_range[0]) * val_variation
        
        # Convert HSV to RGB
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
    
    # Parse the info field - handle both semicolon-separated format
    info = parts[3]
    
    # Extract k-mer - it may be after locus_N; prefix
    kmer_match = re.search(r'kmer=([ACGT]+)', info)
    if not kmer_match:
        return None
    
    kmer = kmer_match.group(1)
    
    # Extract other statistics
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

def create_chromosome_plot(repeats_by_chrom, chrom_sizes, kmer_colors, output_file, 
                          fig_width=16, chrom_height=0.8, chrom_spacing=1.2,
                          title='K-mer Repeat Distribution Across Chromosomes'):
    """Create visualization of repeats on chromosomes"""
    
    # Sort chromosomes by size (largest first) or by name
    sorted_chroms = sorted(chrom_sizes.keys(), key=lambda x: (-chrom_sizes[x], x))
    
    # Calculate figure height
    n_chroms = len(sorted_chroms)
    fig_height = (n_chroms * chrom_spacing) + 2
    
    # Create figure
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Plot each chromosome
    y_positions = {}
    for i, chrom in enumerate(sorted_chroms):
        y_pos = i * chrom_spacing
        y_positions[chrom] = y_pos
        
        # Draw chromosome backbone
        chrom_length = chrom_sizes[chrom]
        chrom_rect = patches.Rectangle((0, y_pos - chrom_height/2), 
                                      chrom_length, chrom_height,
                                      linewidth=1, edgecolor='black', 
                                      facecolor='lightgray', alpha=0.3)
        ax.add_patch(chrom_rect)
        
        # Add chromosome label
        ax.text(-chrom_length * 0.02, y_pos, chrom, 
                ha='right', va='center', fontsize=10, weight='bold')
        
        # Add size label
        size_mb = chrom_length / 1_000_000
        ax.text(chrom_length + chrom_length * 0.01, y_pos, f'{size_mb:.1f} Mb', 
                ha='left', va='center', fontsize=8, alpha=0.7)
    
    # Plot repeats
    for chrom, repeats in repeats_by_chrom.items():
        if chrom not in y_positions:
            continue
            
        y_pos = y_positions[chrom]
        
        for repeat in repeats:
            color = kmer_colors[repeat['canonical_kmer']]
            
            # Scale height by local_tf or entropy
            height_scale = min(1.0, repeat['stats'].get('local_tf', 50) / 200)
            repeat_height = chrom_height * (0.5 + 0.5 * height_scale)
            
            repeat_rect = patches.Rectangle(
                (repeat['start'], y_pos - repeat_height/2),
                repeat['end'] - repeat['start'], repeat_height,
                linewidth=0, facecolor=color, alpha=0.8
            )
            ax.add_patch(repeat_rect)
    
    # Set axis properties
    ax.set_xlim(-max(chrom_sizes.values()) * 0.05, max(chrom_sizes.values()) * 1.05)
    ax.set_ylim(-1, n_chroms * chrom_spacing)
    ax.set_xlabel('Position (bp)', fontsize=12)
    ax.set_title(title, fontsize=16, weight='bold')
    
    # Remove y-axis
    ax.yaxis.set_visible(False)
    
    # Add grid
    ax.grid(True, axis='x', alpha=0.3)
    
    # Format x-axis
    ax.ticklabel_format(style='scientific', axis='x', scilimits=(6, 6))
    
    # Create legend
    create_legend(fig, ax, kmer_colors)
    
    # Add "Created by TRAKD" signature
    fig.text(0.99, 0.01, 'Created by TRAKD', ha='right', va='bottom', 
             fontsize=8, alpha=0.7, style='italic', transform=fig.transFigure)
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def create_legend(fig, ax, kmer_colors, max_items=50):
    """Create a legend showing k-mer to color mapping"""
    # Sort k-mers by frequency (if available) or alphabetically
    sorted_kmers = sorted(kmer_colors.keys())
    
    if len(sorted_kmers) > max_items:
        # Show only the most common ones
        legend_text = f"Showing {max_items} of {len(sorted_kmers)} unique k-mer types"
        ax.text(0.02, 0.98, legend_text, transform=ax.transAxes, 
                fontsize=10, va='top', style='italic')
        sorted_kmers = sorted_kmers[:max_items]
    
    # Create legend patches
    legend_patches = []
    for kmer in sorted_kmers:
        color = kmer_colors[kmer]
        patch = patches.Patch(color=color, label=kmer)
        legend_patches.append(patch)
    
    # Add legend to the right side
    legend = ax.legend(handles=legend_patches, loc='center left', 
                      bbox_to_anchor=(1.02, 0.5), ncol=2,
                      fontsize=8, title='K-mer Signatures')
    
    # Make legend scrollable if too many items
    legend.set_title('K-mer Signatures', prop={'size': 10, 'weight': 'bold'})

def main():
    parser = argparse.ArgumentParser(
        description='Visualize TRAKD k-mer repeats on chromosomes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python visualize_repeats.py detailed_loci.bed genome.fai repeats_plot.png
  python visualize_repeats.py detailed_loci.bed genome.fai repeats_plot.pdf --width 20 --min-tf 1000
        """)
    
    parser.add_argument('bed_file', help='Detailed BED file from LocusBedGeneratorDetailed')
    parser.add_argument('fai_file', help='FASTA index file (.fai)')
    parser.add_argument('output', help='Output image file (PNG, PDF, SVG)')
    parser.add_argument('--width', type=float, default=16, 
                       help='Figure width in inches (default: 16)')
    parser.add_argument('--min-tf', type=int, default=0,
                       help='Minimum k-mer frequency to display (default: 0)')
    parser.add_argument('--chromosomes', nargs='+',
                       help='Specific chromosomes to plot (default: all)')
    parser.add_argument('--title', type=str, default='K-mer Repeat Distribution Across Chromosomes',
                       help='Title for the plot')
    
    args = parser.parse_args()
    
    # Parse FASTA index
    print("Reading FASTA index...")
    chrom_sizes = parse_fai(args.fai_file)
    
    # Filter chromosomes if specified
    if args.chromosomes:
        chrom_sizes = {c: s for c, s in chrom_sizes.items() if c in args.chromosomes}
    
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
        print("Please check:")
        print("1. The BED file path is correct") 
        print("2. The BED file format matches LocusBedGeneratorDetailed output")
        print("3. Try lowering --min-tf parameter if filtering by frequency")
        return
    
    # Generate colors for k-mers
    print("Generating color palette...")
    colors = generate_distinct_colors(len(unique_kmers))
    kmer_colors = dict(zip(sorted(unique_kmers), colors))
    
    # Create visualization
    print("Creating visualization...")
    create_chromosome_plot(repeats_by_chrom, chrom_sizes, kmer_colors, 
                          args.output, fig_width=args.width, title=args.title)
    
    print(f"Visualization saved to: {args.output}")

if __name__ == '__main__':
    main()