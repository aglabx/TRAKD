#!/usr/bin/env python3
"""
TRAKD Circular Repeat Visualization Tool
Creates circular (Circos-style) visualization of k-mer repeats
"""

import sys
import argparse
import re
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Wedge, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.image as mpimg
import numpy as np
import colorsys
import math
import os

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

def parse_chromosome_name(chrom_name):
    """Parse chromosome name to extract number and type"""
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

def organize_chromosomes_circular(chrom_sizes):
    """Organize chromosomes for circular display"""
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

def create_circular_plot(repeats_by_chrom, chrom_sizes, kmer_colors, output_file, 
                        fig_size=14, inner_radius=0.4, outer_radius=0.9, gap_angle=2,
                        center_image=None, show_legend=False, title='Circular K-mer Repeat Map'):
    """Create circular visualization of repeats on chromosomes"""
    
    # Organize chromosomes into pairs
    pairs, sex_chroms = organize_chromosomes_circular(chrom_sizes)
    
    # Sort chromosome pairs by number
    def sort_key(x):
        num_match = re.match(r'(\d+)', x)
        if num_match:
            return (0, int(num_match.group(1)), x)
        else:
            return (1, 0, x)
    
    sorted_pairs = sorted(pairs.keys(), key=sort_key)
    
    # Calculate total size for angle calculation
    total_segments = len(sorted_pairs) + len(sex_chroms)
    angle_per_segment = (360 - gap_angle * total_segments) / total_segments
    
    # Calculate angles for each chromosome pair
    chrom_angles = {}
    current_angle = 0
    
    # Process autosome pairs
    for chr_num in sorted_pairs:
        pair = pairs[chr_num]
        
        # Both maternal and paternal share the same angular position
        if 'mat' in pair:
            chrom_angles[pair['mat']] = {
                'start': current_angle,
                'end': current_angle + angle_per_segment,
                'size': chrom_sizes[pair['mat']],
                'radius': 'outer'  # Maternal on outer ring
            }
        
        if 'pat' in pair:
            chrom_angles[pair['pat']] = {
                'start': current_angle,
                'end': current_angle + angle_per_segment,
                'size': chrom_sizes[pair['pat']],
                'radius': 'inner'  # Paternal on inner ring
            }
        
        current_angle += angle_per_segment + gap_angle
    
    # Process sex chromosomes separately
    for chr_type, chrom_name in sex_chroms.items():
        chrom_angles[chrom_name] = {
            'start': current_angle,
            'end': current_angle + angle_per_segment,
            'size': chrom_sizes[chrom_name],
            'radius': 'outer' if chr_type == 'W' else 'inner'
        }
        current_angle += angle_per_segment + gap_angle
    
    # Create figure
    fig = plt.figure(figsize=(fig_size, fig_size))
    ax = fig.add_subplot(111, projection='polar')
    
    # Define radii for maternal and paternal chromosomes
    mid_radius = (inner_radius + outer_radius) / 2
    ring_width = (outer_radius - inner_radius) / 2.5
    
    maternal_inner = mid_radius + 0.02
    maternal_outer = maternal_inner + ring_width
    paternal_inner = inner_radius
    paternal_outer = paternal_inner + ring_width
    
    # Draw chromosomes
    for chrom_name, angles in chrom_angles.items():
        start_rad = np.radians(angles['start'])
        end_rad = np.radians(angles['end'])
        
        # Determine which ring to use
        if angles['radius'] == 'outer':
            ring_inner = maternal_inner
            ring_outer = maternal_outer
        else:
            ring_inner = paternal_inner
            ring_outer = paternal_outer
        
        # Draw chromosome arc
        theta = np.linspace(start_rad, end_rad, 100)
        ax.fill_between(theta, ring_inner, ring_outer, 
                       color='lightgray', alpha=0.3, edgecolor='black', linewidth=0.5)
    
    # Add labels for chromosome pairs
    labeled_positions = set()
    for chr_num in sorted_pairs:
        pair = pairs[chr_num]
        if 'mat' in pair or 'pat' in pair:
            # Get the angular position (same for both)
            example_chrom = pair.get('mat', pair.get('pat'))
            angles = chrom_angles[example_chrom]
            mid_angle = np.radians((angles['start'] + angles['end']) / 2)
            
            if mid_angle not in labeled_positions:
                labeled_positions.add(mid_angle)
                label_radius = outer_radius + 0.05
                
                # Rotate text for better readability
                rotation = np.degrees(mid_angle) - 90
                if rotation > 90:
                    rotation -= 180
                
                ax.text(mid_angle, label_radius, chr_num,
                        ha='center', va='center', rotation=rotation, 
                        fontsize=10, weight='bold')
    
    # Add labels for sex chromosomes
    for chr_type, chrom_name in sex_chroms.items():
        angles = chrom_angles[chrom_name]
        mid_angle = np.radians((angles['start'] + angles['end']) / 2)
        label_radius = outer_radius + 0.05
        
        rotation = np.degrees(mid_angle) - 90
        if rotation > 90:
            rotation -= 180
            
        ax.text(mid_angle, label_radius, chr_type,
                ha='center', va='center', rotation=rotation, 
                fontsize=10, weight='bold')
    
    # Plot repeats
    for chrom, repeats in repeats_by_chrom.items():
        if chrom not in chrom_angles:
            continue
            
        angles = chrom_angles[chrom]
        chrom_size = angles['size']
        angle_range = angles['end'] - angles['start']
        
        # Determine which ring to use
        if angles['radius'] == 'outer':
            ring_inner = maternal_inner
            ring_outer = maternal_outer
        else:
            ring_inner = paternal_inner
            ring_outer = paternal_outer
        
        ring_height = ring_outer - ring_inner
        
        for repeat in repeats:
            # Calculate angular position
            start_angle = angles['start'] + (repeat['start'] / chrom_size) * angle_range
            end_angle = angles['start'] + (repeat['end'] / chrom_size) * angle_range
            
            # Get color
            color = kmer_colors[repeat['canonical_kmer']]
            
            # Scale height by local_tf
            height_scale = min(1.0, repeat['stats'].get('local_tf', 50) / 200)
            repeat_height = ring_height * (0.3 + 0.7 * height_scale)
            
            # Center the repeat within the ring
            repeat_inner = ring_inner + (ring_height - repeat_height) / 2
            repeat_outer = repeat_inner + repeat_height
            
            # Draw repeat
            start_rad = np.radians(start_angle)
            end_rad = np.radians(end_angle)
            
            theta = np.linspace(start_rad, end_rad, max(2, int(end_angle - start_angle)))
            ax.fill_between(theta, repeat_inner, repeat_outer,
                           color=color, alpha=0.8)
    
    # Customize plot
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)  # Clockwise
    ax.set_ylim(0, outer_radius + 0.15)
    ax.set_axis_off()
    
    # Add title
    plt.title(title, fontsize=16, weight='bold', pad=20)
    
    # Add center image if provided
    if center_image and os.path.exists(center_image):
        try:
            # Read the image
            img = mpimg.imread(center_image)
            
            # Calculate image position and size
            # Make image fit within inner radius
            image_radius = inner_radius * 0.8  # Leave some margin
            
            # Create a new axis for the image in the center
            # Convert polar coordinates to figure coordinates
            ax_image = fig.add_axes([0.5 - image_radius/2, 0.5 - image_radius/2, 
                                    image_radius, image_radius])
            ax_image.imshow(img)
            ax_image.axis('off')
            
        except Exception as e:
            print(f"Warning: Could not load center image: {e}")
    
    # Create legend only if requested
    if kmer_colors and show_legend:
        create_circular_legend(fig, sorted(kmer_colors.items()), show_all=True)
    
    # Add "Created by TRAKD" signature
    fig.text(0.99, 0.01, 'Created by TRAKD', ha='right', va='bottom', 
             fontsize=8, alpha=0.7, style='italic', transform=fig.transFigure)
    
    # Save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

def create_circular_legend(fig, kmer_color_items, show_all=True):
    """Create a compact legend for circular plot"""
    # Create a separate axis for legend - moved further right
    legend_ax = fig.add_axes([0.95, 0.1, 0.15, 0.8])
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.axis('off')
    
    # Title
    legend_ax.text(0.5, 0.98, 'K-mer Types', ha='center', va='top', 
                  fontsize=12, weight='bold')
    
    # Calculate positions
    n_items = len(kmer_color_items)
    
    if show_all:
        # Show all k-mers
        y_start = 0.94
        if n_items > 30:
            # Use smaller font and spacing for many items
            fontsize = max(4, min(8, 240 / n_items))
            y_step = (y_start - 0.05) / n_items
            box_height = y_step * 0.8
            box_width = 0.12
        else:
            fontsize = 8
            y_step = (y_start - 0.05) / n_items
            box_height = y_step * 0.6
            box_width = 0.15
        
        for i, (kmer, color) in enumerate(kmer_color_items):
            y_pos = y_start - (i * y_step)
            
            # Color box
            rect = Rectangle((0.05, y_pos - box_height/2), box_width, box_height,
                            facecolor=color, edgecolor='none')
            legend_ax.add_patch(rect)
            
            # K-mer text
            display_kmer = kmer if len(kmer) <= 10 else kmer[:8] + '..'
            legend_ax.text(0.25, y_pos, display_kmer, ha='left', va='center', 
                          fontsize=fontsize, family='monospace')
    else:
        # Original behavior - limited items
        max_items = 30
        n_items = min(len(kmer_color_items), max_items)
        
        if len(kmer_color_items) > max_items:
            legend_ax.text(0.5, 0.94, f'({max_items} of {len(kmer_color_items)})', 
                          ha='center', va='top', fontsize=9, style='italic')
            y_start = 0.90
        else:
            y_start = 0.94
        
        if n_items == 0:
            return
        y_step = (y_start - 0.05) / n_items
        
        for i, (kmer, color) in enumerate(kmer_color_items[:n_items]):
            y_pos = y_start - (i * y_step)
            
            # Color box
            rect = Rectangle((0.05, y_pos - y_step * 0.4), 0.15, y_step * 0.6,
                            facecolor=color, edgecolor='none')
            legend_ax.add_patch(rect)
            
            # K-mer text
            legend_ax.text(0.25, y_pos, kmer, ha='left', va='center', fontsize=8)

def create_heatmap_plot(repeats_by_chrom, chrom_sizes, kmer_colors, output_file,
                       bin_size=1000000, fig_width=16):
    """Create a heatmap visualization showing repeat density"""
    
    # Sort chromosomes
    sorted_chroms = sorted(chrom_sizes.keys(), key=lambda x: (-chrom_sizes[x], x))
    
    # Calculate bins for each chromosome
    chrom_bins = {}
    max_bins = 0
    
    for chrom in sorted_chroms:
        n_bins = int(np.ceil(chrom_sizes[chrom] / bin_size))
        chrom_bins[chrom] = n_bins
        max_bins = max(max_bins, n_bins)
    
    # Create density matrix
    n_chroms = len(sorted_chroms)
    density_matrix = np.zeros((n_chroms, max_bins))
    
    # Fill density matrix
    for i, chrom in enumerate(sorted_chroms):
        if chrom in repeats_by_chrom:
            for repeat in repeats_by_chrom[chrom]:
                start_bin = int(repeat['start'] / bin_size)
                end_bin = min(int(repeat['end'] / bin_size) + 1, chrom_bins[chrom])
                
                # Add repeat density weighted by local_tf
                weight = repeat['stats'].get('local_tf', 1)
                for bin_idx in range(start_bin, end_bin):
                    if bin_idx < max_bins:
                        density_matrix[i, bin_idx] += weight
    
    # Create figure
    fig_height = max(6, n_chroms * 0.5)
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Plot heatmap
    im = ax.imshow(density_matrix, aspect='auto', cmap='YlOrRd', 
                   interpolation='nearest')
    
    # Set labels
    ax.set_yticks(range(n_chroms))
    ax.set_yticklabels(sorted_chroms)
    
    # Set x-axis to show genomic positions
    x_ticks = np.arange(0, max_bins, max(1, max_bins // 10))
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f'{x * bin_size / 1e6:.0f}' for x in x_ticks])
    ax.set_xlabel('Position (Mb)', fontsize=12)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Repeat Density (weighted by local frequency)', fontsize=10)
    
    # Add title
    ax.set_title('K-mer Repeat Density Heatmap', fontsize=16, weight='bold')
    
    # Save
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Create circular visualization of TRAKD k-mer repeats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python visualize_repeats_circular.py detailed_loci.bed genome.fai circular_plot.png
  python visualize_repeats_circular.py detailed_loci.bed genome.fai circular_plot.pdf --style heatmap
        """)
    
    parser.add_argument('bed_file', help='Detailed BED file from LocusBedGeneratorDetailed')
    parser.add_argument('fai_file', help='FASTA index file (.fai)')
    parser.add_argument('output', help='Output image file (PNG, PDF, SVG)')
    parser.add_argument('--style', choices=['circular', 'heatmap'], default='circular',
                       help='Visualization style (default: circular)')
    parser.add_argument('--size', type=float, default=14, 
                       help='Figure size in inches (default: 14)')
    parser.add_argument('--min-tf', type=int, default=0,
                       help='Minimum k-mer frequency to display (default: 0)')
    parser.add_argument('--bin-size', type=int, default=1000000,
                       help='Bin size for heatmap in bp (default: 1000000)')
    parser.add_argument('--center-image', type=str, default=None,
                       help='Path to image file (PNG/JPG) to place in the center of circular plot')
    parser.add_argument('--show-legend', action='store_true',
                       help='Show k-mer legend panel (default: hidden)')
    parser.add_argument('--title', type=str, default='Circular K-mer Repeat Map',
                       help='Title for the plot (default: Circular K-mer Repeat Map)')
    
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
    print(f"Creating {args.style} visualization...")
    
    if args.style == 'circular':
        create_circular_plot(repeats_by_chrom, chrom_sizes, kmer_colors, 
                           args.output, fig_size=args.size, center_image=args.center_image,
                           show_legend=args.show_legend, title=args.title)
    else:  # heatmap
        create_heatmap_plot(repeats_by_chrom, chrom_sizes, kmer_colors,
                          args.output, bin_size=args.bin_size, fig_width=args.size)
    
    print(f"Visualization saved to: {args.output}")

if __name__ == '__main__':
    main()