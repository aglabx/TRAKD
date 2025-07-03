#!/usr/bin/env python3
"""
TRAKD Repeat Statistics Analysis Tool
Analyzes and summarizes k-mer repeat patterns from BED output
"""

import sys
import argparse
import re
from collections import defaultdict, Counter
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def reverse_complement(seq):
    """Generate reverse complement of a DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def canonical_kmer(kmer):
    """Return the lexicographically smaller of kmer and its reverse complement"""
    rc = reverse_complement(kmer)
    return min(kmer, rc)

def parse_bed_file(bed_file, min_length=0):
    """Parse the detailed BED file and extract all information"""
    repeats = []
    
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('track') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            length = end - start
            
            # Skip if below minimum length
            if length < min_length:
                continue
            
            # Parse the info field - handle both semicolon-separated format
            info = parts[3]
            
            # Extract all fields
            data = {
                'chrom': chrom,
                'start': start,
                'end': end,
                'length': length
            }
            
            # Parse k-mer - it may be after locus_N; prefix
            kmer_match = re.search(r'kmer=([ACGT]+)', info)
            if kmer_match:
                kmer = kmer_match.group(1)
                data['kmer'] = kmer
                data['canonical_kmer'] = canonical_kmer(kmer)
            
            # Parse numeric fields
            for field in ['tf', 'local_tf', 'entropy', 'locus_count', 'max_locus_size', 'unique_kmers']:
                match = re.search(f'{field}=([0-9.]+)', info)
                if match:
                    data[field] = float(match.group(1))
            
            # Parse Jaccard values
            for stat in ['jaccard_min', 'jaccard_median', 'jaccard_max']:
                match = re.search(f'{stat}=([0-9.]+)', info)
                if match:
                    data[stat] = float(match.group(1))
            
            repeats.append(data)
    
    return pd.DataFrame(repeats)

def analyze_kmer_composition(df):
    """Analyze k-mer nucleotide composition"""
    kmer_stats = []
    
    for kmer in df['canonical_kmer'].unique():
        gc_content = (kmer.count('G') + kmer.count('C')) / len(kmer) * 100
        at_content = (kmer.count('A') + kmer.count('T')) / len(kmer) * 100
        
        # Check for homopolymer runs
        max_run = 1
        current_run = 1
        for i in range(1, len(kmer)):
            if kmer[i] == kmer[i-1]:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
        
        kmer_stats.append({
            'kmer': kmer,
            'gc_content': gc_content,
            'at_content': at_content,
            'max_homopolymer_run': max_run,
            'is_homopolymer': max_run == len(kmer)
        })
    
    return pd.DataFrame(kmer_stats)

def create_summary_plots(df, output_prefix):
    """Create summary visualization plots"""
    
    # Set style
    sns.set_style("whitegrid")
    
    # 1. Distribution of repeat lengths
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # Length distribution
    ax = axes[0, 0]
    df['length'].hist(bins=50, ax=ax, edgecolor='black')
    ax.set_xlabel('Repeat Length (bp)')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Repeat Lengths')
    ax.set_yscale('log')
    
    # K-mer frequency distribution
    ax = axes[0, 1]
    df['tf'].hist(bins=50, ax=ax, edgecolor='black', color='orange')
    ax.set_xlabel('K-mer Total Frequency')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of K-mer Frequencies')
    ax.set_xscale('log')
    
    # Entropy distribution
    ax = axes[1, 0]
    df['entropy'].hist(bins=50, ax=ax, edgecolor='black', color='green')
    ax.set_xlabel('Distance Entropy')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of Distance Entropy')
    
    # Local frequency vs total frequency
    ax = axes[1, 1]
    ax.scatter(df['tf'], df['local_tf'], alpha=0.5)
    ax.set_xlabel('Total Frequency')
    ax.set_ylabel('Local Frequency')
    ax.set_title('Local vs Total K-mer Frequency')
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_distributions.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. K-mer composition analysis
    kmer_comp = analyze_kmer_composition(df)
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # GC content distribution
    ax = axes[0, 0]
    kmer_comp['gc_content'].hist(bins=20, ax=ax, edgecolor='black', color='purple')
    ax.set_xlabel('GC Content (%)')
    ax.set_ylabel('Number of K-mer Types')
    ax.set_title('GC Content Distribution of K-mers')
    
    # Top k-mers by total occurrences
    ax = axes[0, 1]
    top_kmers = df.groupby('canonical_kmer')['tf'].first().sort_values(ascending=False).head(20)
    ax.barh(range(len(top_kmers)), top_kmers.values)
    ax.set_yticks(range(len(top_kmers)))
    ax.set_yticklabels(top_kmers.index, fontsize=8)
    ax.set_xlabel('Total Frequency')
    ax.set_title('Top 20 Most Frequent K-mers')
    ax.set_xscale('log')
    
    # Chromosome distribution
    ax = axes[1, 0]
    chrom_counts = df['chrom'].value_counts().sort_index()
    ax.bar(range(len(chrom_counts)), chrom_counts.values)
    ax.set_xticks(range(len(chrom_counts)))
    ax.set_xticklabels(chrom_counts.index, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Number of Repeats')
    ax.set_title('Repeat Distribution Across Chromosomes')
    
    # Entropy vs repeat count relationship
    ax = axes[1, 1]
    grouped = df.groupby('canonical_kmer').agg({
        'entropy': 'mean',
        'kmer': 'count'
    }).rename(columns={'kmer': 'repeat_count'})
    
    ax.scatter(grouped['repeat_count'], grouped['entropy'], alpha=0.5)
    ax.set_xlabel('Number of Repeat Loci')
    ax.set_ylabel('Mean Distance Entropy')
    ax.set_title('Entropy vs Repeat Count per K-mer')
    ax.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_composition.png', dpi=300, bbox_inches='tight')
    plt.close()

def generate_summary_report(df, output_file):
    """Generate a text summary report"""
    
    with open(output_file, 'w') as f:
        f.write("TRAKD Repeat Analysis Summary Report\n")
        f.write("=" * 50 + "\n\n")
        
        # Basic statistics
        f.write("Basic Statistics:\n")
        f.write(f"Total repeat loci: {len(df):,}\n")
        f.write(f"Total unique k-mers: {df['canonical_kmer'].nunique():,}\n")
        f.write(f"Total chromosomes/contigs: {df['chrom'].nunique()}\n")
        f.write(f"Total genome coverage: {df['length'].sum():,} bp\n")
        f.write(f"Mean repeat length: {df['length'].mean():.1f} bp\n")
        f.write(f"Median repeat length: {df['length'].median():.1f} bp\n")
        f.write(f"Largest repeat: {df['length'].max():,} bp\n\n")
        
        # K-mer statistics
        f.write("K-mer Statistics:\n")
        kmer_stats = df.groupby('canonical_kmer').agg({
            'tf': 'first',
            'length': ['count', 'sum', 'mean'],
            'entropy': 'mean'
        }).sort_values(('tf', 'first'), ascending=False)
        
        f.write("Top 10 most frequent k-mers:\n")
        for idx, (kmer, row) in enumerate(kmer_stats.head(10).iterrows()):
            f.write(f"  {idx+1}. {kmer}: ")
            f.write(f"TF={row[('tf', 'first')]:,}, ")
            f.write(f"Loci={row[('length', 'count')]}, ")
            f.write(f"Total_bp={row[('length', 'sum')]:,}, ")
            f.write(f"Mean_entropy={row[('entropy', 'mean')]:.2f}\n")
        
        f.write("\n")
        
        # Chromosome statistics
        f.write("Chromosome Statistics:\n")
        chrom_stats = df.groupby('chrom').agg({
            'length': ['count', 'sum', 'mean'],
            'canonical_kmer': 'nunique'
        }).sort_values(('length', 'sum'), ascending=False)
        
        for chrom, row in chrom_stats.head(10).iterrows():
            f.write(f"  {chrom}: ")
            f.write(f"Repeats={row[('length', 'count')]}, ")
            f.write(f"Total_bp={row[('length', 'sum')]:,}, ")
            f.write(f"Unique_kmers={row[('canonical_kmer', 'nunique')]}\n")
        
        f.write("\n")
        
        # Entropy analysis
        f.write("Entropy Analysis:\n")
        f.write(f"Mean entropy: {df['entropy'].mean():.2f}\n")
        f.write(f"Low entropy repeats (< 2.0): {(df['entropy'] < 2.0).sum()} ({(df['entropy'] < 2.0).sum() / len(df) * 100:.1f}%)\n")
        f.write(f"Medium entropy repeats (2.0-5.0): {((df['entropy'] >= 2.0) & (df['entropy'] < 5.0)).sum()} ({((df['entropy'] >= 2.0) & (df['entropy'] < 5.0)).sum() / len(df) * 100:.1f}%)\n")
        f.write(f"High entropy repeats (>= 5.0): {(df['entropy'] >= 5.0).sum()} ({(df['entropy'] >= 5.0).sum() / len(df) * 100:.1f}%)\n")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze TRAKD repeat patterns and generate statistics',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python analyze_repeats_stats.py detailed_loci.bed --output-prefix repeat_analysis
  python analyze_repeats_stats.py detailed_loci.bed --min-length 1000 --output-prefix filtered_analysis
        """)
    
    parser.add_argument('bed_file', help='Detailed BED file from LocusBedGeneratorDetailed')
    parser.add_argument('--output-prefix', default='repeat_stats',
                       help='Prefix for output files (default: repeat_stats)')
    parser.add_argument('--min-length', type=int, default=0,
                       help='Minimum repeat length to include (default: 0)')
    parser.add_argument('--no-plots', action='store_true',
                       help='Skip generating plots')
    
    args = parser.parse_args()
    
    # Parse BED file
    print("Parsing BED file...")
    df = parse_bed_file(args.bed_file, args.min_length)
    
    if len(df) == 0:
        print("No repeats found matching criteria")
        return
    
    print(f"Found {len(df):,} repeat loci")
    print(f"Unique k-mers: {df['canonical_kmer'].nunique():,}")
    
    # Generate summary report
    print("Generating summary report...")
    generate_summary_report(df, f'{args.output_prefix}_summary.txt')
    
    # Generate plots
    if not args.no_plots:
        print("Creating visualization plots...")
        create_summary_plots(df, args.output_prefix)
    
    # Export detailed statistics
    print("Exporting detailed statistics...")
    
    # K-mer level statistics
    kmer_stats = df.groupby('canonical_kmer').agg({
        'tf': 'first',
        'length': ['count', 'sum', 'mean', 'std', 'min', 'max'],
        'entropy': ['mean', 'std'],
        'local_tf': ['mean', 'std'],
        'locus_count': 'first',
        'max_locus_size': 'first'
    })
    kmer_stats.columns = ['_'.join(col).strip() for col in kmer_stats.columns.values]
    kmer_stats.to_csv(f'{args.output_prefix}_kmer_stats.csv')
    
    # Chromosome level statistics  
    chrom_stats = df.groupby('chrom').agg({
        'length': ['count', 'sum', 'mean', 'std'],
        'canonical_kmer': 'nunique',
        'entropy': 'mean'
    })
    chrom_stats.columns = ['_'.join(col).strip() for col in chrom_stats.columns.values]
    chrom_stats.to_csv(f'{args.output_prefix}_chrom_stats.csv')
    
    print(f"Analysis complete. Output files:")
    print(f"  - {args.output_prefix}_summary.txt")
    if not args.no_plots:
        print(f"  - {args.output_prefix}_distributions.png")
        print(f"  - {args.output_prefix}_composition.png")
    print(f"  - {args.output_prefix}_kmer_stats.csv")
    print(f"  - {args.output_prefix}_chrom_stats.csv")

if __name__ == '__main__':
    main()