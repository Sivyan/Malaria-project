"""
Preliminary Analysis of Falcipain Protein Sequences
This script performs comprehensive MSA analysis including:
1. Conservation score calculation
2. Key residue analysis
3. Visualization of conservation profiles
4. Comparison of conservation across different versions
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from Bio import AlignIO


def calculate_conservation_scores(alignment):
    """Calculate conservation scores for an MSA."""
    sequences = [list(str(record.seq)) for record in alignment]
    df = pd.DataFrame(sequences)
    valid_amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    df = df.applymap(lambda x: x if x in valid_amino_acids else "-")
    return df.apply(lambda col: col.value_counts(normalize=True).max(), axis=0)


def analyze_msa_file(file_path, output_dir):
    """Analyze a single MSA file and create visualization."""
    alignment = AlignIO.read(file_path, "fasta")
    conservation_scores = calculate_conservation_scores(alignment)
    plt.figure(figsize=(12, 6))
    plt.plot(range(1, len(conservation_scores) + 1), conservation_scores, 
             marker='o', linestyle='-', color='blue', alpha=0.6, markersize=2)
    
    residues_to_highlight = [285, 417, 447]
    for pos in residues_to_highlight:
        if pos <= len(conservation_scores):
            plt.axvline(x=pos, color='r', linestyle='--', alpha=0.3)
            plt.plot(pos, conservation_scores[pos-1], 'ro', markersize=8)
    
    plt.title(f'Conservation Scores - {os.path.basename(file_path)}')
    plt.xlabel('Position')
    plt.ylabel('Conservation Score')
    plt.ylim(0.4, 1)
    plt.grid(True, alpha=0.3)
    
    os.makedirs(output_dir, exist_ok=True)
    
    output_plot = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(file_path))[0]}_conservation.png")
    plt.savefig(output_plot)
    plt.close()
    
    highly_conserved = [(i+1, score) for i, score in enumerate(conservation_scores) if score > 0.95]
    
    return highly_conserved, conservation_scores


def analyze_key_residues(alignment, positions=[285, 417, 447]):
    """Analyze conservation at specific residue positions."""
    conservation_scores = calculate_conservation_scores(alignment)
    results = []
    
    for pos in positions:
        if pos <= len(conservation_scores):
            results.append({
                'Position': pos,
                'Conservation_Score': conservation_scores[pos-1]
            })
    
    return pd.DataFrame(results)


def process_folder(folder_path):
    """Process all MSA files in a given folder."""
    version_name = os.path.basename(folder_path)
    output_dir = os.path.join("results", version_name)
    os.makedirs(output_dir, exist_ok=True)
    
    msa_files = glob.glob(os.path.join(folder_path, "*.afa"))
    print(f"\nProcessing {version_name}: Found {len(msa_files)} MSA files")
    
    all_results = []
    for msa_file in sorted(msa_files):
        print(f"\nAnalyzing {os.path.basename(msa_file)}...")
        
        try:
            highly_conserved, conservation_scores = analyze_msa_file(msa_file, output_dir)
            
            alignment = AlignIO.read(msa_file, "fasta")
            key_residues = analyze_key_residues(alignment)
            
            all_results.append({
                'Version': version_name,
                'MSA_file': os.path.basename(msa_file),
                'Key_Residues': key_residues,
                'Highly_Conserved': highly_conserved
            })
            
            output_file = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(msa_file))[0]}_conserved.txt")
            with open(output_file, 'w') as f:
                f.write(f"Highly conserved residues (score > 0.95) in {os.path.basename(msa_file)}:\n")
                f.write("Position\tConservation Score\n")
                for pos, score in highly_conserved:
                    f.write(f"{pos}\t{score:.3f}\n")
            
            print(f"Analysis complete. Results saved to {output_file}")
            
        except Exception as e:
            print(f"Error analyzing {msa_file}: {str(e)}")
    
    return all_results


def main():
    os.makedirs("results", exist_ok=True)
    
    all_version_results = []
    for version in ['v1', 'v2', 'v3']:
        if os.path.exists(version):
            print(f"\nProcessing {version}...")
            results = process_folder(version)
            all_version_results.extend(results)
        else:
            print(f"\nWarning: {version} folder not found")
    
    if all_version_results:
        results_df = pd.DataFrame(all_version_results)
        
        output_file = os.path.join("results", "combined_analysis.tsv")
        results_df.to_csv(output_file, sep='\t', index=False)
        print(f"\nCombined results saved to {output_file}")


if __name__ == "__main__":
    main() 