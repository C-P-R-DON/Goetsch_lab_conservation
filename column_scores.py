#!/usr/bin/env python3
"""
Calculate sum of pairs column scores from a multiple sequence alignment of protein sequences.
"""
import argparse
import os
import sys
import subprocess
import tempfile
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo
from typing import List, Dict, Tuple
import numpy as np

def fetch_sequences(accessions: List[str]) -> List[SeqIO.SeqRecord]:
    """Fetch protein sequences from NCBI using accessions."""
    from Bio import Entrez
    from Bio import SeqIO
    
    print(f"Fetching {len(accessions)} sequences from NCBI...")
    
    # Set your email for NCBI (required)
    Entrez.email = "your.email@example.com"  # Please replace with your email
    
    try:
        # Fetch the sequences
        handle = Entrez.efetch(
            db="protein",
            id=",".join(accessions),
            rettype="fasta",
            retmode="text"
        )
        
        # Parse the sequences
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        
        if len(records) != len(accessions):
            print(f"Warning: Requested {len(accessions)} sequences but got {len(records)} from NCBI")
        
        return records
        
    except Exception as e:
        print(f"Error fetching sequences from NCBI: {e}", file=sys.stderr)
        sys.exit(1)

def read_fasta(fasta_file: str) -> List[SeqIO.SeqRecord]:
    """Read sequence accessions from a file and fetch sequences from NCBI."""
    try:
        # Read accessions from the file (one per line)
        with open(fasta_file, 'r') as f:
            accessions = [line.strip() for line in f if line.strip()]
        
        # Remove '>' if present (in case file is in FASTA format)
        accessions = [acc[1:] if acc.startswith('>') else acc for acc in accessions]
        
        # Fetch sequences from NCBI
        return fetch_sequences(accessions)
        
    except Exception as e:
        print(f"Error reading accessions file: {e}", file=sys.stderr)
        sys.exit(1)

def calculate_sum_of_pairs(alignment: MultipleSeqAlignment) -> List[float]:
    """
    Calculate sum of pairs score for each column in the alignment using BLOSUM62.
    
    The sum of pairs score for a column is calculated as the average of
    all pairwise BLOSUM62 substitution scores in that column.
    """
    from Bio.Align import substitution_matrices
    
    try:
        # Load the BLOSUM62 substitution matrix
        blosum62 = substitution_matrices.load("BLOSUM62")
    except Exception as e:
        print(f"Error loading BLOSUM62 matrix: {e}", file=sys.stderr)
        sys.exit(1)
    
    num_seqs = len(alignment)
    num_columns = alignment.get_alignment_length()
    
    if num_seqs < 2:
        print("Error: At least two sequences are required for alignment scoring", file=sys.stderr)
        return []
    
    def get_blosum62_score(a: str, b: str) -> float:
        """Get the BLOSUM62 score for a pair of amino acids."""
        # Convert to uppercase and handle gaps
        a = a.upper()
        b = b.upper()
        
        # Handle gaps
        if a == '-' or b == '-':
            return blosum62.get((a, b), blosum62.get((b, a), blosum62.get(('-', '-'), -4)))
        
        # Get the score from the matrix
        return blosum62.get((a, b), blosum62.get((b, a), -4))  # Default to -4 if pair not found
    
    column_scores = []
    
    for col in range(num_columns):
        column = alignment[:, col]
        col_score = 0.0
        pairs = 0
        
        # Compare all unique pairs
        for i in range(num_seqs):
            for j in range(i + 1, num_seqs):
                score = get_blosum62_score(column[i], column[j])
                col_score += score
                pairs += 1
        
        # Calculate average score per pair
        if pairs > 0:
            column_scores.append(col_score / pairs)
        else:
            column_scores.append(0.0)
    
    return column_scores

def main():
    parser = argparse.ArgumentParser(
        description='Calculate sum of pairs column scores from a multiple sequence alignment.'
    )
    parser.add_argument(
        'input_file',
        type=str,
        help='Input FASTA file containing protein sequences'
    )
    parser.add_argument(
        '-o', '--output',
        type=str,
        default='output.fasta',
        help='Output file for the multiple sequence alignment (default: output.fasta)'
    )
    
    args = parser.parse_args()
    
    # Read input sequences
    sequences = read_fasta(args.input_file)
    
    if len(sequences) < 2:
        print("Error: At least two sequences are required", file=sys.stderr)
        sys.exit(1)
    
    print(f"Read {len(sequences)} sequences from {args.input_file}")
    
    # Verify we have sequences
    if not sequences:
        print("Error: No sequences to align", file=sys.stderr)
        sys.exit(1)
        
    print(f"Retrieved {len(sequences)} sequences from NCBI")
    
    # MUSCLE path
    muscle_path = r"C:\Users\Connor\muscle-win64.v5.3.exe"
    
    if not os.path.exists(muscle_path):
        print(f"Error: MUSCLE executable not found at {muscle_path}", file=sys.stderr)
        sys.exit(1)
    
    print("Running multiple sequence alignment with MUSCLE...")
    
    # Create temporary files for input and output
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta", encoding='utf-8') as infile:
        SeqIO.write(sequences, infile, "fasta")
        infile_path = infile.name
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta", encoding='utf-8') as outfile:
        outfile_path = outfile.name
    
    try:
        # Try a different command format that works with this version of MUSCLE
        # Create a command string instead of a list to have more control over formatting
        cmd = f'"{muscle_path}" "{infile_path}" "{outfile_path}"'
        print(f"Running MUSCLE with command: {cmd}")
        
        # Run the command with shell=True since we're using a string command
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True,
                shell=True
            )
            print("MUSCLE alignment completed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"MUSCLE failed with error code {e.returncode}")
            print(f"STDOUT: {e.stdout}")
            print(f"STDERR: {e.stderr}")
            
            # Try one more time with a different format
            print("Trying alternative command format...")
            try:
                cmd = f'"{muscle_path}" -align "{infile_path}" -output "{outfile_path}"'
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True,
                    shell=True
                )
                print("MUSCLE alignment completed successfully with alternative format.")
            except subprocess.CalledProcessError as e2:
                print(f"Alternative format also failed with error: {e2.stderr}")
                raise

        # Read the aligned sequences
        alignment = MultipleSeqAlignment(SeqIO.parse(outfile_path, "fasta"))
        print("MUSCLE alignment completed successfully.")
        
    except FileNotFoundError:
        print(f"Error: MUSCLE executable not found at '{muscle_path}'.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE. It returned a non-zero exit code {e.returncode}.", file=sys.stderr)
        print(f"MUSCLE output:\n{e.stderr}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred during MUSCLE alignment: {e}", file=sys.stderr)
        sys.exit(1)
    finally:
        # Clean up the temporary files
        if infile_path and os.path.exists(infile_path):
            os.remove(infile_path)
        if outfile_path and os.path.exists(outfile_path):
            os.remove(outfile_path)
        
    # Calculate sum of pairs scores
    scores = calculate_sum_of_pairs(alignment)
    
    # Print results
    print("\nColumn scores (sum of pairs):")
    for i, score in enumerate(scores, 1):
        print(f"Column {i:4d}: {score:.4f}")
    
    # Save the alignment
    with open(args.output, 'w') as f:
        AlignIO.write(alignment, f, 'fasta')
    print(f"\nAligned sequences saved to {args.output}")

if __name__ == "__main__":
    main()