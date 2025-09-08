#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
align_and_tree_for_codeml.py

A comprehensive pipeline for evolutionary analysis of protein-coding sequences. This script:
1. Fetches protein and corresponding nucleotide sequences from NCBI
2. Performs multiple sequence alignment using MUSCLE
3. Creates codon-aware nucleotide alignments using pal2nal.pl
4. Builds phylogenetic trees
5. Performs selection analysis using PAML's codeml
6. Visualizes results with branch-specific dN/dS ratios

Dependencies:
- BioPython (biopython)
- MUSCLE (for multiple sequence alignment)
- pal2nal.pl (for codon alignment)
- PAML (for dN/dS calculation)
- Matplotlib (for visualization)

Author: Cascade
Date: 2025-07-11
Version: 1.0
"""

import argparse
import os
import sys
import subprocess
import time
import tempfile
from io import StringIO
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo.PAML import codeml
import matplotlib.pyplot as plt


# --- Constants ---
# Delay between NCBI Entrez requests (in seconds) to comply with usage guidelines
# Note: NCBI recommends no more than 3 requests per second without an API key
NCBI_REQUEST_DELAY = 0.5


def fetch_sequences(protein_acc, nucleotide_acc):
    """
    Fetches and validates protein-CDS sequence pairs from NCBI Entrez.
    
    This function retrieves both protein and coding DNA sequences (CDS) from NCBI,
    handles cases where multiple records are returned, and ensures the sequences
    are properly paired and valid for downstream analysis.

    Args:
        protein_acc (str): The protein accession number (e.g., 'XP_12345.1').
        nucleotide_acc (str): The nucleotide accession number (e.g., 'NM_123456').

    Returns:
        tuple: A tuple containing (protein_seqrecord, nucleotide_seqrecord) if successful,
              or (None, None) if any error occurs or no matching pair is found.
              
    Note:
        - Validates that nucleotide sequence length is a multiple of 3
        - Verifies that the translated CDS matches the protein sequence
        - Cleans sequence IDs for consistent use in downstream analysis
    """
    print(f"Fetching sequences for {protein_acc} and {nucleotide_acc}...")
    try:
        # Fetch protein sequence(s)
        time.sleep(NCBI_REQUEST_DELAY)
        handle_prot = Entrez.efetch(db="protein", id=protein_acc, rettype="fasta", retmode="text")
        protein_records = list(SeqIO.parse(handle_prot, "fasta"))
        handle_prot.close()
        if not protein_records:
            print(f"Error: No protein record found for {protein_acc}", file=sys.stderr)
            return None, None

        # Fetch nucleotide sequence(s)
        time.sleep(NCBI_REQUEST_DELAY)
        handle_nucl = Entrez.efetch(db="nuccore", id=nucleotide_acc, rettype="fasta_cds_na", retmode="text")
        nucleotide_records = list(SeqIO.parse(handle_nucl, "fasta"))
        handle_nucl.close()
        if not nucleotide_records:
            print(f"Error: No nucleotide CDS found for {nucleotide_acc}", file=sys.stderr)
            return None, None

        # --- Find the matching pair ---
        protein_record = None
        nucleotide_record = None
        if len(protein_records) == 1 and len(nucleotide_records) == 1:
            # The simple case, one of each.
            protein_record = protein_records[0]
            nucleotide_record = nucleotide_records[0]
        else:
            print(f"Info: Found {len(protein_records)} protein and {len(nucleotide_records)} nucleotide records. Searching for a match...")
            found_match = False
            for prot_rec in protein_records:
                for nuc_rec in nucleotide_records:
                    if len(nuc_rec.seq) % 3 != 0:
                        continue
                    
                    translated_cds = nuc_rec.seq.translate(to_stop=True)
                    if str(prot_rec.seq).strip('*') == str(translated_cds).strip('*'):
                        protein_record = prot_rec
                        nucleotide_record = nuc_rec
                        print(f"Info: Found matching protein/CDS pair for {protein_acc}/{nucleotide_acc}.")
                        found_match = True
                        break
                if found_match:
                    break
            
            if not found_match:
                print(f"Warning: Could not find a matching protein/CDS pair for {protein_acc}/{nucleotide_acc}.", file=sys.stderr)
                return None, None

        # --- Validation and ID cleaning (on the matched pair) ---
        if len(nucleotide_record.seq) % 3 != 0:
            print(f"Error: Nucleotide sequence length for {nucleotide_acc} is not a multiple of 3.", file=sys.stderr)
            return None, None

        translated_cds = nucleotide_record.seq.translate(to_stop=True)
        if str(protein_record.seq).strip('*') != str(translated_cds).strip('*'):
            print(f"Warning: Translated CDS for {nucleotide_acc} does not match protein {protein_acc}. pal2nal will handle this.", file=sys.stderr)

        # Use a clean ID for alignment that matches for both records
        clean_id = protein_acc.split('.')[0]
        protein_record.id = clean_id
        nucleotide_record.id = clean_id
        protein_record.description = ""
        nucleotide_record.description = ""

        return protein_record, nucleotide_record

    except Exception as e:
        print(f"Error fetching sequences for {protein_acc}/{nucleotide_acc}: {e}", file=sys.stderr)
        return None, None


def run_muscle_alignment(protein_records, muscle_path):
    """
    Performs multiple sequence alignment of protein sequences using MUSCLE.
    
    This function handles the alignment process using temporary files to ensure
    robustness, especially with large sequence sets. It cleans up all temporary
    files after alignment is complete.
    
    Args:
        protein_records (list): List of Bio.SeqRecord objects to be aligned.
        muscle_path (str): Path to the MUSCLE executable.
        
    Returns:
        Bio.Align.MultipleSeqAlignment: The aligned sequences, or None if alignment fails.
        
    Raises:
        FileNotFoundError: If MUSCLE executable is not found at the specified path.
        subprocess.CalledProcessError: If MUSCLE encounters an error during execution.
    """
    print("Running MUSCLE for protein alignment using temporary files...")
    infile_path, outfile_path = None, None
    try:
        # Create temporary files for input and output
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta", encoding='utf-8') as infile:
            SeqIO.write(protein_records, infile, "fasta")
            infile_path = infile.name
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta", encoding='utf-8') as outfile:
            outfile_path = outfile.name

        # Use the standard '-align' and '-output' flags for MUSCLE v5.
        command = [muscle_path, "-align", infile_path, "-output", outfile_path]

        subprocess.run(command, capture_output=True, text=True, check=True)

        # The output of a MUSCLE alignment is in aligned FASTA format.
        alignment = MultipleSeqAlignment(SeqIO.parse(outfile_path, "fasta"))
        print("MUSCLE alignment completed successfully.")
        return alignment

    except FileNotFoundError:
        print(f"Error: MUSCLE executable not found at '{muscle_path}'.", file=sys.stderr)
        return None
    except subprocess.CalledProcessError as e:
        print(f"Error running MUSCLE. It returned a non-zero exit code {e.returncode}.", file=sys.stderr)
        print(f"MUSCLE output:\n{e.stdout or e.stderr}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred during MUSCLE alignment: {e}", file=sys.stderr)
        return None
    finally:
        # Clean up the temporary files
        if infile_path and os.path.exists(infile_path):
            os.remove(infile_path)
        if outfile_path and os.path.exists(outfile_path):
            os.remove(outfile_path)


def run_pal2nal(protein_alignment, nucleotide_records, pal2nal_path):
    """
    Creates codon-aware nucleotide alignments from protein alignments using pal2nal.pl.
    
    This function bridges the gap between protein and nucleotide sequence spaces
    by using the protein alignment as a guide to align the corresponding
    nucleotide sequences in a codon-aware manner.
    
    Args:
        protein_alignment (MultipleSeqAlignment): Aligned protein sequences from MUSCLE.
        nucleotide_records (dict): Dictionary mapping sequence IDs to unaligned
                                 nucleotide SeqRecords (must match protein records).
        pal2nal_path (str): Path to the pal2nal.pl Perl script.
        
    Returns:
        Bio.Align.MultipleSeqAlignment: Codon-aligned nucleotide sequences in the
                                      same order as the input protein alignment.
                                      Returns None if alignment fails.
                                       
    Note:
        - Requires Perl to be installed and in the system PATH
        - The 'nogap' option is used to remove columns with only gaps
        - Handles temporary file creation and cleanup automatically
    """
    print("Running pal2nal.pl to generate codon alignment...")
    prot_aln_file, nuc_fasta_file = None, None
    try:
        # Write protein alignment to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta", encoding='utf-8') as temp_prot_aln:
            SeqIO.write(protein_alignment, temp_prot_aln, "fasta")
            prot_aln_file = temp_prot_aln.name

        # Write unaligned nucleotide sequences to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".fasta", encoding='utf-8') as temp_nuc_fasta:
            SeqIO.write(nucleotide_records.values(), temp_nuc_fasta, "fasta")
            nuc_fasta_file = temp_nuc_fasta.name

        # Construct the pal2nal.pl command
        # NOTE: The -frameshift option is not supported in pal2nal v14 and has been removed.
        # For frameshift correction, consider upgrading to the latest version of pal2nal.pl.
        command = [
            "perl", pal2nal_path,
            prot_aln_file,
            nuc_fasta_file,
            "-output", "fasta",  # Output in FASTA format for easy parsing
            "-nogap"  # Remove columns with only gaps
        ]
        
        # Run pal2nal.pl and capture its output
        result = subprocess.run(command, capture_output=True, text=True, check=True)

        # Parse the FASTA output from stdout
        codon_alignment_records = list(SeqIO.parse(StringIO(result.stdout), "fasta"))
        
        if not codon_alignment_records:
            print("Error: pal2nal.pl produced an empty output.", file=sys.stderr)
            print(f"Stderr:\n{result.stderr}", file=sys.stderr)
            return None

        print("pal2nal.pl completed successfully.")
        return MultipleSeqAlignment(codon_alignment_records)

    except FileNotFoundError:
        print(f"Error: 'perl' or '{pal2nal_path}' not found.", file=sys.stderr)
        print("Please ensure Perl is installed and in your system's PATH,", file=sys.stderr)
        print(f"and that the path to pal2nal.pl is correct.", file=sys.stderr)
        return None
    except subprocess.CalledProcessError as e:
        print(f"Error running pal2nal.pl. Return code: {e.returncode}", file=sys.stderr)
        print(f"Stderr:\n{e.stderr}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred during pal2nal.pl execution: {e}", file=sys.stderr)
        return None
    finally:
        # Clean up temporary files
        if prot_aln_file and os.path.exists(prot_aln_file):
            os.remove(prot_aln_file)
        if nuc_fasta_file and os.path.exists(nuc_fasta_file):
            os.remove(nuc_fasta_file)


def write_phylip_for_paml(alignment, filename):
    """
    Writes a sequence alignment in PHYLIP format compatible with PAML.
    
    This function handles the specific formatting requirements of PAML, including:
    - Strict 10-character limit on sequence IDs
    - Two spaces between ID and sequence
    - Proper handling of sequence wrapping
    
    Args:
        alignment (MultipleSeqAlignment): The alignment to be written.
        filename (str): Path to the output file.
        
    Note:
        - Truncates sequence IDs to 10 characters if necessary
        - Pads shorter IDs with spaces to maintain fixed-width format
        - Wraps sequences to 50 characters per line (standard PHYLIP format)
    """
    print(f"Writing alignment to {filename} in PAML-compatible PHYLIP format...")
    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    with open(filename, "w") as f_out:
        # Write the header
        f_out.write(f" {num_sequences} {alignment_length}\n")
        
        # Write each sequence
        for record in alignment:
            # PAML requires a strict 10-character name.
            # Truncate or pad with spaces to make it exactly 10 characters.
            seq_id = record.id[:10].ljust(10) 
            
            # The sequence data
            seq_data = str(record.seq)
            
            # Write the formatted line: 10-char name, TWO spaces, sequence
            f_out.write(f"{seq_id}  {seq_data}\n")
    print("Custom PHYLIP file written successfully.")


def main():
    """
    Main execution function for the alignment and phylogenetic analysis pipeline.
    
    This function orchestrates the entire workflow:
    1. Parses command-line arguments
    2. Fetches sequences from NCBI
    3. Performs protein sequence alignment
    4. Creates codon-aware nucleotide alignments
    5. Builds phylogenetic trees
    6. Performs selection analysis using PAML
    7. Generates visualizations
    
    Command-line Arguments:
        --input: TSV file with protein and nucleotide accessions
        --output: Base name for output files
        --email: NCBI email (required)
        --entrez-api-key: NCBI API key (optional)
        --muscle-path: Path to MUSCLE executable
        --pal2nal-path: Path to pal2nal.pl script
        --codeml-path: Path to PAML codeml executable
    
    Output Files:
        {output}.phy: PHYLIP format alignment
        {output}.nwk: Newick format tree
        {output}_omega_tree.png: Visualization of tree with dN/dS values
        {output}.mlc: PAML codeml output
    """
    parser = argparse.ArgumentParser(
        description="Aligns proteins and uses pal2nal.pl to create codon alignments in PHYLIP format."
    )
    parser.add_argument("--input", required=True, help="Input TSV file with protein and nucleotide accessions.")
    parser.add_argument("--output", required=True, help="Output file for the codon alignment in PHYLIP format (e.g., 'codon_alignment.phy').")
    parser.add_argument("--email", required=True, help="Your email for NCBI Entrez.")
    parser.add_argument("--entrez-api-key", help="Optional NCBI API key for Entrez.")
    parser.add_argument("--muscle-path", default="muscle", help="Path to the MUSCLE executable.")
    parser.add_argument("--pal2nal-path", default="pal2nal.pl", help="Path to the pal2nal.pl script.")
    parser.add_argument("--codeml-path", default="codeml.exe", help="Path to the codeml executable.")

    args = parser.parse_args()

    # Setup Entrez
    Entrez.email = args.email
    if args.entrez_api_key:
        Entrez.api_key = args.entrez_api_key

    # --- Step 1 & 2: Parse input and fetch sequences ---
    protein_records = []
    nucleotide_records = {}
    try:
        with open(args.input, 'r') as f_in:
            # next(f_in) # Uncomment if your TSV has a header
            for line in f_in:
                parts = line.strip().split('\t')
                if len(parts) < 4:
                    continue
                
                protein_acc = parts[2]
                nucleotide_acc = parts[3]

                prot_rec, nuc_rec = fetch_sequences(protein_acc, nucleotide_acc)
                if prot_rec and nuc_rec:
                    protein_records.append(prot_rec)
                    nucleotide_records[prot_rec.id] = nuc_rec
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input}", file=sys.stderr)
        sys.exit(1)

    if not protein_records:
        print("Error: No valid sequences could be fetched. Exiting.", file=sys.stderr)
        sys.exit(1)

    # --- Step 3: Align proteins with MUSCLE ---
    protein_alignment = run_muscle_alignment(protein_records, args.muscle_path)
    if not protein_alignment:
        sys.exit(1)

    # --- Step 4: Translate protein alignment to codon alignment using pal2nal.pl ---
    codon_alignment = run_pal2nal(protein_alignment, nucleotide_records, args.pal2nal_path)
    if not codon_alignment:
        print("Failed to generate codon alignment. Exiting.", file=sys.stderr)
        sys.exit(1)

    # --- Step 5: Write output in PAML-compatible PHYLIP format ---
    write_phylip_for_paml(codon_alignment, args.output + ".phy")

    print("Processing complete.")

    # Load the PHYLIP-formatted alignment
    with open(args.output + ".phy", "r") as f_in:
        alignment = AlignIO.read(f_in, "phylip")

    # Compute the distance matrix and build the tree
    calculator = DistanceCalculator("identity")
    dm = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    # Strip branch lengths and bootstrap/internal node labels
    for clade in tree.find_clades():
        clade.branch_length = None
        clade.confidence = None  # Clears bootstrap/confidence values

    # Save the cleaned tree in Newick format
    Phylo.write([tree], args.output + ".nwk", "newick")

    # run codeml to annotate the tree
    codeml_executable = os.path.abspath(args.codeml_path)
    print(f"Attempting to use codeml executable at: {codeml_executable}")

    if not os.path.isfile(codeml_executable):
        print(f"Error: The specified codeml executable was not found or is not a file: '{codeml_executable}'.", file=sys.stderr)
        print("Please check the path provided via the --codeml-path argument.", file=sys.stderr)
        sys.exit(1)
    
    try:
        cml = codeml.Codeml()
        cml.codeml = codeml_executable # Set attribute for consistency
        cml.read_ctl_file("codeml.ctl")
    #    try:
        print("Running codeml...")
        # Explicitly pass the command to the run method for robustness
        results = cml.run(command=codeml_executable, verbose=True, parse=True)
        print("codeml analysis completed.")
        
        # Read the MLC file to get omega values
        def parse_mlc_for_omegas(mlc_file):
            omega_values = {}
            with open(mlc_file, 'r') as f:
                in_branch_section = False
                for line in f:
                    line = line.strip()
                    if 'TREE' in line and 'w' in line and 'dN/dS' in line:
                        in_branch_section = True
                        continue
                    if in_branch_section and line.startswith('('):
                        # This is the tree line with omega values
                        # Example format: (XP_015733360.1:0.0163[&&NHX:w=0.00000],XP_015739798.1:0.0155[&&NHX:w=0.00000])
                        import re
                        # Extract all omega values and their corresponding branch names
                        pattern = r'([A-Za-z0-9_.]+):[^\[\]]*\[&&NHX:w=([0-9.]+)\]'
                        matches = re.findall(pattern, line)
                        for branch_name, omega in matches:
                            omega_values[branch_name] = float(omega)
                        break
            return omega_values
        
        # Get omega values from MLC file
        mlc_file = f"{args.output}.mlc"
        omega_values = parse_mlc_for_omegas(mlc_file)
        
        if not omega_values:
            print("Warning: No omega values found in MLC file. Trying alternative parsing method.")
            # Alternative method using the results object
            if 'branches' in results:
                for branch, params in results['branches'].items():
                    if 'w' in params:
                        omega_values[branch] = params['w']
        
        if not omega_values and 'parameters' in results and 'omega' in results['parameters']:
            # If still no values, use global omega
            global_omega = results['parameters']['omega']
            print(f"Using global omega value: {global_omega}")
            # We'll handle this after reading the tree
        
        # Read the tree from the MLC file
        tree = None
        with open(mlc_file, 'r') as f:
            for line in f:
                if line.startswith('('):
                    # This is the tree line
                    import io
                    tree = Phylo.read(io.StringIO(line), 'newick')
                    break
        
        if tree is None:
            raise ValueError("Could not parse tree from MLC file")
        
        # If we only have a global omega, apply it to all branches
        if not omega_values and 'parameters' in results and 'omega' in results['parameters']:
            for clade in tree.find_clades():
                if clade.name:
                    omega_values[clade.name] = results['parameters']['omega']
        
        # Prepare branch labels with omega values
        for clade in tree.find_clades():
            if clade.name in omega_values:
                clade.branch_length = omega_values[clade.name]
        
        # Create and save the tree visualization
        import matplotlib.pyplot as plt
        plt.figure(figsize=(16, 10))
        ax = plt.gca()
        
        # Draw the tree with omega values as branch labels
        Phylo.draw(
            tree,
            axes=ax,
            do_show=False,
            branch_labels=lambda c: f"{c.branch_length:.3f}" if c.branch_length is not None else "",
            label_func=lambda x: x.name if x.name else ""
        )
        
        # Add title and save
        plt.title("Phylogenetic Tree with dN/dS (ω) Values")
        plt.tight_layout()
        output_png = f"{args.output}_omega_tree.png"
        plt.savefig(output_png, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Tree visualization saved to {output_png}")
        
    except (FileNotFoundError, OSError) as e:  # Catch more general OS errors
        print(f"Error executing codeml: {e}", file=sys.stderr)
        print(f"Failed to run codeml from path: '{codeml_executable}'.", file=sys.stderr)
        print("Please ensure the PAML suite is installed correctly, the path is correct,", file=sys.stderr)
        print("and the file has execute permissions.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing results: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
