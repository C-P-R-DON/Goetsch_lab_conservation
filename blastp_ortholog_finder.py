#!/usr/bin/env python3

import argparse
import sys
import time
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW, NCBIXML
import re
import os
import subprocess

# --- Configuration ---\n# IMPORTANT: Always tell NCBI who you are
Entrez.email = "cdonahue@mtu.edu"  # Replace with your email address
Entrez.api_key = None  # Replace with your NCBI API key if you have one, otherwise leave as None or remove

# NCBI Request settings
NCBI_REQUEST_DELAY = 1  # seconds between Entrez queries
BLAST_RETRY_LIMIT = 3
BLAST_RETRY_DELAY = 5  # seconds
E_VALUE_THRESHOLD = 0.05

# --- Helper Functions ---

def validate_nucleotide_accession(accession):
    """Validate if the accession is a valid nucleotide accession using NCBI Entrez."""
    try:
        time.sleep(NCBI_REQUEST_DELAY)
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = handle.read()
        handle.close()
        if not record or not record.startswith(">"):
            return False
        return True
    except Exception as e:
        print(f"Error validating accession {accession}: {e}", file=sys.stderr)
        return False

def get_nucleotide_sequence_and_organism(accession):
    """Retrieve nucleotide sequence and organism from NCBI."""
    try:
        time.sleep(NCBI_REQUEST_DELAY)
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        if not records:
            return None, None
        
        seq_record = records[0]
        sequence = seq_record.get("GBSeq_sequence", "").upper()
        organism = seq_record.get("GBSeq_organism", "N/A")
        return sequence, organism
    except Exception as e:
        print(f"Error fetching sequence/organism for {accession}: {e}", file=sys.stderr)
        return None, None

def get_taxonomy_details(species_name):
    """Fetch phylum for a given species name."""
    if not species_name or species_name == "N/A":
        return "N/A"
    try:
        time.sleep(NCBI_REQUEST_DELAY)
        # Search for the TaxID
        handle = Entrez.esearch(db="taxonomy", term=species_name)
        record = Entrez.read(handle)
        handle.close()
        if not record["IdList"]:
            print(f"Warning: Could not find TaxID for species '{species_name}'.", file=sys.stderr)
            return "N/A"
        tax_id = record["IdList"][0]

        time.sleep(NCBI_REQUEST_DELAY)
        # Fetch the taxonomy record for the TaxID
        handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        if records and records[0].get("LineageEx"):
            for taxon in records[0]["LineageEx"]:
                if taxon.get("Rank") == "phylum":
                    return taxon.get("ScientificName", "N/A")
        print(f"Warning: Phylum not found for species '{species_name}' (TaxID: {tax_id}).", file=sys.stderr)
        return "N/A"
    except Exception as e:
        print(f"Error fetching taxonomy for {species_name}: {e}", file=sys.stderr)
        return "N/A"


def get_linked_protein_from_nucleotide(nucleotide_accession):
    """
    Finds a protein sequence linked to the given nucleotide accession.
    Resolves GI numbers to full RefSeq accessions (e.g., XP_, NP_).
    """
    print(f"Attempting to find linked protein for nucleotide accession: {nucleotide_accession}...")
    link_names_to_try = ["nuccore_protein_cds", "nuccore_protein"]
    protein_ids = []

    for link_name in link_names_to_try:
        try:
            time.sleep(NCBI_REQUEST_DELAY)
            handle = Entrez.elink(dbfrom="nuccore", db="protein", id=nucleotide_accession, linkname=link_name)
            results = Entrez.read(handle)
            handle.close()
            if results and results[0].get("LinkSetDb") and results[0]["LinkSetDb"][0].get("Link"):
                protein_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
                if protein_ids:
                    break
        except Exception as e:
            print(f"Error during elink for {nucleotide_accession} with linkname {link_name}: {e}", file=sys.stderr)

    if not protein_ids:
        print(f"No linked protein IDs found for {nucleotide_accession}.", file=sys.stderr)
        return None, None

    try:
        time.sleep(NCBI_REQUEST_DELAY)
        handle = Entrez.esummary(db="protein", id=",".join(protein_ids))
        summaries = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Error fetching protein summaries for IDs {protein_ids}: {e}", file=sys.stderr)
        return None, None

    # Extract all valid accession numbers from the summaries
    accessions = [s.get("AccessionVersion") for s in summaries if s.get("AccessionVersion")]
    
    protein_accession = None
    # Prioritize RefSeq accessions
    for prefix in ["XP_", "NP_", "WP_"]:
        for acc in accessions:
            if acc and acc.startswith(prefix):
                protein_accession = acc
                break
        if protein_accession:
            break
    
    # Fallback to the first accession if no preferred one is found
    if not protein_accession and accessions:
        protein_accession = accessions[0]
        print(f"Warning: No preferred RefSeq protein accession found. Falling back to '{protein_accession}'.", file=sys.stderr)

    if not protein_accession:
        print(f"Could not resolve a valid protein accession from IDs: {protein_ids}", file=sys.stderr)
        return None, None

    print(f"Found linked protein accession: {protein_accession}. Fetching its sequence...")
    try:
        time.sleep(NCBI_REQUEST_DELAY)
        handle_fasta = Entrez.efetch(db="protein", id=protein_accession, rettype="fasta", retmode="text")
        fasta_record_str = handle_fasta.read()
        handle_fasta.close()
        if not fasta_record_str or not fasta_record_str.startswith(">"):
            return protein_accession, None
        protein_sequence = "".join(fasta_record_str.splitlines()[1:]).strip()
        return protein_accession, protein_sequence
    except Exception as e:
        print(f"Error fetching sequence for protein {protein_accession}: {e}", file=sys.stderr)
        return protein_accession, None

def run_blast_search(input_sequence, evalue_threshold, blast_type, database, target_phylum_taxid=None, blast_path=None):
    """Submit BLAST search and retrieve results, either via NCBIWWW or local command line."""
    if database and blast_path:
        # Local BLAST execution
        print(f"Running local {blast_type} against database: {database}")
        temp_input_file = "temp_blast_input.fasta"
        with open(temp_input_file, "w") as f:
            f.write(f">query\n{input_sequence}")

        temp_output_file = "temp_blast_output.xml"
        blast_executable = os.path.join(blast_path, blast_type)

        cmd = [
            blast_executable,
            "-query", temp_input_file,
            "-db", database,
            "-out", temp_output_file,
            "-outfmt", "5",  # XML format
            "-evalue", str(evalue_threshold),
            "-max_target_seqs", "5000"
        ]
        if target_phylum_taxid:
            cmd.extend(["-taxids", str(target_phylum_taxid)])

        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            with open(temp_output_file, "r") as result_handle:
                blast_records = list(NCBIXML.parse(result_handle))
            os.remove(temp_input_file)
            os.remove(temp_output_file)
            return blast_records
        except subprocess.CalledProcessError as e:
            print(f"Local BLAST search failed: {e}", file=sys.stderr)
            print(f"Stderr: {e.stderr}", file=sys.stderr)
            os.remove(temp_input_file)
            if os.path.exists(temp_output_file):
                os.remove(temp_output_file)
            return None
        except FileNotFoundError:
            print(f"Error: BLAST executable not found at '{blast_executable}'. Please check the --blast_path.", file=sys.stderr)
            return None

    # NCBIWWW execution (fallback)
    for attempt in range(BLAST_RETRY_LIMIT):
        try:
            print("Submitting BLAST search...")
            result_handle = NCBIWWW.qblast(
                program=blast_type,
                database=database,
                sequence=input_sequence,
                expect=evalue_threshold,
                hitlist_size=5000, # Get more hits to ensure species diversity and better filtering
                entrez_query=f"txid{target_phylum_taxid}[Organism]" if target_phylum_taxid else None
            )
            if target_phylum_taxid:
                print(f"BLAST search is limited to Phylum TaxID: {target_phylum_taxid}")
            
            print("BLAST search submitted. Waiting for results...")
            # NCBIXML.parse returns an iterator of Blast records.
            # For typical qblast, it's usually one record.
            blast_records = list(NCBIXML.parse(result_handle)) 
            result_handle.close()
            return blast_records
        except Exception as e:
            print(f"BLAST search attempt {attempt + 1} failed: {e}", file=sys.stderr)
            if attempt < BLAST_RETRY_LIMIT - 1:
                print(f"Retrying in {BLAST_RETRY_DELAY} seconds...", file=sys.stderr)
                time.sleep(BLAST_RETRY_DELAY)
            else:
                print("BLAST search failed after multiple retries.", file=sys.stderr)
                return None
    return None

def get_nucleotide_accession_for_protein(protein_accession):
    """
    Attempt to find a linked nucleotide accession for a protein accession.
    Resolves GI numbers to full RefSeq accessions (e.g., NM_, XM_).
    """
    if not protein_accession:
        return "N/A"
    
    link_names_to_try = ["protein_nuccore_mrna", "protein_nuccore_cds", "protein_nuccore"]
    nucleotide_ids = []
    for link_name in link_names_to_try:
        try:
            time.sleep(NCBI_REQUEST_DELAY)
            handle = Entrez.elink(dbfrom="protein", db="nuccore", id=protein_accession, linkname=link_name)
            results = Entrez.read(handle)
            handle.close()

            if results and results[0].get("LinkSetDb") and results[0]["LinkSetDb"][0].get("Link"):
                nucleotide_ids = [link["Id"] for link in results[0]["LinkSetDb"][0]["Link"]]
                if nucleotide_ids:
                    break
        except Exception as e:
            print(f"Warning: elink failed for {protein_accession} with linkname {link_name}: {e}", file=sys.stderr)

    if not nucleotide_ids:
        print(f"Warning: Could not retrieve nucleotide IDs for protein {protein_accession}.", file=sys.stderr)
        return "N/A"

    try:
        time.sleep(NCBI_REQUEST_DELAY)
        summary_handle = Entrez.esummary(db="nuccore", id=",".join(nucleotide_ids))
        summaries = Entrez.read(summary_handle)
        summary_handle.close()
    except Exception as e:
        print(f"Error fetching nucleotide summaries for IDs {nucleotide_ids}: {e}", file=sys.stderr)
        return "N/A"

    accessions = [s.get("AccessionVersion") for s in summaries if s.get("AccessionVersion")]
    
    nucleotide_accession = None
    # Prioritize RefSeq mRNA accessions
    for prefix in ["NM_", "XM_"]:
        for acc in accessions:
            if acc and acc.startswith(prefix):
                nucleotide_accession = acc
                break
        if nucleotide_accession:
            break
    
    if not nucleotide_accession and accessions:
        nucleotide_accession = accessions[0]
        print(f"Warning: No preferred RefSeq nucleotide accession found. Falling back to '{nucleotide_accession}'.", file=sys.stderr)

    if nucleotide_accession:
        return nucleotide_accession
    else:
        print(f"Warning: Could not resolve a valid nucleotide accession for protein {protein_accession}.", file=sys.stderr)
        return "N/A"

# --- Main Script Logic ---

def main():
    global E_VALUE_THRESHOLD
    parser = argparse.ArgumentParser(
        description="Performs NCBI BLAST search using a nucleotide or protein accession depending on type of BLAST search run, "
                    "translates it (if running blastp), and retrieves significant hits per species."
    )
    parser.add_argument("seed_accession", help="NCBI accession of the seed nucleotide sequence (e.g., NM_001301339)")
    parser.add_argument("seed_gene_name", help="Common name of the seed gene (e.g., BRCA1)")
    parser.add_argument("blast_type", help="Type of BLAST search to run ('blastn', 'tblastn', 'blastp', 'blastx', 'tblastx')")
    parser.add_argument("--email", help="Email address for NCBI Entrez queries. Overrides script default.")
    parser.add_argument("--apikey", help="NCBI API key. Overrides script default.")
    parser.add_argument("--evalue", type=float, default=E_VALUE_THRESHOLD, help=f"E-value threshold for BLAST (default: {E_VALUE_THRESHOLD})")
    parser.add_argument("--phylum_taxid", type=int, default=None, help="Optional: NCBI Taxonomy ID of the target phylum to restrict BLAST search (e.g., 6656 for Arthropoda, 7711 for Chordata).")
    parser.add_argument("--blast_db", help="Path to the local BLAST database. If provided, command-line BLAST will be used.")
    parser.add_argument("--blast_path", help="Path to the directory containing BLAST+ executables.")

    args = parser.parse_args()

    E_VALUE_THRESHOLD = args.evalue # Update global E_VALUE_THRESHOLD from argument

    if args.email:
        Entrez.email = args.email
    if args.apikey:
        Entrez.api_key = args.apikey
    
    if not Entrez.email or Entrez.email == "YOUR_EMAIL_HERE":
        print("Error: NCBI Entrez email not set. Please edit the script or use --email argument.", file=sys.stderr)
        sys.exit(1)

    print(f"Using NCBI Email: {Entrez.email}")
    if Entrez.api_key and Entrez.api_key != "YOUR_API_KEY_HERE" and Entrez.api_key is not None: # Check if it's not None
        print(f"Using NCBI API Key: {'*' * (len(Entrez.api_key) - 4) + Entrez.api_key[-4:]}")
    else:
        print("NCBI API Key not set or default. Using default NCBI limits (3 requests/second without key, 10 with key).")

    # 1. Input Validation
    print(f"Validating seed accession: {args.seed_accession}...")
    if not validate_nucleotide_accession(args.seed_accession):
        print(f"Error: Invalid or non-nucleotide seed accession '{args.seed_accession}'.", file=sys.stderr)
        sys.exit(1)
    
    if not args.seed_gene_name.strip():
        print("Error: Seed gene name cannot be empty.", file=sys.stderr)
        sys.exit(1)
    print("Input validation successful.")

    # 2. Organism Retrieval and Linked Protein Sequence Fetching
    print(f"Retrieving organism information for {args.seed_accession}...")
    nucleotide_sequence, seed_organism_raw = get_nucleotide_sequence_and_organism(args.seed_accession)
    if not nucleotide_sequence or not seed_organism_raw:
        print(f"Error: Could not retrieve sequence or organism for {args.seed_accession}. Exiting.", file=sys.stderr)
        sys.exit(1)
    
    print(f"Fetching linked protein sequence for {args.seed_accession}...")
    linked_protein_accession, protein_sequence = get_linked_protein_from_nucleotide(args.seed_accession)
    if not protein_sequence:
        print(f"Error: Failed to retrieve a linked protein sequence for {args.seed_accession}. Exiting.", file=sys.stderr)
        sys.exit(1)
    print(f"Using protein {linked_protein_accession} (Length: {len(protein_sequence)} amino acids) for BLAST search.")

    # 3. Seed Species Identification
    print(f"Identifying seed species: {seed_organism_raw}...")
    seed_phylum = get_taxonomy_details(seed_organism_raw)
    seed_species_name = seed_organism_raw 

    # 4. BLAST Search
    # Determine input sequence and database for BLAST
    if args.blast_type in ['blastp', 'tblastn']:
        input_sequence_for_blast = protein_sequence
    elif args.blast_type in ['blastn', 'blastx', 'tblastx']:
        input_sequence_for_blast = nucleotide_sequence
    else:
        print(f"Error: Unsupported BLAST type '{args.blast_type}' provided.", file=sys.stderr)
        sys.exit(1)

    # Determine database
    if args.blast_db:
        database = args.blast_db
    elif args.blast_type in ['blastp', 'blastx', 'tblastn']:
        database = 'nr'
    elif args.blast_type in ['blastn', 'tblastx']:
        database = 'nt'
    else:
        print(f"Error: Could not determine database for BLAST type '{args.blast_type}'.", file=sys.stderr)
        sys.exit(1)
    blast_records_list = run_blast_search(input_sequence_for_blast, args.evalue, args.blast_type, database, args.phylum_taxid, args.blast_path)
    if not blast_records_list:
        print("Error: BLAST search did not return results or failed. Exiting.", file=sys.stderr)
        sys.exit(1)
    
    # 5. Result Parsing
    print("Parsing BLAST results...")
    output_data_list = []
    # species_name -> (evalue, bit_score, hit_data_dict)
    best_hits_per_species = {} 

    for blast_record in blast_records_list: # Should be only one record from qblast typically
        if not blast_record.alignments:
            print("Warning: BLAST record contains no alignments.", file=sys.stderr)
            continue
        for alignment in blast_record.alignments:
            hit_title = alignment.title
            hit_species_from_title = "N/A"
            # A more robust, multi-step approach to parse species names.
            # 1. Try to find a species name in brackets, e.g., [Homo sapiens]
            bracket_match = re.search(r'\[([A-Za-z]+\s[a-z]+(?:[\sA-Za-z0-9.-]+)?)\]', hit_title)
            if bracket_match:
                hit_species_from_title = bracket_match.group(1).strip()
            else:
                # 2. If not, try to find "Genus species" or "Genus sp." after "PREDICTED: "
                # This is common in RefSeq mRNA titles.
                species_pattern = r'([A-Z][a-z]+[ \t]+(?:[a-z]+|sp\.))'
                predicted_match = re.search(r'PREDICTED: ' + species_pattern, hit_title)
                if predicted_match:
                    hit_species_from_title = predicted_match.group(1).strip()
                else:
                    # 3. If not, try to find "Genus species" after a pipe and space, common in other DBs
                    pipe_match = re.search(r'\|\s*' + species_pattern, hit_title)
                    if pipe_match:
                        hit_species_from_title = pipe_match.group(1).strip()
                    else:
                        # 4. As a last resort, take the first two words after "PREDICTED: "
                        # This is a bit of a guess but can catch formats missed above.
                        last_resort_match = re.search(r'PREDICTED: ([A-Za-z]+\s[A-Za-z\d.-]+)', hit_title)
                        if last_resort_match:
                            hit_species_from_title = last_resort_match.group(1).strip()
                        else:
                            print(f"Warning: Could not parse species from title: '{hit_title}'.", file=sys.stderr)

            if not alignment.hsps:
                continue
            hsp = alignment.hsps[0]
            e_value = hsp.expect
            bit_score = float(hsp.score)

            # Filter: Keep only the most significant hit per species
            current_key = hit_species_from_title
            if current_key in best_hits_per_species:
                prev_evalue, prev_bit_score, _ = best_hits_per_species[current_key]
                # If the new hit is worse (higher E-value, or same E-value and lower/equal bit score), skip it.
                if e_value > prev_evalue or (e_value == prev_evalue and bit_score <= prev_bit_score):
                    continue
            
            # This is the best hit for this species so far. Process based on the type of hits returned by the BLAST search.
            if args.blast_type in ['blastn', 'tblastn', 'tblastx']:
                # These BLAST types return nucleotide hits.
                hit_nucleotide_accession = alignment.accession
                print(f"Processing best hit for species: '{current_key}' (Nucleotide Acc: {hit_nucleotide_accession}, E-value: {e_value:.2e}, Bit-score: {bit_score})")
                hit_protein_accession, _ = get_linked_protein_from_nucleotide(hit_nucleotide_accession)
                if not hit_protein_accession:
                    hit_protein_accession = "N/A"
            elif args.blast_type in ['blastp', 'blastx']:
                # These BLAST types return protein hits.
                hit_protein_accession = alignment.accession
                print(f"Processing best hit for species: '{current_key}' (Protein Acc: {hit_protein_accession}, E-value: {e_value:.2e}, Bit-score: {bit_score})")
                hit_nucleotide_accession = get_nucleotide_accession_for_protein(hit_protein_accession)
            else:
                print(f"Warning: Unhandled BLAST type '{args.blast_type}' in result parsing.", file=sys.stderr)
                continue

            hit_phylum = get_taxonomy_details(current_key)

            hit_data_dict = {
                "Seed species": seed_species_name,
                "Common gene symbol of seed sequence": args.seed_gene_name,
                "Hit protein sequence accession": hit_protein_accession,
                "Hit nucleotide sequence accession": hit_nucleotide_accession,
                "Hit species": current_key,
                "Hit phylum": hit_phylum
            }
            best_hits_per_species[current_key] = (e_value, bit_score, hit_data_dict)

    # Collect final data from best_hits_per_species
    for _, _, data_dict in best_hits_per_species.values():
        output_data_list.append(data_dict)

    if not output_data_list:
        print("Warning: No suitable BLAST hits found after filtering.", file=sys.stderr)
        # Create an empty TSV with headers as per spec
    
    # 6. Output
    phylum_tag = args.phylum_taxid if args.phylum_taxid is not None else "all"
    output_filename = f"{args.seed_accession}_{args.seed_gene_name}_{phylum_tag}_{args.blast_type}_hits.tsv"
    print(f"Writing results to {output_filename}...")
    
    df = pd.DataFrame(output_data_list)
    
    # Define columns in the specified order for the output TSV
    output_columns = [
        "Seed species",
        "Common gene symbol of seed sequence",
        "Hit protein sequence accession",
        "Hit nucleotide sequence accession",
        "Hit species",
        "Hit phylum"
    ]

    if df.empty: # If no data, create an empty DataFrame with correct columns
        df = pd.DataFrame(columns=output_columns)
    else:
        # Reorder columns to match specification if DataFrame is not empty
        df = df[output_columns] 

    try:
        df.to_csv(output_filename, sep='\t', index=False, na_rep='N/A')
        print(f"Script finished successfully. Output written to {output_filename}")
    except Exception as e:
        print(f"Error writing TSV file {output_filename}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()