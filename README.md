Current progress cosists of four files: blastp_ortholog_finder.py, align_and_tree_for_codeml.py, codeml.ctl, column_scores.py

blastp_ortholog_finder.py - Summary
Functionality
This script is designed to find orthologous genes across different species using BLAST searches. It takes a seed gene (nucleotide sequence) as input, finds its protein sequence, and performs BLAST searches to identify orthologs in other species. When there are multiple hits per species, only the best hit is recorded in the output. The script can handle different types of BLAST searches (blastn, tblastn, blastp, blastx, tblastx) and can be configured to search within specific taxonomic groups.

	Input
	The script requires the following command-line arguments:

	seed_accession: NCBI accession of the seed nucleotide sequence (e.g., NM_001301339)
	seed_gene_name: Common name of the seed gene (e.g., BRCA1)
	blast_type: Type of BLAST search to run ('blastn', 'tblastn', 'blastp', 'blastx', 'tblastx')
	Optional Parameters:
	--email: Email address for NCBI Entrez queries (required if not set in the script)
	--apikey: NCBI API key for higher request limits
	--evalue: E-value threshold for BLAST (default: 0.05)
	--phylum_taxid: NCBI Taxonomy ID to restrict BLAST search to a specific phylum
	--blast_db: Path to a local BLAST database
	--blast_path: Path to BLAST+ executables
	Output
	The script generates a TSV (tab-separated values) file with the following columns:

	Seed species
	Common gene symbol of seed sequence
	Hit protein sequence accession
	Hit nucleotide sequence accession
	Hit species
	Hit phylum
	The output filename follows the pattern: {seed_accession}_{seed_gene_name}_{phylum_tag}_{blast_type}_hits.tsv

align_and_tree_for_codeml.py- Summary
Functionality
This script automates the process of creating codon-aligned sequences and performing phylogenetic analysis. It takes protein and nucleotide sequence accessions as input, performs multiple sequence alignment, and generates phylogenetic trees with dN/dS (ω) ratios. The workflow includes:

	Fetching protein and nucleotide sequences from NCBI
	Performing multiple sequence alignment using MUSCLE
	Creating codon-based alignments using pal2nal
	Building phylogenetic trees
	Calculating dN/dS ratios using PAML's codeml
	Visualizing phylogenetic trees with branch-specific dN/dS values
	Input
	The script requires the following command-line arguments:

	Required:
	--input: Input TSV file containing protein and nucleotide accessions
	--output: Base name for output files
	--email: Your email for NCBI Entrez
	Optional:
	--entrez-api-key: NCBI API key for higher request limits
	--muscle-path: Path to MUSCLE executable (default: "muscle")
	--pal2nal-path: Path to pal2nal.pl script (default: "pal2nal.pl")
	--codeml-path: Path to codeml executable (default: "codeml.exe")
	Output
	The script generates several output files:

	{output}.phy: Codon alignment in PHYLIP format
	{output}.nwk: Phylogenetic tree in Newick format
	{output}_omega_tree.png: Visualization of the phylogenetic tree with dN/dS values
	{output}.mlc: PAML codeml output file with detailed analysis results


codeml.ctl-summary
Functionality -- serves as input file for codeml. Contains command to execute free ratio model, which estimates omega values for each branch 

column_scores.py - Summary
Functionality
This script is designed to assess the conservation of amino acid positions across a set of related protein sequences. It automates the process of fetching sequences from NCBI, performing a multiple sequence alignment with MUSCLE, and then calculating a sum-of-pairs score for each column in the alignment. The score is based on the BLOSUM62 substitution matrix, where a higher score for a column indicates a more highly conserved position.

	Input
	The script requires the following command-line arguments:
	
	input_file: A required text file containing a list of NCBI protein accession numbers, with one accession per line.
	
	-o or --output: The desired name for the output alignment file (default: output.fasta).
	
	Note: The user must manually edit the script to provide two hardcoded values:
	
	A valid email address for NCBI Entrez queries.
	
	The correct file path to the MUSCLE executable on their local machine.
	
	Output
	The script generates two forms of output:
	
	Console Output: A list of sum-of-pairs scores is printed directly to the terminal, with each line showing the column number and its calculated conservation score.
	
	Alignment File: A FASTA file containing the complete multiple sequence alignment. The filename defaults to output.fasta but can be specified with the --output argument.


Dependencies
	Dependencies List For blastp_ortholog_finder.py
		Python Packages:
			argparse (standard library)
			sys (standard library)
			time (standard library)
			pandas (install with pip install pandas)
			BioPython (install with pip install biopython)
			re (standard library)
			os (standard library)
			subprocess (standard library)
		External Tools:
			NCBI BLAST+ tools (if using local BLAST) Download from: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
			NCBI Access: Valid email address (required)
			Optional: NCBI API key (for higher request limits)
	Dependencies List For align_and_tree_for_codeml.py
		Python Packages:
			argparse (standard library)
			os (standard library)
			sys (standard library)
			subprocess (standard library)
			time (standard library)
			tempfile (standard library)
			io (standard library)
			BioPython (install with pip install biopython)
			matplotlib (install with pip install matplotlib)
		External Tools:
			MUSCLE (for multiple sequence alignment) Download from: https://www.drive5.com/muscle/downloads.htm
				Add to system PATH or provide path via --muscle-path
			pal2nal.pl (for codon alignment) Download from: http://www.bork.embl.de/pal2nal/
				Requires Perl to be installed
				Add to system PATH or provide path via --pal2nal-path
			PAML (for dN/dS calculation) Download from: http://abacus.gene.ucl.ac.uk/software/paml.html
				Specifically requires codeml executable
				Add to system PATH or provide path via --codeml-path
			NCBI Access: Valid email address (required)
			Optional: NCBI API key (for higher request limits)
	Dependencies List For column_scores.py
		Python Packages:
			argparse (standard library)
			os (standard library)
			sys (standard library)
			subprocess (standard library)
			tempfile (standard library)
			BioPython (install with pip install biopython)
			numpy (install with pip install numpy)
		External Tools:
			MUSCLE (for multiple sequence alignment)
				Download from: https://www.drive5.com/muscle/downloads.htm
				The path to the executable must be hardcoded into the script.
			NCBI Access: A valid email address is required for Entrez queries and must be hardcoded into the script.
	# Install all Python dependencies
	pip install pandas biopython matplotlib
