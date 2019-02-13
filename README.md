# Cut-and-run
Cut-and-Run: A New Mechanism by which V(D)J Recombination causes Genome Instability

Python scripts used in the analysis of LAM-HTGTS data in the above paper

1. RSS script (Find_Translocations_to_RSSs.py)
Finds RSS motifs (using RSS site, https://www.itb.cnr.it/rss/index.html) within a specified window of translocations. The script takes all TLX files generated using the Translocation Sequencing Pipeline (Alt Lab; http://robinmeyers.github.io/transloc_pipeline/) within the script directory as the input and returns a csv with the number of translocations containing an RSS of the specified type per experiment. The window size around the translocation breakpoint to be searched for RSSs can be adjusted with the -w/--window flag (default = 45 bp) and the RSS type (12 or 23) selected with the -RSS/--RSS_type flag (default = 12). The location of a genome assembly (e.g. hg19) in fasta format must also be supplied using the -g/--genome flag.   

2. Translocations in genes script (Find_translocations_to_genes.py)
Finds the number of translocations observed in a supplied list of genes. The script requires a list of gene symbols (one gene name per line), a genome annotation file in gtf format and a list of translocation in TLX format (all TLX files in the script directory will be processed). The script outputs a csv detailing the number of translocations observed in each of the specified genes per experiment. Translocation counts can be normalised to the total number of translocations using the -n/--normalise flag (default is False).     

3. Chromosome distribution script (Find_chromosome_distribution.py)
Finds the number of translations on each chromosome. This script analyses all TLX files within the directory and returns the number of translocations observed on each chromosome in a csv file and a pie chart describing the distribution.

4. Patient breakpoint overlap script (Find_translocations_to_breakpoint.py)
Finds overlaps between a list of breakpoints (in bed format) and translocations. The script processes all TLX files in the script directory and returns the breakpoints overlapped by translocations in csv format. The breakpoints file must be specified first and the window around the translocation breakpoints to be considered for overlap with the supplied breakpoints can be adjusted with the -ow/--overlap_window flag (default = 50 bp). 
