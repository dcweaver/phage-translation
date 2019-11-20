# To do:
1. Edit gff3 parsing to allow for list of chromosomes
2. Edit "Comparative_analysis.ipynb" to write separate ".tsv" files for each host/virus
    * Run on all hosts/viruses currently available
    * Save as taxon_id.tsv in the same folder as the gffs/fastas currently reside
3. Write/move "analysis" code to a new notebook
    * Reads in the previously created ".tsv" files and does all statistics / plotting
    * Add in the FDR p-value method
    * Basically we're separating out the creation of the data it's analysis
4. Write code to be more restrictive in our analyses on the host side of things
    * Remove genes that aren't a multiple of 3
    * Remove "duplicate" locus_tags from the .tsv's
    * Remove potential "phage" genes from hosts
    * Add in codon usage bias measurement for each gene for each organism
    * Write as a new table titled something like "taxon_id.clean.tsv"
5. Write code to perform codon bias corrected statistics for all cases
6. Write E. coli specific code in a new notebook
    * Analyze only essential genes
    * Analyze only highly expressed proteins
7. Make pretty figure/s for publication

**Longer term things to consider**
* Add/use/compare Phanotate phage annotations?
* Use prodigal predicted host genomes?
* Comparison of prodigal/phanotate annotation accuracies on some well known phages?
* Stricter phage clustering?
* Incorporate secondary structure strength?
* Null "min_SD" distribution calculated on like the center of a gene?


