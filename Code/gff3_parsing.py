'''
to-do:
1. Tests to make sure everythings works(start codons, lengths, etc.)
2. Test to compare gbk vs gff3 to make sure everything matches 100%
3. Adapting for more than one chromosome contig
4. Function to filter genes we might not want to analyze(genes not %3, letters that aren't ATGC, really long/short genes) 
5. Cleaning up genes that have the same locus tag in multiple rows (ribosomal slippage also)
6. Better commenting
7. Consider genes that overlap/straddle boundaries

'''
import pandas as pd
from Bio import SeqIO

def compile_sequences(gff3_location, fasta_location, upstream_len = 10):
    """
    Inputs: 
    gff3_location - location of gff3 file on your computer as a string
    fasta_location - location of fasta file on your computer as a string
    upstream_len - integer specifiying the length of the upstream sequence
    
    Outputs:
    
    
    description of funtion
    """
    

    #setting up dataframe
    df = pd.read_csv(gff3_location, sep='\t', skiprows=2, header=None)
    df.columns = ["genome_id", "source", "type", "start", "stop", "idk", "strand", "trash", "qualifiers"]
    df = df[df["type"]=='CDS']


    #read in genome 
    genome = list(SeqIO.parse(fasta_location, "fasta"))
    assert len(genome) == 1
    assert len(set(df["genome_id"])) == 1
    genome = genome[0]
    
    assert(df.iloc[0]["genome_id"] == genome.id)

    #getting coding sequences

    for index in df.index:
        #Getting coding sequence
        temp_start = df.loc[index]["start"]
        temp_stop = df.loc[index]["stop"]
        coding_seq = genome.seq[temp_start-1: temp_stop]

        if df.loc[index]["strand"] == "-":
            coding_seq = coding_seq.reverse_complement()


        coding_seq = str(coding_seq)

        df.at[index, "coding_sequence"] = coding_seq

        #getting upstream seequence

        if df.loc[index]["strand"] == "+":

            upstream_seq = genome.seq[temp_start - upstream_len - 1: temp_start -1]

        elif df.loc[index]["strand"] == "-":

            upstream_seq = genome.seq[temp_stop : temp_stop + upstream_len]
            upstream_seq = upstream_seq.reverse_complement()

        upstream_seq = str(upstream_seq)

        df.at[index, "upstream_sequence"] = upstream_seq
        
    return df, genome


