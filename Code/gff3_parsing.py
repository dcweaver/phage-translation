import pandas as pd
from Bio import SeqIO

def compile_sequences(gff3_location, fasta_location, upstream_len = 10):
    """
    
    """
    

    #setting up dataframe
    df = pd.read_csv(gff3_location, sep='\t', skiprows=2, header=None)
    df.columns = ["genome_id", "source", "type", "start", "stop", "idk", "strand", "trash", "qualifiers"]
    df = df[df["type"]=='CDS']


    #read in genome 
    genome = list(SeqIO.parse(fasta_location, "fasta"))
    assert len(genome) == 1
    genome = genome[0]


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


