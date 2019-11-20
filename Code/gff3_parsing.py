import pandas as pd
from Bio import SeqIO

def compile_sequences(gff3_location_list, fasta_location_list, upstream_len = 10):
    '''
    Inputs: 
    gff3_location - location of gff3 file on your computer as a string
    fasta_location - location of fasta file on your computer as a string
    upstream_len - integer specifiying the length of the upstream sequence
    
    Outputs:
    df - returns pandas dataframe with information about the genes in the genome that are coding sequences
    genome - returns a list object that contains the entire raw genome sequence from the selected fasta file
    
    Description of funtion:
    
        Reads in a gff3 file and a fasta file and compiles all of the data in the gff3 file into a dataframe that contains information about the genes inside of the genome such as the strand, location of start/stop codons, etc. and filters out all of the non-coding genes. Then, using a fasta file for the same genome, the coding sequence and the region of nucleotides upstream from the start codon are gathered and appended into new columns.  
    '''
    
    df_list = []
    
    for gff3_location, fasta_location in zip(gff3_location_list, fasta_location_list):
        

        #setting up dataframe
        df = pd.read_csv(gff3_location, sep='\t', comment = "#", header=None)
        df.columns = ["genome_id", "source", "type", "start", "stop", "idk", "strand", "trash", "qualifiers"]
        df = df[df["type"]=='CDS']
        df['genome_id'] = df['genome_id'].astype(str)


        #read in genome 
        genome = list(SeqIO.parse(fasta_location, "fasta"))
        #assert statement to ensure only 1 whole genome is being read in, otherwise an error will occur
        assert len(genome) == 1
        #Assert statement to ensure that all of the IDs for the genome are identical, or else an error will be raised
        assert len(set(df["genome_id"])) == 1
        genome = genome[0]
        #Assert statement to ensure that the genome ID from the gff3 matches the genome ID from the fasta file
        assert(df.iloc[0]["genome_id"] == genome.id)

        #getting coding sequences
        for index in df.index:
            #Getting coding sequence
            temp_start = df.loc[index]["start"]
            temp_stop = df.loc[index]["stop"]
            coding_seq = genome.seq[temp_start-1: temp_stop]

            if df.loc[index]["strand"] == "-":
                #must use reverse compliment of bottom strand since start codon is on top strand
                coding_seq = coding_seq.reverse_complement()

            #turn coding sequence into string to be added to DF
            coding_seq = str(coding_seq)
            #Creating new column in dataframe with coding sequence
            df.at[index, "coding_sequence"] = coding_seq

            #getting upstream sequences
            if df.loc[index]["strand"] == "+":

                upstream_seq = genome.seq[temp_start - upstream_len - 1: temp_start -1]

            elif df.loc[index]["strand"] == "-":

                upstream_seq = genome.seq[temp_stop : temp_stop + upstream_len]
                upstream_seq = upstream_seq.reverse_complement()

            upstream_seq = str(upstream_seq)
            #adding upstream sequence column
            df.at[index, "upstream_sequence"] = upstream_seq
            
        df_list.append(df)
        
    if len(df_list) == 1:
        combined_df = df_list[0]
        
    else:
        combined_df = pd.concat(df_list, ignore_index = True) 
        
    return combined_df, genome