from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC 
from Bio.SeqFeature import FeatureLocation, ExactPosition, SeqFeature
import sys 
import pandas as pd
import os

####################
# uses a bed cov file and fasta file to pick out contiguos regions of coverage
#####################


cov_table = sys.argv[1]
sequence = sys.argv[2]
outfile = str(os.getcwd()+"/covered_seqeunces_")
#print outfile

# coordinate grabber and fasta outputer 

contig_list = []

def contig_id_grabber(fasta, contigs):
        #pulls full seqeunces from a fast file, requires an empty list
        for seq_record in SeqIO.parse(open(fasta, "r"), "fasta"):
                contigs.append(seq_record.id)
        print contigs    


def seq_grabber(fasta,coord,seq_id):
        covered_regions = []
        contig_count = 1
        for seq_record in SeqIO.parse(open(fasta, "r"), "fasta"):
                if seq_record.id == seq_id:
                        for seq_start, seq_end in zip(coord[::2], coord[1::2]):
                                covered_regions.append(SeqRecord(seq_record.seq[seq_start:seq_end], id = str("%s_%d" % (seq_record.id, contig_count)), description = str("bases %d -> %d" % (seq_start, seq_end))))
                                #print covered_regions
                                contig_count += 1
                else:
                        print "%s seqeunce does not contain same sequence IDs and coverage file" % seq_record.id                                        
                SeqIO.write(covered_regions, str("%s_%s.fasta" % (outfile+seq_record.id)), "fasta")

# cov parse

contig_id_grabber(sequence,contig_list)

coord = []
#print coord
with open(cov_table) as cov_table:
        cov_df = pd.read_csv(cov_table, index_col=False, header=None, names = ['chrm', 'base', 'cov'], sep ="\t", lineterminator='\n')
        for ids in contig_list:
                for index, row in cov_df.iterrows():
                        if cov_df.iloc[index]['chrm'] == ids:
                                if cov_df.iloc[index]['cov'] == 0 and cov_df.iloc[index-1]['cov'] == 0:
                                        pass
                                elif cov_df.iloc[index]['cov'] > 0 and cov_df.iloc[index-1]['cov'] > 0:
                                        pass
                                elif index == 0 and cov_df.iloc[index-1]['cov'] > 0:
                                        coord.append(row['base'])
                                else:
                                        coord.append(row['base'])       
                        else:
                                pass
                                #print "fasta seqeunce does not contain same sequence IDs and coverage file"            
                print coord
                if not coord:
                        pass
                else:    
                        seq_grabber(sequence,coord,ids)
                coord = []                              

