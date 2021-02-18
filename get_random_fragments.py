import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


parser = argparse.ArgumentParser()
parser.add_argument("input", nargs='?', help="Input fasta file")
parser.add_argument('--fragment_len_min', help='Minimum fragment length, default = 300', default=300)
parser.add_argument('--fragment_len_max', help='Maximum fragment length, default = 1000', default=1000)
parser.add_argument('--fragments_per_record', help='Number of random fragment for a single fasta record, default = 10', default=10)
parser.add_argument('--output', help='Output fasta file')

args = parser.parse_args()

fasta_filename = args.input
output_fasta_filename = args.output
fragment_len_min=int(args.fragment_len_min)
fragment_len_max=int(args.fragment_len_max)
fragments_per_record = int(args.fragments_per_record)

if fragment_len_min>fragment_len_max:
	print("fragment_len_min>fragment_len_max")
	exit()
	
fasta_io = SeqIO.parse(fasta_filename, "fasta")

fragments=[]

for record in fasta_io:
    current_id = record.id
    genome=record.seq.upper()
    genome_length = len(str(genome))
    record_fragments=[]
    i = 0
    iteration_limit = fragments_per_record*100

    while len(record_fragments) < fragments_per_record:
        # prevent endless loop
        i+=1
        if i >= iteration_limit:
            break
            
        # select fragment from random position with random length
        fragment_position = random.randint(0, genome_length-1)    
        if genome_length - fragment_position < fragment_len_min:
            continue
        
        fragment_length = min(random.randint(fragment_len_min, fragment_len_max),genome_length - fragment_position)
        fragment = genome[fragment_position:fragment_position+fragment_length]
        
        # discard fragments with Ns
        if 'N' in fragment: 
            continue
        
        record_fragments.append(SeqRecord(fragment, id=current_id+' '+str(fragment_position),description = ''))
    
    # append fragments to output
    fragments+=record_fragments

# write output
output_handle = open(output_fasta_filename,"w")
SeqIO.write(fragments,output_handle,"fasta")        
print("Result is in " + output_fasta_filename)