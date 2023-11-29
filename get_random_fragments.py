import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument("input", help="Input fasta file")
parser.add_argument('--random_length', action='store_true', help='Generate fragments with random length, default: disabled')
parser.add_argument('--allow_duplicates', action='store_true', help='Allow fragments duplicates, default: disabled')
parser.add_argument('--allow_n_content', action='store_true', help=', default: disabled')
parser.add_argument('--length_list', help="List of fragment length values, used if --random_length is not defined, default=[]", nargs='+', default=[])
parser.add_argument('--fragment_length_min', help='Minimum fragment length, --random_length is required, default = 300', default=300)
parser.add_argument('--fragment_length_max', help='Maximum fragment length --random_length is required, default = 1000', default=1000)
parser.add_argument('--fragments_per_record', help='Number of random fragment for a single fasta record, default = 10', default=10)
parser.add_argument('--output', help='Output fasta file', required=True)

args = parser.parse_args()

if not args.random_length and not args.length_list:
    print("Please provide --random_length option or a list of length values with --length_list")
    exit()
elif args.random_length and args.length_list: 
    print("Please select --random_length or a predefined list of length values")
    exit()

fasta_filename = args.input
output_fasta_filename = args.output
length_list = [int(l) for l in args.length_list]

if args.random_length:
    fragment_length_min=int(args.fragment_length_min)
else:
    fragment_length_min = min(length_list)

if args.random_length:
    fragment_length_max=int(args.fragment_length_max)
else:
    fragment_length_max = max(length_list)

fragments_per_record = int(args.fragments_per_record)

if fragment_length_min>fragment_length_max:
	print("fragment_len_min>fragment_length_max")
	exit()
	
fasta_io = SeqIO.parse(fasta_filename, "fasta")

fragments=[]
n_bp_list = ['R','Y','S','W','K','M','B','D','H','V','N']

for record in fasta_io:
    current_id = record.id
    genome=record.seq.upper()
    genome_length = len(str(genome))
    record_fragments=[]
    str_record_fragments=[]
    i = 0
    iteration_limit = fragments_per_record*100

    while len(record_fragments) < fragments_per_record:
        # prevent endless loop
        i+=1
        if i >= iteration_limit:
            break
            
        # select fragment from random position with random length
        fragment_position = random.randint(0, genome_length-1)    
        if genome_length - fragment_position < fragment_length_min:
            continue
        
        if args.random_length:
            fragment_length = min(random.randint(fragment_length_min, fragment_length_max),genome_length - fragment_position)
        else:
            fragment_length = length_list[random.randint(0,len(length_list)-1)]

        fragment = genome[fragment_position:fragment_position+fragment_length]
        
        # discard fragments with Ns
        if not args.allow_n_content and any([bp in n_bp_list for bp in fragment]): 
            continue
        
        if not args.allow_duplicates and (fragment in str_record_fragments):
            continue

        record_fragments.append(SeqRecord(fragment, id=current_id+'_'+str(fragment_position),description = ''))
        str_record_fragments.append(fragment)
    
    # append fragments to output
    fragments+=record_fragments

# write output
with open(output_fasta_filename,"w") as output_handle:
    SeqIO.write(fragments,output_handle,"fasta")        
    print("Result is in " + output_fasta_filename)
