Select random fragments for each record in fasta file

usage: get_random_fragments.py [-h] [--fragment_len_min FRAGMENT_LEN_MIN] [--fragment_len_max FRAGMENT_LEN_MAX] [--fragments_per_record FRAGMENTS_PER_RECORD] [--output OUTPUT] [input]

positional arguments:

  input Input fasta file
  
optional arguments:

  -h, --help  show this help message and exit
  
  --fragment_len_min FRAGMENT_LEN_MIN Minimum fragment length, default = 300
  
  --fragment_len_max FRAGMENT_LEN_MAX Maximum fragment length, default = 1000
  
  --fragments_per_record FRAGMENTS_PER_RECORD Number of random fragment for a single fasta record, default = 10
  
  --output OUTPUT Output fasta file
  
