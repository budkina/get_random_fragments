Select random fragments for each record in fasta file

usage: get_random_fragments.py [-h] [--random_length] [--allow_duplicates] [--allow_n_content] [--length_list LENGTH_LIST [LENGTH_LIST ...]] [--fragment_length_min FRAGMENT_LENGTH_MIN] [--fragment_length_max FRAGMENT_LENGTH_MAX] [--fragments_per_record FRAGMENTS_PER_RECORD] --output OUTPUT [input]

positional arguments:

  input                 Input fasta file


optional arguments:

  -h, --help            show this help message and exit

  --random_length       Generate fragments with random length, default: disabled

  --allow_duplicates    Allow fragments duplicates, default: disabled

  --allow_n_content     , default: disabled

  --length_list LENGTH_LIST [LENGTH_LIST ...]       List of fragment length values, used if --random_length is not defined, default=[]

  --fragment_length_min FRAGMENT_LENGTH_MIN       Minimum fragment length, --random_length is required, default = 300

  --fragment_length_max FRAGMENT_LENGTH_MAX       Maximum fragment length --random_length is required, default = 1000

  --fragments_per_record FRAGMENTS_PER_RECORD       Number of random fragment for a single fasta record, default = 10

  --output OUTPUT       Output fasta file
