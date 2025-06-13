import gzip
import dask.dataframe as dd
import os
import matplotlib.pyplot as plt
from .constants import ref_mutations

def run_stat_hotspot(args):
  """
  Parse hotspot.bed.gz records and generate statistics/plots
  The input file must be in gz format and must be able to fit into memory
  """
  
  os.makedirs(os.path.abspath(args.output), exist_ok = True)
  
  # Read in hotspot.bed.gz
  with gzip.open(args.hotspot, 'rt') as f:
    header = f.readline()
    if header.startswith('#'):
      header = header.strip().lstrip('#').split('\t')
    else:
      header = None

    if len(header) != 23:
      print(len(header))
      raise ValueError("The input hotspot file must have 23 columns!")
    
    if header is not None:
      df = dd.read_csv(args.hotspot, sep = '\t', names = header, header = None, blocksize = None, skiprows = 1)
    else:
      df = dd.read_csv(args.hotspot, sep = '\t', header = None, blocksize = None)
  
    # Hotspot length distribution
    df['length'] = df['ref_end'] - df['ref_start']
    lengths = df['length'].compute()
    plt.figure(figsize = (8, 6))
    plt.hist(lengths, bins = 50, color = 'skyblue', edgecolor = 'black')
    plt.title('Hotspot Length Distribution', fontsize = 14)
    plt.xlabel('Length', fontsize = 12)
    plt.ylabel('Frequency', fontsize = 12)
    plt.grid(True, linestyle='--', alpha = 0.6)
    plt.tight_layout()
    plt.savefig(args.output + "/hotspot_length_distribution.svg")
    plt.close()
    
    # Per-read hotspot number distribution
    per_read_hotspots = df.groupby('read_id').size().compute()
    counts = per_read_hotspots.value_counts().sort_index()
    plt.figure(figsize = (8, 6))
    plt.bar(counts.index.astype(str), counts.values, color = 'skyblue', edgecolor = 'black')
    plt.title('Hotspot Number per Read', fontsize = 14)
    plt.xlabel('Number of Hotspots', fontsize = 12)
    plt.ylabel('Read Count', fontsize = 12)
    plt.xticks(rotation = 0)
    plt.tight_layout()
    plt.savefig(args.output + "/hotspot_number_per_read.svg")
    plt.close()

    # Number of distinct mutation types per hotspot
    mutation_types = list(ref_mutations.keys())
    df['num_distinct_mutation_types'] = df[mutation_types].apply(
        lambda row: (row > 0).sum(), axis = 1, meta = ('num_distinct_mutation_types', 'int64')
    )
    counts = df['num_distinct_mutation_types'].compute().value_counts().sort_index()
    plt.figure(figsize = (8, 6))
    plt.bar(counts.index.astype(str), counts.values, color = 'skyblue', edgecolor = 'black')
    plt.title('Number of Distinct Mutation Types per Hotspot', fontsize = 14)
    plt.xlabel('Number of Mutation Types', fontsize = 12)
    plt.ylabel('Hotspot Count', fontsize = 12)
    plt.xticks(rotation = 0)
    plt.tight_layout()
    plt.savefig(args.output + "/distinct_mutation_type_per_hotspot.svg")
    plt.close()

  
  
