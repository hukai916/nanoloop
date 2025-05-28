import os
from .utils import ensure_tsv_index
import gzip
import pandas as pd
def run_tsv_to_bed(args):
  """
  Given tsv.gz from bam_to_tsv, parse out bed-like format for MACS2 peak calling.
  """
  
  # Check input
  if not os.path.exists(args.tsv):
    raise FileNotFoundError(f"The input TSV file {args.tsv} does not exist.")
  ensure_tsv_index(args.tsv)

  if args.type == 'nt_qual':
    with gzip.open(os.path.abspath(args.output), 'wt') as out:
      with gzip.open(args.tsv, 'rt') as f:
        header = f.readline()[1:].strip().split('\t')
        chunk_iter = pd.read_csv(f, sep = '\t', header = None, names = header, chunksize = 10000)
        for chunk in chunk_iter:
          for _, row in chunk.iterrows():
            qual_avg = row['qual_avg']
            qual_avg = 60 - qual_avg # 60 is the max quality score theoretically
            qual_avg = qual_avg * args.scale
            qual_avg = int(qual_avg)
            if qual_avg < 0:
              qual_avg = 0
            bed_line = f"{row['chr']}\t{row['start']}\t{row['end']}\n"
            for _ in range(qual_avg):
              out.write(bed_line)

  elif args.type == 'nt_count':
    with gzip.open(os.path.abspath(args.output), 'wt') as out:
      with gzip.open(args.tsv, 'rt') as f:
        header = f.readline()[1:].strip().split('\t')
        chunk_iter = pd.read_csv(f, sep = '\t', header = None, names = header, chunksize = 10000)
        for chunk in chunk_iter:
          for _, row in chunk.iterrows():
            if row['ref_nt'] == args.ref_nt:
              row["non_x_fraction"] = (1 - row[args.ref_nt] / (row['C'] + row['T'] + row['G'] + row['A'] + row['N'] + 0.01)) * 100 * args.scale 
              row["non_x_fraction"] = int(row["non_x_fraction"])
              
              bed_line = f"{row['chr']}\t{row['start']}\t{row['end']}\n"
              for _ in range(row["non_x_fraction"]):
                out.write(bed_line)