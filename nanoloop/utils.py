import pysam
import os
import io
from statistics import mean
import pandas as pd
import numpy as np
from statistics import stdev
import re
from collections import Counter
import gzip

def validate_range(user_range):
    range_pattern = r"^(\w+):(\d+)-(\d+)$"
    match = re.match(range_pattern, user_range)
    if match == None:
      raise ValueError(f"--range format is incorrect. Expected format: chrname:start-end")
    
    _, start, end = match.groups()
    start, end = int(start), int(end)
    if start >= end:
      raise ValueError(f"--range: start position must be less than end position.")

def ensure_bam_index(bam_file):
    bam_index = bam_file + '.bai'
    if not os.path.exists(bam_index):
        print(f"Index not found for {bam_file}. Creating index ...")
        pysam.index(bam_file)
        
def ensure_tsv_index(tsv_file):
    tsv_index = tsv_file + '.tbi'
    if not os.path.exists(tsv_index):
        print(f"Index not found for {tsv_file}. Creating index ...")
        pysam.tabix_index(tsv_file, preset = "bed", force = True)
        
def parse_filename(path):
    filename = os.path.basename(path)
    match = re.match(r"([a-zA-Z0-9]+)_([\d]+)_([\d]+)", filename)
    if match:
        chr = match.group(1)
        start = int(match.group(2))
        return chr, start
    return None, None

def chunk_bam(ref_name, ref_length, num_chunks):
  chunk_size = (ref_length + num_chunks - 1) // num_chunks 
  chunks = [(ref_name, start, min(start + chunk_size, ref_length)) 
            for start in range(0, ref_length, chunk_size)]  
  return chunks

def extract_range_from_tsv(file_path, region):
  with gzip.open(file_path, 'rt') as f:
    header = f.readline().strip() 
    if header.startswith('#'):
      header = header[1:].split('\t')
    else:
      n_col = len(header.split('\t'))
      header = [f'col{i + 1}' for i in range(n_col)]

  with pysam.TabixFile(file_path) as tabix_file:
    print(file_path)
    lines = list(tabix_file.fetch(region = region))
    if not lines:
      raise ValueError('Warnning: no data found in the specified region.')
      
    data = '\n'.join(lines)

    df = pd.read_csv(io.StringIO(data), sep = '\t', header = None)
    df.columns = header 
    
    return df

def process_bam_chunk_nt_qual(ref_name, start_pos, end_pos, bam_path, ref_path, temp_dir):
  print(f"Processing chunk: {ref_name}_{start_pos}_{end_pos}")

  results = []
  with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(ref_path) as ref:
    temp_file = os.path.join(temp_dir, '{}_{}_{}.tsv'.format(ref_name, start_pos, end_pos))
  
    for pileup_column in bam.pileup(ref_name, start_pos, end_pos, min_base_quality = 0, max_depth = 2_147_483_647):
      pos = pileup_column.pos
      if pos < start_pos or pos >= end_pos:
        continue 
      
      ref_base = ref.fetch(ref_name, pos, pos + 1).upper()
      ref_position = pos # 0-based
      
      quals = []
      for pileup_read in pileup_column.pileups:
        read_pos = pileup_read.query_position
        read_qual = pileup_read.alignment.query_qualities[read_pos] if read_pos is not None else 0 
        quals.append(read_qual)

      bins = [0, 10, 20, 30, 40, np.inf]
      bin_labels = ['qual_0_10', 'qual_10_20', 'qual_20_30', 'qual_30_40', 'qual_40_above']
      bin = pd.cut(quals, bins = bins, labels = bin_labels, right = False)
      
      results.append(
          {
            'chr': ref_name,
            'start': ref_position,
            'end': ref_position + 1,
            'ref_base': ref_base,
            bin_labels[0]: bin.value_counts()[bin_labels[0]],
            bin_labels[1]: bin.value_counts()[bin_labels[1]],
            bin_labels[2]: bin.value_counts()[bin_labels[2]],
            bin_labels[3]: bin.value_counts()[bin_labels[3]],
            bin_labels[4]: bin.value_counts()[bin_labels[4]],
            'qual_avg': np.mean(quals)
          }
      )
      
      if len(results) >= 1000000:
        df = pd.DataFrame(results)
        if not pd.io.common.file_exists(temp_file):
            df.to_csv(temp_file, sep = '\t', index = False, mode = 'w', header = False)
        else:
            df.to_csv(temp_file, sep = '\t', index = False, mode = 'a', header = False)
        results = []

  df = pd.DataFrame(results)
  if not pd.io.common.file_exists(temp_file):
    df.to_csv(temp_file, sep = '\t', index= False, mode = 'w', header = False)
  else:
    df.to_csv(temp_file, sep = '\t', index = False, mode = 'a', header = False)
  
  return (temp_file)

def process_bam_chunk_nt_count(ref_name, start_pos, end_pos, bam_path, ref_path, temp_dir):
  print(f"Processing chunk: {ref_name}_{start_pos}_{end_pos}")

  results = []
  with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(ref_path) as ref:
    temp_file = os.path.join(temp_dir, '{}_{}_{}.tsv'.format(ref_name, start_pos, end_pos))
  
    for pileup_column in bam.pileup(ref_name, start_pos, end_pos, min_base_quality = 0, max_depth = 2_147_483_647):
      pos = pileup_column.pos
      if pos < start_pos or pos >= end_pos:
        continue 
      
      ref_base = ref.fetch(ref_name, pos, pos + 1).upper()
      ref_position = pos # 0-based
      
      nt = []
      for pileup_read in pileup_column.pileups:
        read_pos = pileup_read.query_position
        read_nt = pileup_read.alignment.query_sequence[read_pos].upper() if read_pos is not None else 'N'
        nt.append(read_nt)

      bins = ['A', 'T', 'C', 'G', 'N']
      counts = Counter(nt)
      bin = {nt: counts.get(nt, 0) for nt in bins}
            
      results.append(
          {
            'chr': ref_name,
            'start': ref_position,
            'end': ref_position + 1,
            'ref_base': ref_base,
            bins[0]: bin[bins[0]],
            bins[1]: bin[bins[1]],
            bins[2]: bin[bins[2]],
            bins[3]: bin[bins[3]],
            bins[4]: bin[bins[4]]
          }
      )
      
      if len(results) >= 1000000:
        df = pd.DataFrame(results)
        if not pd.io.common.file_exists(temp_file):
            df.to_csv(temp_file, sep = '\t', index = False, mode = 'w', header = False)
        else:
            df.to_csv(temp_file, sep = '\t', index = False, mode = 'a', header = False)
        results = []

  df = pd.DataFrame(results)
  if not pd.io.common.file_exists(temp_file):
    df.to_csv(temp_file, sep = '\t', index= False, mode = 'w', header = False)
  else:
    df.to_csv(temp_file, sep = '\t', index = False, mode = 'a', header = False)
  
  return (temp_file)
