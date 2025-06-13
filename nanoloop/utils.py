import pysam
import os
import io
import pandas as pd
import numpy as np
import re
from collections import Counter
import gzip
from .constants import ref_mutations
import json

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
    temp_file = os.path.join(temp_dir, '{}_{}_{}.tsv'.format(ref_name.replace("_", "."), start_pos, end_pos))
  
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
    temp_file = os.path.join(temp_dir, '{}_{}_{}.tsv'.format(ref_name.replace("_", "."), start_pos, end_pos))
  
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

# Functions for bam_to_json.py


def process_per_read_alignment(read, ref_name, ref_seq, chunk_start): 
  ref_start = read.reference_start
  ref_end = read.reference_end
  aligned_bases = {
    'read_id': read.query_name,
    'ref_chr': ref_name,
    'ref_start': ref_start,
    'ref_end': ref_end,
    'ref_seq': ref_seq,
    'ref_A_count': ref_seq.count('A'),
    'ref_T_count': ref_seq.count('T'),
    'ref_C_count': ref_seq.count('C'),
    'ref_G_count': ref_seq.count('G')
  }
  
  # Initialize ref_mutations to 0: CtoA, CtoG, CtoT, CtoN, AtoC, AtoG, AtoT, AtoN, GtoC, GtoA, GtoT, GtoN, TtoC, TtoA, TtoG, TtoN, etc.
  for ref_mut in ref_mutations:
    aligned_bases[ref_mut] = {
      'ref_pos': [],
      'read_pos': [],
      'base_quality': []
    }
  
  for read_pos, ref_pos in read.get_aligned_pairs(matches_only = True):
    if read_pos is None or ref_pos is None:
      continue
    
    read_nt = read.query_sequence[read_pos]
    read_qual = read.query_qualities[read_pos]
    ref_index = ref_pos - chunk_start
    if ref_index < 0 or ref_index >= len(ref_seq):
      continue # skip rare edge cases
    ref_nt = ref_seq[ref_index]
    if ref_nt == 'A':
      if read_nt == 'C':
        aligned_bases['AtoC']['ref_pos'].append(ref_pos)
        aligned_bases['AtoC']['read_pos'].append(read_pos)
        aligned_bases['AtoC']['base_quality'].append(read_qual)
      elif read_nt == 'G':
        aligned_bases['AtoG']['ref_pos'].append(ref_pos)
        aligned_bases['AtoG']['read_pos'].append(read_pos)
        aligned_bases['AtoG']['base_quality'].append(read_qual)
      elif read_nt == 'T':
        aligned_bases['AtoT']['ref_pos'].append(ref_pos)
        aligned_bases['AtoT']['read_pos'].append(read_pos)
        aligned_bases['AtoT']['base_quality'].append(read_qual)
      elif read_nt == 'N':
        aligned_bases['AtoN']['ref_pos'].append(ref_pos)
        aligned_bases['AtoN']['read_pos'].append(read_pos)
        aligned_bases['AtoN']['base_quality'].append(read_qual)
    elif ref_nt == 'C':
      if read_nt == 'A':
        aligned_bases['CtoA']['ref_pos'].append(ref_pos)
        aligned_bases['CtoA']['read_pos'].append(read_pos)
        aligned_bases['CtoA']['base_quality'].append(read_qual)
      elif read_nt == 'G':
        aligned_bases['CtoG']['ref_pos'].append(ref_pos)
        aligned_bases['CtoG']['read_pos'].append(read_pos)
        aligned_bases['CtoG']['base_quality'].append(read_qual)
      elif read_nt == 'T':
        aligned_bases['CtoT']['ref_pos'].append(ref_pos)
        aligned_bases['CtoT']['read_pos'].append(read_pos)
        aligned_bases['CtoT']['base_quality'].append(read_qual)
      elif read_nt == 'N':
        aligned_bases['CtoN']['ref_pos'].append(ref_pos)
        aligned_bases['CtoN']['read_pos'].append(read_pos)
        aligned_bases['CtoN']['base_quality'].append(read_qual)
    elif ref_nt == 'G':
      if read_nt == 'A':
        aligned_bases['GtoA']['ref_pos'].append(ref_pos)
        aligned_bases['GtoA']['read_pos'].append(read_pos)
        aligned_bases['GtoA']['base_quality'].append(read_qual)
      elif read_nt == 'C':
        aligned_bases['GtoC']['ref_pos'].append(ref_pos)
        aligned_bases['GtoC']['read_pos'].append(read_pos)
        aligned_bases['GtoC']['base_quality'].append(read_qual)
      elif read_nt == 'T':
        aligned_bases['GtoT']['ref_pos'].append(ref_pos)
        aligned_bases['GtoT']['read_pos'].append(read_pos)
        aligned_bases['GtoT']['base_quality'].append(read_qual)
      elif read_nt == 'N':
        aligned_bases['GtoN']['ref_pos'].append(ref_pos)
        aligned_bases['GtoN']['read_pos'].append(read_pos)
        aligned_bases['GtoN']['base_quality'].append(read_qual)
    elif ref_nt == 'T':
      if read_nt == 'A':
        aligned_bases['TtoA']['ref_pos'].append(ref_pos)
        aligned_bases['TtoA']['read_pos'].append(read_pos)
        aligned_bases['TtoA']['base_quality'].append(read_qual)
      elif read_nt == 'C':
        aligned_bases['TtoC']['ref_pos'].append(ref_pos)
        aligned_bases['TtoC']['read_pos'].append(read_pos)
        aligned_bases['TtoC']['base_quality'].append(read_qual)
      elif read_nt == 'G':
        aligned_bases['TtoG']['ref_pos'].append(ref_pos)
        aligned_bases['TtoG']['read_pos'].append(read_pos)
        aligned_bases['TtoG']['base_quality'].append(read_qual)
      elif read_nt == 'N':
        aligned_bases['TtoN']['ref_pos'].append(ref_pos)
        aligned_bases['TtoN']['read_pos'].append(read_pos)
        aligned_bases['TtoN']['base_quality'].append(read_qual)

  return aligned_bases

# chunk, args.bam, args.ref, temp_dir, result_queue
def process_bam_chunk_json(ref_name, start_pos, end_pos, bam_path, ref_path, temp_dir, result_queue):
  try:
    print(f"Processing chunk: {ref_name}_{start_pos}_{end_pos}")

    with pysam.AlignmentFile(bam_path, "rb") as bam, pysam.FastaFile(ref_path) as ref:
      results = []
      ref_seq = ref.fetch(ref_name, start_pos, end_pos) # chunk ref sequence
      for read in bam.fetch(ref_name, start_pos, end_pos):
        if read.is_unmapped:
          continue
        if read.query_sequence is None:
          continue
        if read.reference_start < start_pos:
          continue # avoid duplicating reads 
        
        aligned_bases = process_per_read_alignment(read, ref_name, ref_seq, start_pos)
        results.append(aligned_bases)
    result_queue.put(results)

  except Exception as e:
    raise e
    # result_queue.put(e)

# functions for filter_json.py
def chunk_json_file(file, chunk_size = 1000):
  """Generator to yield chunks of lines from a file"""
  chunk = []
  for line in file:
    chunk.append(line)
    if len(chunk) >= chunk_size:
      yield chunk
      chunk = []
  if chunk:
    yield chunk

def process_filter_json(chunk, by, count_cutoff, frac_cutoff, base_quality_cutoff, queue):
  """Process a chunk of NDJSON records and put filtered results in queue"""
  
  try:
    filtered_records = []
    for line in chunk:
      record = json.loads(line)
      
      # print(by, count_cutoff, frac_cutoff, base_quality_cutoff)
      for ref_mut in ref_mutations:
        # for each ref_mut, filter out mutations with base quality below base_quality_cutoff, do it for ref_pos, read_pos, base_quality
        record[ref_mut]['ref_pos'] = [x for x in record[ref_mut]['ref_pos'] if record[ref_mut]['base_quality'][record[ref_mut]['ref_pos'].index(x)] >= base_quality_cutoff]
        record[ref_mut]['read_pos'] = [x for x in record[ref_mut]['read_pos'] if record[ref_mut]['base_quality'][record[ref_mut]['read_pos'].index(x)] >= base_quality_cutoff]
        record[ref_mut]['base_quality'] = [x for x in record[ref_mut]['base_quality'] if x >= base_quality_cutoff]

      total_count = 0
      if by == "count":
        for ref_mut in ref_mutations:
          total_count += len(record[ref_mut]['ref_pos'])
        if total_count >= count_cutoff:
          filtered_records.append(record)
      elif by == "frac":
        for ref_mut in ref_mutations:
          total_count += len(record[ref_mut]['ref_pos'])
        if total_count >= frac_cutoff * len(record['ref_seq']):
          filtered_records.append(record)
    queue.put(filtered_records)
  except Exception as e:
    raise e

# Functions for json_to_hotspot.py
# chunk, 
  # How to define a hotspot:
  # mutation_type: "all, "CtoT", "CtoT|CtoG", ...
  # window_size: 25
  # window_step: 5
  # mutation_frac_cutoff: 0.5
    # if in a window, at least 50% of the bases are mutated, that window will be treated as a "hotspot" window. The fraction is calculated by dividing the number of interested muations (mutation_type) by the total number of that base in the reference window. For example, if mutation_type is "all", then it is calculated by dividing the number of mutations by the total number of bases in the reference window; if mutation_type is "CtoT|GtoA", then it is calculated by dividing the number of CtoT and GtoA mutations by the total number of C and G bases in the reference window.
def process_json_to_hotspot(chunk, mutation_type, window_size, window_step, mutation_frac_cutoff, include_read_id, include_ref_seq, include_mutation_details, queue): 
  print("Processing chunk_id", id(chunk), "...")

  if mutation_type == "all":
    mutation_types = ref_mutations.keys()
  else:
    mutation_types = mutation_type.split("|")
  ref_bases = list(set([x.split("to")[0] for x in mutation_types]))
  
  hotspots_all = []
  try:
    for line in chunk:
      hotspots_per_read = []
      hotspot_windows = []
      record = json.loads(line)
      
      # Find all mutated positions
      mutated_positions = []
      for ref_mut in mutation_types:
        for ref_pos in record[ref_mut]['ref_pos']:
          mutated_positions.append(ref_pos)

      # Use a sliding window along the reference sequence to find hotspots in each read
      for ref_pos in range(0, len(record['ref_seq']), window_step):
        ref_seq_window = record['ref_seq'][ref_pos:ref_pos + window_size]
        ref_s, ref_d = ref_pos + record['ref_start'], ref_pos + record['ref_start'] + window_size
        ref_base_count = sum([ref_seq_window.count(base) for base in ref_bases])
        ref_chr = record['ref_chr']

        # count how many read mutations in the window
        mutation_count = 0
        for ref_mut in mutation_types:
          for ref_pos in record[ref_mut]['ref_pos']:
            if ref_pos in range(ref_s, ref_d):
              mutation_count += 1
        
        mutation_frac = mutation_count / ref_base_count if ref_base_count > 0 else 0 
        if mutation_frac >= mutation_frac_cutoff:
          hotspot_windows.append([ref_chr, ref_s, ref_d])

      # Merge hotspot windows if within window_step
      for hotspot in hotspot_windows:
        if not hotspots_per_read:
          hotspots_per_read.append(hotspot)
        else:
          if hotspot[1] - hotspots_per_read[-1][2] < window_step:
            hotspots_per_read[-1][2] = hotspot[2]
          else:
            hotspots_per_read.append(hotspot)
      
      # Add read_id, ref_seq to each hotspot
      if include_read_id:
        for hotspot in hotspots_per_read:
          hotspot.append(record['read_id'])
      if include_ref_seq:
        for hotspot in hotspots_per_read:
          hotspot_seq = record['ref_seq'][hotspot[1] - record['ref_start']:hotspot[2] - record['ref_start']]
          hotspot.append(hotspot_seq)
      if include_mutation_details:
        for hotspot in hotspots_per_read:
          for ref_mut in ref_mutations:
            mutation_count = 0
            for ref_pos in record[ref_mut]['ref_pos']:
              if ref_pos in range(hotspot[1], hotspot[2]):
                mutation_count += 1
            hotspot.append(mutation_count)
          # Also add aligned reference range:
          hotspot.append(record['ref_start'])
          hotspot.append(record['ref_end'])
      
      hotspots_all.append(hotspots_per_read)

    if not hotspots_all == []:
      queue.put(hotspots_all)
  except Exception as e:
    raise e
    
