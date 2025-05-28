import os
import dask.dataframe as dd
import pandas as pd

def merge_overlapping_intervals(df):
  if df.empty:
    return pd.DataFrame(columns=['start', 'end'])
  
  df = df.reset_index(drop = True)
  merged = []
  current_start, current_end = df.loc[0, 'start'], df.loc[0, 'end']
  for i in range(1, len(df)):
    row = df.loc[i]
    if row['start'] <= current_end:
      current_end = max(current_end, row['end'])
    else:
      merged.append((current_start, current_end))
      current_start, current_end = row['start'], row['end']
  merged.append((current_start, current_end))
  return pd.DataFrame(merged, columns = ['start', 'end'])

def merge_nearby_intervals(df, window_size):
  if df.empty:
    return pd.DataFrame(columns = ['start', 'end'])

  merged = []
  current_start, current_end = df.loc[0, 'start'], df.loc[0, 'end']
  
  for i in range(1, len(df)):
    row = df.loc[i]
    if row['start'] - current_end <= window_size:
      current_end = max(current_end, row['end'])
    else:
      merged.append((current_start, current_end))
      current_start, current_end = row['start'], row['end']

  merged.append((current_start, current_end))
  
  return pd.DataFrame(merged, columns = ['start', 'end'])

def group_func_nt_qual(x, low_qual_cutoff = 20, frac_cutoff = 0.3, window_size = 20, merge_nearby_peaks = True):
  '''
  Each chr is a group, perform sliding window approach within each group.
  '''
  
  chr = x['chr'].iloc[0]
  print("Processing chromosome: ", chr)
  # min, max = x['start'].min(), x['start'].max()
  if low_qual_cutoff == 10:
    low_group = x['qual_0_10']
  elif low_qual_cutoff == 20:
    low_group = x['qual_0_10'] + x['qual_10_20']
  elif low_qual_cutoff == 30:
    low_group = x['qual_0_10'] + x['qual_10_20'] + x['qual_20_30']
  elif low_qual_cutoff == 40:
    low_group = x['qual_0_10'] + x['qual_10_20'] + x['qual_20_30'] + x['qual_30_40']
    
  x['low_frac'] = (low_group) / (x['qual_0_10'] + x['qual_10_20'] + x['qual_20_30'] + x['qual_30_40'] + x['qual_40_above'])
  x = x.sort_values('start').reset_index(drop = True)
  x['low_frac_rolling_avg'] = x['low_frac'].rolling(window = window_size, min_periods = 1).mean() 
  x = x[x['low_frac_rolling_avg'] >= frac_cutoff]
  df = merge_overlapping_intervals(x)
  if merge_nearby_peaks:
    df = merge_nearby_intervals(df, window_size * 2)

  df['chr'] = chr
  df = df[['chr', 'start', 'end']]
  return df

def group_func_nt_count(x, conversion_cutoff = 0.1, window_size = 20, merge_nearby_peaks = True):
    '''
    Each chr is a group, perform sliding window approach within each group.
    '''
    
    chr = x['chr'].iloc[0]
    print("Processing chromosome: ", chr)
    
    x = x.copy()
    c_positions = x[x['ref_nt'] == 'C'].index.tolist()
    
    # For each non-C position, find the next C position and use its values
    for i in range(len(x)):
      if x.loc[i, 'ref_nt'] != 'C':
        next_c_pos = None
        for c_pos in c_positions:
          if c_pos > i:
            next_c_pos = c_pos
            break
    
        if next_c_pos is not None:
          x.loc[i, 'C'] = x.loc[next_c_pos, 'C']
          x.loc[i, 'T'] = x.loc[next_c_pos, 'T']
          
    # x = x[x['ref_nt'] == 'C']
    x['conversion_frac'] = x['T'] / (x['T'] + x['C'] + 0.01)
    x['conversion_frac_rolling_avg'] = x['conversion_frac'].rolling(window = window_size, min_periods = 1).mean() 
    x = x[x['conversion_frac_rolling_avg'] >= conversion_cutoff]
    
    df = merge_overlapping_intervals(x)
    if merge_nearby_peaks:
      df = merge_nearby_intervals(df, window_size * 2)
    df['chr'] = chr
    df = df[['chr', 'start', 'end']]
    return df

def run_tsv_to_peak(args):
  """
  Given tsv.gz from bam_to_tsv, call peaks using a sliding window approach.
  """
  
  # Check input
  if not os.path.exists(args.tsv):
    raise FileNotFoundError(f"The input TSV file {args.tsv} does not exist.")

  import dask
  dask.config.set(scheduler = 'processes', num_workers = args.ncpus) # parallel groupby may take large memory
  
  if args.type == 'nt_qual':
    print("Reading tsv.gz ...")
    cols = ['chr', 'start', 'end', 'ref_nt', 'qual_0_10', 'qual_10_20', 'qual_20_30', 'qual_30_40', 'qual_40_above', 'qual_avg']
    dtypes = {'chr': str, 'start': int, 'end': int, 'ref_nt': str, 'qual_0_10': int, 'qual_10_20': int, 'qual_20_30': int, 'qual_30_40': int, 'qual_40_above': int, 'qual_avg': float}
    
    if args.tsv.endswith('.gz'):
      df = dd.read_csv(args.tsv, compression = 'gzip', sep = '\t', names = cols, dtype = dtypes, header = None, blocksize = None, skiprows = 1)
    elif args.tsv.endswith('.tsv'):
      df = dd.read_csv(args.tsv, sep = '\t', names = cols, dtype = dtypes, header = None, blocksize = '500MB', skiprows = 1)
    
    # Group by chromosome and apply the rolling average window approach to call peaks
    df_res = df.groupby('chr').apply(group_func_nt_qual, low_qual_cutoff = args.low_qual_cutoff, frac_cutoff = args.frac_cutoff, window_size = args.window_size, merge_nearby_peaks = args.merge_nearby_peaks, meta = {'chr': str, 'start': int, 'end': int}).compute()
  elif args.type == 'nt_count':
    print("Reading tsv.gz ...")
    cols = ['chr', 'start', 'end', 'ref_nt', 'A', 'T', 'C', 'G', 'N']
    dtypes = {'chr': str, 'start': int, 'end': int, 'ref_nt': str, 'A': int, 'T': int, 'C': int, 'G': int, 'N': int}
    
    if args.tsv.endswith('.gz'):
      df = dd.read_csv(args.tsv, compression = 'gzip', sep = '\t', names = cols, dtype = dtypes, header = None, blocksize = None, skiprows = 1)
    elif args.tsv.endswith('.tsv'):
      df = dd.read_csv(args.tsv, sep = '\t', names = cols, dtype = dtypes, header = None, blocksize = '500MB', skiprows = 1)
    
    # Group by chromosome and apply the rolling average window approach to call peaks
    df_res = df.groupby('chr').apply(group_func_nt_count, conversion_cutoff = args.conversion_cutoff, window_size = args.window_size, merge_nearby_peaks = args.merge_nearby_peaks, meta = {'chr': str, 'start': int, 'end': int}).compute()

  # Drop peaks that are too short 
  df_res = df_res[df_res['end'] - df_res['start'] >= args.min_peak_length]
  df_res.to_csv(os.path.abspath(args.output), sep = '\t', index = False, compression = 'gzip', header = False)
  print("Peaks called and saved to ", os.path.abspath(args.output), sep = '')