import gzip
import dask.dataframe as dd
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.colors import ListedColormap
from .constants import ref_mutations

def summarize_read(group, start, end):
  length = end - start + 1

  # Initialize results array as 0 as not covered by read
  results = np.zeros(length, dtype=np.uint8)

  read_ref_start = int(group['read_ref_start'].iloc[0])
  read_ref_end = int(group['read_ref_end'].iloc[0])

  # Step 1: Mark read-covered region as 1
  cov_start = max(read_ref_start, start)
  cov_end = min(read_ref_end, end)

  if cov_start <= cov_end:
    results[(cov_start - start):(cov_end - start + 1)] = 1

  # Step 2: Get all hotspot ranges as arrays
  hs_starts = np.maximum(group['ref_start'].astype(int).values, start)
  hs_ends = np.minimum(group['ref_end'].astype(int).values, end)

  # Step 3: Build a mask array for hotspots
  for hs_start, hs_end in zip(hs_starts, hs_ends):
    if hs_start <= hs_end:
      results[(hs_start - start):(hs_end - start + 1)] = 2

  return results.tolist()

def run_cluster_hotspot(args):
  """
  Given hotspot.bed.gz and a reference range, cluster reads by its hotspot pattern.
  Two types of heatmaps are generated:
  1. Per-hotspot mutation (by count or by fraction; global or along a given reference range) pattern heatmap
  2. Per-read hotspot pattern heatmap along a given reference range
  """
  
  os.makedirs(os.path.abspath(args.output), exist_ok = True)
  
  # Read in hotspot.bed.gz
  # Skip first line if it starts with #
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
  
    # Heatmap 1: global hotspot mutation pattern by count 
    print("Generating global hotspot mutation pattern by count heatmap ...")
    # create a heatmap using ref_mutations.keys() as columns and hotspot mutation counts as rows
    col_mutations = list(ref_mutations.keys())
    df1 = df[col_mutations].compute()
    if df1.empty:
      print(" No hotspot-containing reads, skipped!")
    else:
      # Create heatmap without row labels and with column clustering
      plt.figure(figsize = (12, 8))
      g = sns.clustermap(df1,
                    cmap = 'YlGnBu',
                    col_cluster = True,
                    row_cluster = True,
                    yticklabels = False,
                    xticklabels = True,
                    cbar_pos = (0.02, 0.8, 0.03, 0.18),
                    rasterized = True)
      g.ax_heatmap.set_facecolor('white')
      for _, spine in g.ax_heatmap.spines.items():
        spine.set_visible(False)
      g.ax_heatmap.grid(False)
      plt.suptitle('Heatmap of Mutation Counts Per Hotspot', fontsize = 14)
      plt.savefig(args.output + "/heatmap_per_hotspot_global_count.svg")
      plt.close()
      print("  Done!")
 
    # Heatmap 2: global hotspot mutation pattern by fraction 
    print("Generating global hotspot mutation pattern by fraction heatmap ...")
    region_lengths = (df['ref_end'] - df['ref_start']).compute()
    df2 = df1.div(region_lengths, axis = 0)
    if df2.empty:
      print(" No hotspot-containing reads, skipped!")
    else:
      plt.figure(figsize = (12, 8))
      g = sns.clustermap(df2,
                     cmap = 'YlGnBu',
                     col_cluster = True,
                     row_cluster = True,
                     yticklabels = False,
                     xticklabels = True,
                     cbar_pos = (0.02, 0.8, 0.03, 0.18),
                     rasterized = True)
      g.ax_heatmap.set_facecolor('white')
      for _, spine in g.ax_heatmap.spines.items():
        spine.set_visible(False)
      g.ax_heatmap.grid(False)
      plt.suptitle('Heatmap of Mutation Counts Per Hotspot', fontsize = 14)
      plt.savefig(args.output + "/heatmap_per_hotspot_global_fraction.svg")
      plt.close()
      print("  Done!")

    # Heatmap 3: by-range hotspot mutation pattern by count 
    print("Generating by-range hotspot mutation pattern by count heatmap ...")
    # extract hotspots that fall into the desired ref range
    chr, ref_range = args.range.strip().split(':')
    start, end = map(int, ref_range.split('-'))
    df3 = df[
      (df['ref_chr'] == chr) &
      (df['ref_end'] > start) &
      (df['ref_start'] < end)
    ]

    if df3.compute().empty:
      print(" No hotspot-containing reads fall into the desired range, skipped!")
    else:
      plt.figure(figsize = (12, 8))
      g = sns.clustermap(df3[col_mutations].compute(),
                     cmap = 'YlGnBu',
                     col_cluster = True,
                     row_cluster = True,
                     yticklabels = False,
                     xticklabels = True,
                     cbar_pos = (0.02, 0.8, 0.03, 0.18),
                     rasterized = True)
      g.ax_heatmap.set_facecolor('white')
      for _, spine in g.ax_heatmap.spines.items():
        spine.set_visible(False)
      g.ax_heatmap.grid(False)
      plt.suptitle('Heatmap of Mutation Counts Per Hotspot', fontsize = 14)
      plt.savefig(args.output + "/heatmap_per_hotspot_" + chr + "_" + str(start) + "_" + str(end) + "_count.svg")
      plt.close()
      print("  Done!")
  
    # Heatmap 4: by-range hotspot mutation pattern by fraction 
    print("Generating by-range hotspot mutation pattern by fraction heatmap ...")
    region_lengths = (df3['ref_end'] - df3['ref_start']).compute()
    df4 = df3[col_mutations].compute().div(region_lengths, axis = 0)
    if df4.empty:
      print(" No hotspot-containing reads fall into the desired range, skipped!")
    else:
      plt.figure(figsize = (12, 8))
      g = sns.clustermap(df4,
                     cmap = 'YlGnBu',
                     col_cluster = True,
                     row_cluster = True,
                     yticklabels = False,
                     xticklabels = True,
                     cbar_pos = (0.02, 0.8, 0.03, 0.18),
                     rasterized = True)
      g.ax_heatmap.set_facecolor('white')
      for _, spine in g.ax_heatmap.spines.items():
        spine.set_visible(False)
      g.ax_heatmap.grid(False)
      plt.suptitle('Heatmap of Mutation Counts Per Hotspot', fontsize = 14)
      plt.savefig(args.output + "/heatmap_per_hotspot_" + chr + "_" + str(start) + "_" + str(end) + "_fraction.svg")
      plt.close()
      print("  Done!")
    
  # Heatmap 5: per-read hotspot mutation pattern by range 
    # To use groupby, convert dask to pandas df, and leverage its groupby function, caveat is that it will load the entire df into memory, for this case, it should not be a problem
    print("Generating per-read hotspot mutation pattern by range heatmap ...")
    pd_df = df3.compute()
    result = pd_df.groupby('read_id').apply(summarize_read, start, end)
    
    if result.empty:
      print("No hotspot-containing reads fall into the desired range, skipped!")
    else:
      heatmap_df = pd.DataFrame(result.tolist(), index = result.index)
      heatmap_df.columns = list(range(int(start), int(end) + 1))

      cmap = ListedColormap(['lightgrey', 'lightblue', 'red'])

      g = sns.clustermap(
        heatmap_df,
        cmap = cmap,
        row_cluster = True,
        col_cluster = False,
        xticklabels = 50,
        yticklabels = False,
        figsize = (15, 8),
        cbar_pos = (0.02, 0.8, 0.02, 0.15),
        cbar_kws = {"ticks": [0, 1, 2]},
        rasterized = True # make the figure smaller
      )

      # Optional: Simplify appearance
      g.ax_heatmap.set_facecolor('white')
      for _, spine in g.ax_heatmap.spines.items():
        spine.set_visible(False)
      g.ax_heatmap.grid(False)

      colorbar = g.cax
      colorbar.set_yticks([0, 1, 2])
      colorbar.set_yticklabels(['no read', 'non-hotspot', 'mutation hotspot'])

      g.ax_heatmap.set_xlabel("Position")
      g.ax_heatmap.set_ylabel("Reads")
      g.ax_heatmap.set_title("Per-read Mutation Hotspot Coverage Heatmap")

      plt.savefig(args.output + f"/heatmap_per_read_{chr}_{start}_{end}.svg")
      plt.close()
      print("  Done!")