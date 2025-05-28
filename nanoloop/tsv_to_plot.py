import os
from .utils import ensure_tsv_index, extract_range_from_tsv
import seaborn as sns
import numpy as np  
import matplotlib.pyplot as plt
def run_tsv_to_plot(args):
  """
  Given tsv.gz from bam_to_tsv, parse out desired region and plot.
  """
  
  # Check input
  if not os.path.exists(args.tsv):
    raise FileNotFoundError(f"The input TSV file {args.tsv} does not exist.")
  ensure_tsv_index(args.tsv)

  # Extract desired region 
  df = extract_range_from_tsv(args.tsv, args.range)
  
  if args.type == 'nt_count':
    # Calculate ref_nt C converted to all other nucleotides ratio
    if (args.mode == "ratio"):
      y_label = 'Allele frequency'
      df['C_to_T_ratio'] = df['T'] / (df['C'] + df['T'] + 0.01)
      df['G_to_A_ratio'] = df['A'] / (df['G'] + df['A'] + 0.01) 
    elif (args.mode == "count"):
      y_label = 'Allele count'

    df_C = df[df['ref_nt'] == 'C']
    df_G = df[df['ref_nt'] == 'G']

    _, ax = plt.subplots()
    plot_kwargs = dict(edgecolor = 'black', linewidth = 0.5, s = 12, alpha = 0.75)
    
    if args.mode == "ratio":
      sns.scatterplot(data = df_C, ax = ax, x = 'start', y = 'C_to_T_ratio', label = 'T/(C + T) at ref C sites', **plot_kwargs, color = '#0072B2',)
      sns.scatterplot(data = df_G, ax = ax, x = 'start', y = 'G_to_A_ratio', label = 'A/(A + G) at ref G sites', **plot_kwargs, color = '#D55E00')
    elif args.mode == "count":
      sns.scatterplot(data = df_C, ax = ax, x = 'start', y = 'T', label = 'T count at ref C sites', **plot_kwargs, color = '#0072B2')
      sns.scatterplot(data = df_G, ax = ax, x = 'start', y = 'A', label = 'A count at ref G sites', **plot_kwargs, color = '#D55E00')
    
    # Add moving average line to ax
    window_size = 10
    df_C = df_C.copy() # avoid SettingWithCopyWarning
    df_G = df_G.copy()
    if args.mode == "ratio":
      df_C['C_to_T_ratio_rolling_avg'] = df_C['C_to_T_ratio'].rolling(window = window_size, min_periods = 1).mean()
      df_G['G_to_A_ratio_rolling_avg'] = df_G['G_to_A_ratio'].rolling(window = window_size, min_periods = 1).mean()
      sns.lineplot(data = df_C, ax = ax, x = 'start', 
                  y = 'C_to_T_ratio_rolling_avg',
                  label = 'T/(C + T) at ref C sites (rolling avg)',
                  color = '#0072B2', linewidth = 0.4, alpha = 0.6)
      sns.lineplot(data = df_G, ax = ax, x = 'start', 
                  y = 'G_to_A_ratio_rolling_avg',
                  label = 'A/(A + G) at ref G sites (rolling avg)',
                  color = '#D55E00', linewidth = 0.4, alpha = 0.6)
    elif args.mode == "count":
      df_C['T_rolling_avg'] = df_C['T'].rolling(window = window_size, min_periods = 1).mean()
      df_G['A_rolling_avg'] = df_G['A'].rolling(window = window_size, min_periods = 1).mean()
      sns.lineplot(data = df_C, ax = ax, x = 'start', 
                  y = 'T_rolling_avg', 
                  label = 'T count at ref C sites (rolling avg)',
                  color = '#0072B2', linewidth = 0.4, alpha = 0.6)
      sns.lineplot(data = df_G, ax = ax, x = 'start', 
                  y = 'A_rolling_avg', 
                  label = 'A count at ref G sites (rolling avg)',
                  color = '#D55E00', linewidth = 0.4, alpha = 0.6)
    
    # Add GC% line to ax2 if args.add_gc is True
    if (args.add_gc):
      window_size = 25
      ax2 = ax.twinx()
      
      df['is_GC'] = df['ref_nt'].isin(['G', 'C'])
      df['smoothed_GC_percent'] = df['is_GC'].rolling(window = window_size, min_periods = 1).mean()
    
      sns.lineplot(data = df, ax = ax2, x = 'start', y = 'smoothed_GC_percent', color = '#000000', linewidth = 0.3, alpha = 0.45)
      ax2.set_ylabel('GC%')
    
    ax.set_xlabel('Position')
    ax.set_ylabel(y_label)
    plt.savefig(args.output)
  elif args.type == 'nt_qual':
    # Calculate qual ratio for all levels
    if args.mode == "ratio":
      y_label = 'Qualilty frequency'
      df['0_10_ratio'] = df['qual_0_10'] / (df['qual_0_10'] + df['qual_10_20'] + df['qual_20_30'] + df['qual_30_40'] + df['qual_40_above'])
      df['10_20_ratio'] = df['qual_10_20'] / (df['qual_0_10'] + df['qual_10_20'] + df['qual_20_30'] + df['qual_30_40'] + df['qual_40_above'])
      df['20_30_ratio'] = df['qual_20_30'] / (df['qual_0_10'] + df['qual_10_20'] + df['qual_20_30'] + df['qual_30_40'] + df['qual_40_above'])
      df['30_40_ratio'] = df['qual_30_40'] / (df['qual_0_10'] + df['qual_10_20'] + df['qual_20_30'] + df['qual_30_40'] + df['qual_40_above'])
      df['40_above_ratio'] = df['qual_40_above'] / (df['qual_0_10'] + df['qual_10_20'] + df['qual_20_30'] + df['qual_30_40'] + df['qual_40_above'])
      y_cols = ['0_10_ratio', '10_20_ratio', '20_30_ratio', '30_40_ratio', '40_above_ratio']
    elif args.mode == 'count':
      y_label = 'Quality count'
      y_cols = ['qual_0_10', 'qual_10_20', 'qual_20_30', 'qual_30_40', 'qual_40_above']
    
    _, ax = plt.subplots()
    plot_kwargs = dict(edgecolor = 'black', linewidth = 0.5, s = 12, alpha = 0.75)
    n_bins = 5
    cmap = plt.cm.get_cmap('summer')
    colors = [cmap(i) for i in np.linspace(0, 1, n_bins)][::-1]
    labels = ['0-10', '10-20', '20-30', '30-40', '40+']
    
    bottom = np.zeros(len(df))
    for y_col, color, label in zip(y_cols, colors, labels):
      ax.bar(df['start'], df[y_col], bottom = bottom, color = color, label = label, width = 1.0)
      bottom += df[y_col].values

    ax.set_xlabel('Position')
    ax.set_ylabel(y_label)
    ax.legend(title='Quality Bin')
    
    if (args.add_gc):
      window_size = 25
      ax2 = ax.twinx()
      
      df['is_GC'] = df['ref_nt'].isin(['G', 'C'])
      df['smoothed_GC_percent'] = df['is_GC'].rolling(window = window_size, min_periods = 1).mean()
    
      sns.lineplot(data = df, ax = ax2, x = 'start', y = 'smoothed_GC_percent', color = '#000000', linewidth = 0.3, alpha = 0.45)
      ax2.set_ylabel('GC%')
      
    elif args.add_qual_avg:
      window_size = 25
      ax2 = ax.twinx()
      df['qual_avg_rolling_avg'] = df['qual_avg'].rolling(window = 25, min_periods = 1).mean()
      sns.lineplot(data = df, ax = ax2, x = 'start', y = 'qual_avg_rolling_avg', color = '#000000', linewidth = 0.3, alpha = 0.45)
      ax2.set_ylabel('Qual Avg (rolling avg)')
    
    plt.savefig(os.path.abspath(args.output))