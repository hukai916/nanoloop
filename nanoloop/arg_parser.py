import argparse
from .utils import validate_range

def str2bool(v):
  if isinstance(v, bool):
    return v
  if v.lower() in ('yes', 'true', 't', 'y', '1'):
    return True
  elif v.lower() in ('no', 'false', 'f', 'n', '0'):
    return False
  else:
    raise argparse.ArgumentTypeError('Boolean value expected.')

def create_args():
  parser = argparse.ArgumentParser(description = 'nanoloop')
  subparsers = parser.add_subparsers(dest = 'command')
  
  # Add version command
  try:
    from . import __version__
    version = __version__
  except ImportError:
    version = "unknown"
  parser.add_argument('--version',
                     action = 'version',
                     version = version)

  # subcommand: bam_to_tsv
  parser_bam_to_tsv = subparsers.add_parser('bam_to_tsv',
                                            help = 'Parse BAM file and output TSV file')
  parser_bam_to_tsv.add_argument('--bam', 
                                  type = str, 
                                  required = True,
                                  help = 'Path to input BAM file')
  parser_bam_to_tsv.add_argument('--ref', 
                                  type = str, 
                                  required = True,
                                  help = 'Path to reference fasta file')
  parser_bam_to_tsv.add_argument('--output', 
                                  type = str, 
                                  required = True,
                                  help='Path to output tsv.gz file')
  parser_bam_to_tsv.add_argument('--type', 
                                  type = str, 
                                  required = True,
                                  choices = ['nt_count', 'nt_qual'],
                                  help = 'Type of tsv file to output: "nt_count" or "nt_qual", the former is for nucleotide species count at each site and the latter is for nucleotide read quality count at each site')
  parser_bam_to_tsv.add_argument('--ncpus', 
                                  type = int, 
                                  default = 1,
                                  help = 'Number of CPUs to use')
  
  # subcommand: tsv_to_plot
  parser_tsv_to_plot = subparsers.add_parser('tsv_to_plot',
                                            help = 'Parse TSV file and output plot')
  parser_tsv_to_plot.add_argument('--tsv', 
                                  type = str, 
                                  required = True,
                                  help = 'Path to input TSV file')
  parser_tsv_to_plot.add_argument('--type', 
                                  type = str, 
                                  required = True,
                                  choices = ['nt_count', 'nt_qual'],
                                  help = 'Type of input TSV file being used: "nt_count" or "nt_qual", the former is for nucleotide species count at each site and the latter is for nucleotide read quality count at each site')
  parser_tsv_to_plot.add_argument('--range', 
                                  type = str, 
                                  required = True,
                                  help = 'Range to plot: e.g. "chr1:1000000-1100000"')
  parser_tsv_to_plot.add_argument('--mode', 
                                  type = str, 
                                  choices = ['ratio', 'count'],
                                  default = 'ratio', 
                                  help = 'Plot mode: "ratio" or "count", the former is for C_to_T fraction and the latter is for raw count at each site')
  parser_tsv_to_plot.add_argument('--add_gc', 
                                  type = str2bool, 
                                  default = False,
                                  help = 'Add GC content to plot, note that GC content is calculated from the TSV file, which may not cover all reference positions')
  parser_tsv_to_plot.add_argument('--add_qual_avg', 
                                  type = str2bool, 
                                  default = False,
                                  help = 'Relevant when --type is "nt_qual", add rolling average of quality score to plot, mutual exclusive to --add_gc')
  parser_tsv_to_plot.add_argument('--output', 
                                  type = str, 
                                  required = True,
                                  help='Path to output pdf plot file, must be ended with ".pdf" or ".png"')
  
  # subcommand: tsv_to_bed
  parser_tsv_to_bed = subparsers.add_parser('tsv_to_bed',
                                            help = 'Convert TSV file to BED format for MACS3 peak calling. For "--type nt_qual": regions with lower quality scores will generate more simulated read tags, creating peaks in those regions. For "--type nt_count": the number of simulated read tags is proportional to the fraction of C converted to T at each reference C position')
  parser_tsv_to_bed.add_argument('--tsv', 
                                  type = str, 
                                  required = True,
                                  help = 'Path to input TSV file')
  parser_tsv_to_bed.add_argument('--output',  
                                  type = str, 
                                  required = True,
                                  help = 'Path to output bed.gz file')
  parser_tsv_to_bed.add_argument('--scale', 
                                  type = float, 
                                  default = 1.0,
                                  help = 'Scale factor for quality score, value larger than 1 will make the read coverage difference simulated from quality score larger, and the resulting bed will have more simulated read tags in it')
  parser_tsv_to_bed.add_argument('--type', 
                                  type = str, 
                                  required = True,
                                  choices = ['nt_count', 'nt_qual'],
                                  help = 'Type of input TSV file being used: "nt_count" or "nt_qual", the former is for nucleotide species count at each reference C site and the latter is for nucleotide read quality count at each reference site')
  
  # subcommand: tsv_to_peak
  parser_tsv_to_peak = subparsers.add_parser('tsv_to_peak',
                                            help = 'Call peaks using a sliding window approach')
  parser_tsv_to_peak.add_argument('--tsv', 
                                  type = str, 
                                  required = True,
                                  help = 'Path to input TSV file, both .tsv and .tsv.gz are supported, if .tsv.gz, the file size must be less than your memory size')
  parser_tsv_to_peak.add_argument('--type', 
                                  type = str, 
                                  required = True,
                                  choices = ['nt_qual', 'nt_count'],
                                  help = 'Type of input TSV file being used: "nt_qual" or "nt_count", the former is for nucleotide read quality count at each reference site and the latter is for nucleotide species count at each reference C site')
  parser_tsv_to_peak.add_argument('--low_qual_cutoff', 
                                  type = int, 
                                  default = 30,
                                  choices = [10, 20, 30, 40],
                                  help = 'Low quality cutoff, relevant when --type is "nt_qual", choose from 10, 20, 30, and 40, quality score below this cutoff is considered low quality, default to 30')
  parser_tsv_to_peak.add_argument('--frac_cutoff', 
                                  type = float, 
                                  default = 0.4,
                                  help = 'Fraction cutoff for peak calling, relevant when --type is "nt_qual"; given a window size, the fraction of low quality reads in the window is calculated, if the fraction is larger than the cutoff, the window is considered a peak, default to 0.4')
  parser_tsv_to_peak.add_argument('--conversion_cutoff', 
                                  type = float, 
                                  default = 0.05,
                                  help = 'Fraction cutoff for miscalled T at reference C site, relevant when --type is "nt_count", used to call peaks at reference C sites')
  parser_tsv_to_peak.add_argument('--min_peak_length', 
                                  type = int, 
                                  default = 0,
                                  help = 'Minimum length of peaks to keep, default to 50, detected peaks shorter than this value will be discarded')
  parser_tsv_to_peak.add_argument('--window_size', 
                                  type = int, 
                                  default = 50,
                                  help = 'Window size for peak calling, default to 50') 
  parser_tsv_to_peak.add_argument('--merge_nearby_peaks', 
                                  type = str2bool, 
                                  default = True,
                                  help = 'Merge nearby peaks, peaks that are within the (2 * window size) will be merged into one peak, default to True')
  parser_tsv_to_peak.add_argument('--output', 
                                  type = str, 
                                  required = True,
                                  help = 'Path to output bed.gz file')
  parser_tsv_to_peak.add_argument('--ncpus', 
                                  type = int, 
                                  default = 1,
                                  help = 'Number of CPUs to use, larger number may increase the memory usage significantly')
  
  args = parser.parse_args()
  
  if args.command == 'bam_to_tsv':
    if not args.output.endswith('.tsv.gz'):
      parser.error('The --output path must end with ".tsv.gz"!')
  
  if args.command == 'tsv_to_plot':
    if not args.tsv.endswith('.tsv.gz'):
      parser.error('The --tsv path must end with .tsv.gz and must be created with bgzip!')
    validate_range(args.range)
    if args.type == 'nt_qual':
      if args.add_gc and args.add_qual_avg:
        parser.error('Only one of --add_gc and --add_qual_avg can be specified at a time!')
    if args.type == 'nt_count' and args.add_qual_avg:
      parser.error('--add_qual_avg is not applicable for --type nt_count!')

  if args.command == 'tsv_to_bed':
    if not args.tsv.endswith('.tsv.gz'):
      parser.error('The --tsv path must end with .tsv.gz and must be created with bgzip!')
    if not args.output.endswith('.bed.gz'):
      parser.error('The --output path must end with .bed.gz!')
  
  if args.command == 'tsv_to_peak':
    if not args.output.endswith('.bed.gz'):
      parser.error('The --output path must end with .bed.gz!')

  return args