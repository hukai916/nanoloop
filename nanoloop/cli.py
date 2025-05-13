from .arg_parser import create_args
from .bam_to_tsv import run_bam_to_tsv
from .tsv_to_plot import run_tsv_to_plot
from .tsv_to_bed import run_tsv_to_bed
from .tsv_to_peak import run_tsv_to_peak

def main():
  args = create_args()

  if args.command == 'bam_to_tsv':
    run_bam_to_tsv(args)
  elif args.command == 'tsv_to_plot':
    run_tsv_to_plot(args)
  elif args.command == 'tsv_to_bed':
    run_tsv_to_bed(args)
  elif args.command == 'tsv_to_peak':
    run_tsv_to_peak(args)

if __name__ == '__main__':
  main()