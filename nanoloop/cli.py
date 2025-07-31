from .arg_parser import create_args

def main():
  args = create_args()

  if args.command == 'bam_to_tsv':
    from .bam_to_tsv import run_bam_to_tsv
    run_bam_to_tsv(args)
  elif args.command == 'tsv_to_plot':
    from .tsv_to_plot import run_tsv_to_plot
    run_tsv_to_plot(args)
  elif args.command == 'tsv_to_bed':
    from .tsv_to_bed import run_tsv_to_bed
    run_tsv_to_bed(args)
  elif args.command == 'tsv_to_peak':
    from .tsv_to_peak import run_tsv_to_peak
    run_tsv_to_peak(args)
  elif args.command == 'bam_to_json':
    from .bam_to_json import run_bam_to_json
    run_bam_to_json(args)
  elif args.command == 'filter_json':
    from .filter_json import run_filter_json
    run_filter_json(args)
  elif args.command == 'json_to_hotspot':
    from .json_to_hotspot import run_json_to_hotspot
    run_json_to_hotspot(args)
  elif args.command == 'stat_hotspot':
    from .stat_hotspot import run_stat_hotspot
    run_stat_hotspot(args)
  elif args.command == 'cluster_hotspot':
    from .cluster_hotspot import run_cluster_hotspot
    run_cluster_hotspot(args)