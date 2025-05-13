import os
import pysam
import random
import tempfile
import multiprocessing
from functools import partial
from .utils import ensure_bam_index, chunk_bam, process_bam_chunk_nt_qual, process_bam_chunk_nt_count, parse_filename
import glob
                            
def run_bam_to_tsv(args):
  """
  Given a BAM file and ref fasta, convert it to a TSV file with nucleotide quality counts or nucleotide variant counts depending on args.type.
  Output a TSV file with the following columns (--type nt_qual):
    - chr
    - start (0-based)
    - end
    - ref_nt
    - qual_0_10: nt qualities falling into this range [)
    - qual_10_20 [)
    - qual_20_30 [)
    - qual_30_40 [)
    - qual_40_above [)
  Output a TSV file with the following columns (--type nt_count):
    - chr
    - start (0-based)
    - end
    - ref_nt
    - A (number of 'A' nt)
    - T (number of 'T' nt)
    - C (number of 'C' nt)
    - G (number of 'G' nt)
    - N (number of all other nt)
  """
  
  # Check input
  if not os.path.exists(args.bam):
    raise FileNotFoundError(f"The input BAM file {args.bam} does not exist.")
  ensure_bam_index(args.bam)
  if not os.path.exists(args.ref):
    raise FileNotFoundError(f"The reference FASTA file {args.ref} does not exist.")

  # Chunk bam for multiprocessing
  chunks = []
  with pysam.AlignmentFile(args.bam, "rb") as bamfile:
    ref_names = (bamfile.references)
    ref_lengths = [bamfile.get_reference_length(x) for x in ref_names]

    for ref_name, ref_length in zip(ref_names, ref_lengths):
      chunks.extend(chunk_bam(ref_name, ref_length, args.ncpus))
    random.shuffle(chunks) # ensure each parallel job load is even (some chr region may have more coverage)
    
  # Process each chunk in parallel
  with tempfile.TemporaryDirectory() as temp_dir:
    print('Temp folder path: ', temp_dir)
    with multiprocessing.Pool(processes = args.ncpus) as pool:
      if args.type == 'nt_qual':
        header = '\t'.join(['#chr', 'start', 'end', 'ref_nt', 'qual_0_10', 'qual_10_20', 'qual_20_30', 'qual_30_40', 'qual_40_above', 'qual_avg']) + '\n'
        func = partial(process_bam_chunk_nt_qual, bam_path = args.bam, ref_path = args.ref, temp_dir = temp_dir)
      elif args.type == 'nt_count':
        header = '\t'.join(['#chr', 'start', 'end', 'ref_nt', 'A', 'T', 'C', 'G', 'N']) + '\n'
        func = partial(process_bam_chunk_nt_count, bam_path = args.bam, ref_path = args.ref, temp_dir = temp_dir)
      pool.starmap(func, chunks)
      
      tsv_files = glob.glob(os.path.join(temp_dir, '*.tsv'))
      tsv_sorted = sorted(tsv_files, key = lambda f: (parse_filename(f)[0], parse_filename(f)[1]))
      if not os.path.exists(os.path.dirname(args.output)):
        os.makedirs(os.path.dirname(args.output))
      with pysam.BGZFile(args.output, 'w') as out:
        out.write(header.encode())

        for path in tsv_sorted:
          with open(path) as f:
            for line in f:
                if line.strip():
                  out.write(line.encode())