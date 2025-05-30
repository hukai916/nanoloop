# convert bam to json containing per read information

import pysam
import json
import gzip
from collections import Counter
import os
import random
import tempfile
import multiprocessing  # for Process and Queue
from multiprocessing import Process, Queue
from .utils import process_bam_chunk_json, ensure_bam_index, chunk_bam

def writer_process(queue, output_path):
  """Dedicated process to write results to gzipped NDJSON file"""
  with gzip.open(output_path, "wt") as gz:
    while True:
      item = queue.get()
      if item is None:  # Sentinel to end
        break
      # Assume item is a dict or list of dicts
      if isinstance(item, list):
        for record in item:
          gz.write(json.dumps(record) + "\n")
      else:
        gz.write(json.dumps(item) + "\n")

def run_bam_to_json(args):
  """
  Given a BAM file and ref fasta, convert it to a JSON file with per read information
  Example json output:
    {
      "read_id": "read12345",
      "ref_chr": "chr1",
      "ref_start": 1,
      "ref_end": 10,
      "ref_seq": "AAATTCCG",
      "ref_A_count": 30,
      "ref_T_count": 40,
      "ref_C_count": 50,
      "ref_G_count": 30,
      "CtoA": {
        "ref_pos": [],
        "read_pos": [],
        "base_quality": []
        },
      "CtoT": {
        "ref_pos": [],
        "read_pos": [],
        "base_quality": []
      }, ...
      }
    }
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
  random.shuffle(chunks)

  os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok = True)
    
  # Process chunks in batches based on number of CPUs
  with tempfile.TemporaryDirectory() as temp_dir:
      print('Temp folder path: ', temp_dir)
      
      # Create a queue for results and start writer process
      result_queue = Queue()
      writer = Process(target = writer_process, args = (result_queue, args.output))
      writer.start()
      
      # Process chunks in batches of args.ncpus
      for i in range(0, len(chunks), args.ncpus):
        batch = chunks[i:i + args.ncpus]
        processes = []
        
        for chunk in batch:
          p = Process(target = process_bam_chunk_json,
                      args = (chunk[0], chunk[1], chunk[2], args.bam, args.ref, temp_dir, result_queue)
)
          processes.append(p)
          p.start()
      
        for p in processes:
          p.join()
      
      # Signal writer to end and wait for it to finish
      result_queue.put(None)
      writer.join()


