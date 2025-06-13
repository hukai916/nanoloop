# Filter json files 

import json
import gzip
from multiprocessing import Process, Queue
import os
from .utils import process_filter_json, chunk_json_file

def writer_process(queue, output_path):
  """Dedicated process to write filtered results to gzipped NDJSON file"""
  with gzip.open(output_path, "wt") as gz:
    while True:
      item = queue.get()
      if item is None:  # Sentinel to end
        break
      if isinstance(item, Exception):
        raise item
      # Write filtered records
      for record in item:
        gz.write(json.dumps(record) + "\n")

def run_filter_json(args):
  """
  Filter NDJSON records based on specified criteria
  """
  os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok = True)
  
  result_queue = Queue()
  writer = Process(target = writer_process, args = (result_queue, args.output))
  writer.start()
  
  # Process input file in chunks
  active_processes = []
  
  with gzip.open(args.json, 'rt') as f:
    for chunk in chunk_json_file(f):
      while len(active_processes) >= args.ncpus:
        p = active_processes.pop(0)
        p.join()
      
      p = Process(target = process_filter_json, args = (chunk, args.by, args.count_cutoff, args.frac_cutoff, args.base_quality_cutoff, result_queue))
      active_processes.append(p)
      p.start()
    
    for p in active_processes:
      p.join()
  
  result_queue.put(None)
  writer.join()