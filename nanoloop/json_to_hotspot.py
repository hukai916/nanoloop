# Parse json looking for "mutation hotspot" in each read, output a BED-like file

import gzip
from multiprocessing import Process, Queue
import os
from .utils import process_json_to_hotspot, chunk_json_file, ref_mutations

def writer_process(queue, include_read_id, include_ref_seq, include_mutation_details, output_path):
  """Dedicated process to write filtered results to gzipped NDJSON file"""
  with gzip.open(output_path, "wt") as gz:
    header = ["#ref_chr", "ref_start", "ref_end"]
    if include_read_id:
      header.append("read_id")
    if include_ref_seq:
      header.append("ref_seq")
    if include_mutation_details:
      header.extend(list(ref_mutations.keys()) + ["read_ref_start", "read_ref_end"])
    gz.write("\t".join(header))
    gz.write("\n")  
    offset = 0
    if include_read_id:
      offset += 1
    if include_ref_seq:
      offset += 1
    while True:
      result = queue.get() # do not reassign queue.get() to queue, otherwise queue will be turned into a list
      if result is None:  # Sentinel to end
        break
      if isinstance(result, Exception):
        raise result
      for hotspot_per_read in result:
        for hotspot in hotspot_per_read:
          line = '\t'.join(map(str, hotspot[:3]))
          if include_read_id:
            line += '\t' + str(hotspot[3])
          if include_ref_seq and include_read_id:
            line += '\t' + str(hotspot[4])
          if include_ref_seq and not include_read_id:
            line += '\t' + str(hotspot[3]) 
          if include_mutation_details:
            line += '\t' + '\t'.join(map(str, hotspot[3 + offset:])) 
          gz.write(line + "\n")

def run_json_to_hotspot(args):
  """
  Parse NDJSON records looking for "mutation hotspot" in each read, output a BED-like file
  """
  os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok = True)
  
  result_queue = Queue()
  writer = Process(target = writer_process, args = (result_queue, args.include_read_id, args.include_ref_seq, args.include_mutation_details, args.output))
  writer.start()
  
  # Process input file in chunks
  active_processes = []
  
  with gzip.open(args.json, 'rt') as f:
    for chunk in chunk_json_file(f):
      while len(active_processes) >= args.ncpus:
        p = active_processes.pop(0)
        p.join()
      
      p = Process(target = process_json_to_hotspot, args = (chunk, args.mutation_type, args.window_size, args.window_step, args.mutation_frac_cutoff, args.include_read_id, args.include_ref_seq, args.include_mutation_details, result_queue))
      active_processes.append(p)
      p.start()
    
    for p in active_processes:
      p.join()
  
  result_queue.put(None)
  writer.join()