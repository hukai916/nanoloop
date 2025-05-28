# convert bam to json containing per read information
# json structure is like below:
{
  "read_name": "read12345",
  "ref_base_counts": {
    "A": 102,
    "T": 93,
    "C": 84,
    "G": 101
  },
  "bases": [
    {
      "read_pos": 0,
      "ref_pos": 123456,
      "read_base": "A",
      "ref_base": "G",
      "quality": 38
    },
    {
      "read_pos": 1,
      "ref_pos": 123457,
      "read_base": "C",
      "ref_base": "C",
      "quality": 39
    },
    ...
  ]
}


Example code snippet:
  
  import pysam
import json
from collections import Counter

def process_read(read, ref_seq):
    base_info = []
    ref_bases = []
    counts = Counter()

    # get alignment pairs: [(read_pos, ref_pos), ...]
    # with_seq=True gives you read base and ref base
    for read_pos, ref_pos, read_base, ref_base in read.get_aligned_pairs(with_seq=True, matches_only=True):
        if read_pos is None or ref_pos is None:
            continue  # skip deletions or soft clips

        qual = read.query_qualities[read_pos]
        base_info.append({
            "read_pos": read_pos,
            "ref_pos": ref_pos,
            "read_base": read_base,
            "ref_base": ref_base,
            "quality": qual
        })

        if ref_base in "ATCG":
            counts[ref_base] += 1

    return {
        "read_name": read.query_name,
        "ref_base_counts": dict(counts),
        "bases": base_info
    }

def extract_all(bam_path, ref_fasta_path, output_json_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref = pysam.FastaFile(ref_fasta_path)

    with open(output_json_path, "w") as out:
        for read in bam:
            if read.is_unmapped:
                continue
            chrom = bam.get_reference_name(read.reference_id)
            ref_seq = ref.fetch(chrom, read.reference_start, read.reference_end)

            record = process_read(read, ref_seq)
            out.write(json.dumps(record) + "\n")  # NDJSON format


# multiprocessing 
def process_chunk(args):
    chrom, start, end = args
    bam = pysam.AlignmentFile("your.bam", "rb")
    for read in bam.fetch(chrom, start, end):
        if read.reference_start < start:
            continue  # Prevent duplication
        # process read...

