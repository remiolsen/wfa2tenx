import sys
import re
import gzip
import io
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import cycle
from multiprocessing import Process

WFA = re.compile(" ([ATGC]{20})")
# 10X index kit
SIP02F8 = ["CATGAACA","TCACTCGC","AGCTGGAT","GTGACTTG"]
# "Random" 7 bp oligo
oligo = "TTGCGAG"
r1_qual = ''.join(["A"] * 23)
i1_qual = ''.join(["A"] * 8)


### START DRY-principle violating code-block
def write_read1(wfamap, tenx, read, prefix):
    idx_loop = cycle(SIP02F8)
    with gzip.open(read, 'r') as fz:
        f = io.TextIOWrapper(fz)
        with gzip.open(prefix+"R1_001.fastq.gz", 'wb') as oz:
            for title, seq, qual in FastqGeneralIterator(f):
                tarr = WFA.split(title)
                outln = "@{} 1:N:0:{}\n".format(tarr[0], next(idx_loop))
                outln += "{}{}{}\n".format(tenx[wfamap[tarr[1]]], oligo, seq)
                outln += "+\n"
                outln += "{}{}\n".format(r1_qual, qual)
                oz.write(outln.encode('utf-8'))

def write_read2(wfamap, tenx, read, prefix):
    idx_loop = cycle(SIP02F8)
    with gzip.open(read, 'r') as fz:
        f = io.TextIOWrapper(fz)
        with gzip.open(prefix+"R2_001.fastq.gz", 'wb') as oz:
            for title, seq, qual in FastqGeneralIterator(f):
                tarr = WFA.split(title)
                outln = "@{} 2:N:0:{}\n".format(tarr[0], next(idx_loop))
                outln += "{}\n".format(seq)
                outln += "+\n"
                outln += "{}\n".format(qual)
                oz.write(outln.encode('utf-8'))

def write_i1(wfamap, tenx, read, prefix):
    idx_loop = cycle(SIP02F8)
    with gzip.open(read, 'r') as fz:
        f = io.TextIOWrapper(fz)
        with gzip.open(prefix+"I1_001.fastq.gz", 'wb') as oz:
            for title, seq, qual in FastqGeneralIterator(f):
                tarr = WFA.split(title)
                idx_i = next(idx_loop)
                outln = "@{} 1:N:0:{}\n".format(tarr[0], idx_i)
                outln += "{}\n".format(idx_i)
                outln += "+\n"
                outln += "{}\n".format(i1_qual)
                oz.write(outln.encode('utf-8'))
### END DRY-principle violating code-block

def main(tenxfile, r1, r2, prefix):
    no_reads = 0
    idx = 0
    TENX_BC = []
    wfamap = {} 
    with open(tenxfile, 'r') as f:
        for line in f:
            TENX_BC.append(line.strip())

    with gzip.open(r1, 'r') as fz:
        f = io.TextIOWrapper(fz)
        for title, _, _ in FastqGeneralIterator(f):

            tarr = WFA.split(title)
            if tarr[1] not in wfamap.keys():
                wfamap[tarr[1]] = idx
                idx = idx + 1
            no_reads = no_reads + 1 
        f.close()
    
    p1 = Process(target=write_read1, args=(wfamap, TENX_BC, r1, prefix))
    p2 = Process(target=write_read2, args=(wfamap, TENX_BC, r2, prefix))
    p3 = Process(target=write_i1, args=(wfamap, TENX_BC, r1, prefix))
    p1.start()
    p2.start()
    p3.start()
    p1.join()
    p2.join()
    p3.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate WFA barcodes to 10X barcodes")
    parser.add_argument('--tenx-bc-file', type=str, default="data/4M-with-alts-february-2016.txt", help="Path to text file with 10X barcodes (eg. 4M-with-alts-february-2016.txt)")
    parser.add_argument('--wfa-r1', type=str, required=True, help="Path to read 1")
    parser.add_argument('--wfa-r2', type=str, required=True, help="Path to read 2")
    parser.add_argument('--out-prefix', type=str, default="WFA_OUT_S1_L001_", help="Prefix of the output fastq files (default: WFA_OUT_S1_L001_)")
    args = parser.parse_args()
    sys.exit(main(args.tenx_bc_file, args.wfa_r1, args.wfa_r2, args.out_prefix))

