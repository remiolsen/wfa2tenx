import sys
import re
import io
import os
import subprocess
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from itertools import cycle

WFA = re.compile(" ([ATGCN]{20})")
# 10X index kit
SIP02F8 = ["CATGAACA","TCACTCGC","AGCTGGAT","GTGACTTG"]
# "Random" 7 bp oligo
oligo = "TTGCGAG"
r1_nbc = ''.join(["N"] * 16)
r1_qual = ''.join(["A"] * 23)
i1_qual = ''.join(["A"] * 8)
gz_buf = 131072
try:
    subprocess.Popen(["pigz"], stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
except OSError as e:
    if e.errno == os.errno.ENOENT:
        print("Could not find 'pigz' command in system environment")
        sys.exit(1)
    else:
        raise


### START DRY-principle violating code-block
def write_read1(wfamap, tenx, read, prefix, procs):
    idx_loop = cycle(SIP02F8)
    p_split = str(int(int(procs) / 2))
    with subprocess.Popen(["pigz", "-d", "-c", "-p", p_split, read],
            stdout=subprocess.PIPE, bufsize=gz_buf) as fzi:
        fi = io.TextIOWrapper(fzi.stdout, write_through=True)
        with open(prefix+"R1_001.fastq.gz", 'wb') as ofile:
            with subprocess.Popen(["pigz", "-c", "-p", p_split],
                    stdin=subprocess.PIPE, stdout=ofile, bufsize=gz_buf, close_fds=False) as oz:
                for title, seq, qual in FastqGeneralIterator(fi):
                    tarr = title.split()
                    try:
                        wfa_bc = tenx[wfamap[tarr[1]]]
                    except KeyError:
                        wfa_bc = r1_nbc
                        #print("No barcode in file {} in read {}, inserting Ns".format(read, title), file=sys.stderr)
                    except IndexError:
                        raise

                    outln = "@{} 1:N:0:{}\n".format(tarr[0], next(idx_loop))
                    outln += "{}{}{}\n".format(wfa_bc, oligo, seq)
                    outln += "+\n"
                    outln += "{}{}\n".format(r1_qual, qual)
                    oz.stdin.write(outln.encode('utf-8'))

def write_read2(wfamap, tenx, read, prefix, procs):
    idx_loop = cycle(SIP02F8)
    p_split = str(int(int(procs) / 2))
    with subprocess.Popen(["pigz", "-d", "-c", "-p", p_split, read],
            stdout=subprocess.PIPE, bufsize=gz_buf) as fzi:
        fi = io.TextIOWrapper(fzi.stdout, write_through=True)
        with open(prefix+"R2_001.fastq.gz", 'wb') as ofile:
            with subprocess.Popen(["pigz", "-c", "-p", p_split],
                    stdin=subprocess.PIPE, stdout=ofile, bufsize=gz_buf, close_fds=False) as oz:
                for title, seq, qual in FastqGeneralIterator(fi):
                    tarr = title.split()
                    outln = "@{} 2:N:0:{}\n".format(tarr[0], next(idx_loop))
                    outln += "{}\n".format(seq)
                    outln += "+\n"
                    outln += "{}\n".format(qual)
                    oz.stdin.write(outln.encode('utf-8'))

def write_i1(wfamap, tenx, read, prefix, procs):
    idx_loop = cycle(SIP02F8)
    p_split = str(int(int(procs) / 2))
    with subprocess.Popen(["pigz", "-d", "-c", "-p", p_split, read],
            stdout=subprocess.PIPE, bufsize=gz_buf) as fzi:
        fi = io.TextIOWrapper(fzi.stdout, write_through=True)
        with open(prefix+"I1_001.fastq.gz", 'wb') as ofile:
            with subprocess.Popen(["pigz", "-c", "-p", p_split],
                    stdin=subprocess.PIPE, stdout=ofile, bufsize=gz_buf, close_fds=False) as oz:
                for title, seq, qual in FastqGeneralIterator(fi):
                    tarr = title.split()
                    idx_i = next(idx_loop)
                    outln = "@{} 1:N:0:{}\n".format(tarr[0], idx_i)
                    outln += "{}\n".format(idx_i)
                    outln += "+\n"
                    outln += "{}\n".format(i1_qual)
                    oz.stdin.write(outln.encode('utf-8'))
### END DRY-principle violating code-block

def main(tenxfile, r1, r2, prefix, total_processes, minbc):
    idx = 0
    tenx_c = 0
    TENX_BC = []
    wfamap = {}
    wfamap_c = {}
    with open(tenxfile, 'r') as f:
        for line in f:
            TENX_BC.append(line.strip())
            tenx_c += 1

    with subprocess.Popen(["pigz", "-d", "-c", "-p", total_processes, r1],
            stdout=subprocess.PIPE) as fz:
        with io.TextIOWrapper(fz.stdout, write_through=True) as f:
            for title, _, _ in FastqGeneralIterator(f):
                tarr = WFA.split(title)
                if tarr[1] not in wfamap.keys():
                    #wfamap[tarr[1]] = idx
                    wfamap_c[tarr[1]] = wfamap_c.get(tarr[1], 0) + 1
                    #idx = idx + 1

    nrbc = len(wfamap_c)
    for key, count in wfamap_c.items():
        assert idx <= tenx_c, "Found more barcodes that available for 10X ({}), try using --min-bc argument".format(tenx_c)
        if count >= minbc:
            wfamap[key] = idx
            idx += 1
    
    print("Found:\t{} WFA barcodes".format(len(wfamap_c.keys())))
    print("Made:\t{} 10X barcodes".format(idx))
    write_read1(wfamap, TENX_BC, r1, prefix, total_processes)
    write_read2(wfamap, TENX_BC, r2, prefix, total_processes)
    write_i1(wfamap, TENX_BC, r1, prefix, total_processes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Translate WFA barcodes to 10X barcodes")
    parser.add_argument('--tenx-bc-file', '-t',type=str, default="{}/supernova/tenkit/lib/python/tenkit/barcodes/4M-with-alts-february-2016.txt".format(os.path.dirname(os.path.realpath(__file__))), help="Path to text file with 10X barcodes (eg. 4M-with-alts-february-2016.txt)")
    parser.add_argument('--wfa-r1', '-1', type=str, required=True, help="Path to read 1")
    parser.add_argument('--wfa-r2', '-2', type=str, required=True, help="Path to read 2")
    parser.add_argument('--out-prefix', '-o', type=str, default="WFA_OUT_S1_L001_", help="Prefix of the output fastq files (default: WFA_OUT_S1_L001_)")
    parser.add_argument('--processes', '-p', type=str, default=2, help="Number of processes to spawn (default: 2)")
    parser.add_argument('--min-bc', '-m', type=int, default=1, help="Minumum barcode multiplicity to include it")
    args = parser.parse_args()
    sys.exit(main(args.tenx_bc_file, args.wfa_r1, args.wfa_r2, args.out_prefix, args.processes, args.min_bc))

