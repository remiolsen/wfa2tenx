# wfa2tenx

Highly experimental code to translate WFA to tenx barcodes

# Requirements

* [pigz](http://www.zlib.net/pigz/)
* python >= 2.7
* biopython 
  * `pip install biopython`
* tenx barcode file `4M-with-alts-february-2016.txt` 
  * included as in supernova as a submodule, ie. `git clone --recursive wfa2tenx.git`

# Usage

```bash
usage: wfa2tenx.py [-h] [--tenx-bc-file TENX_BC_FILE] --wfa-r1 WFA_R1 --wfa-r2
                   WFA_R2 [--out-prefix OUT_PREFIX] [--processes PROCESSES]
                   [--min-bc MIN_BC]

Translate WFA barcodes to 10X barcodes

optional arguments:
  -h, --help            show this help message and exit
  --tenx-bc-file TENX_BC_FILE, -t TENX_BC_FILE
                        Path to text file with 10X barcodes (eg. 4M-with-alts-
                        february-2016.txt)
  --wfa-r1 WFA_R1, -1 WFA_R1
                        Path to read 1
  --wfa-r2 WFA_R2, -2 WFA_R2
                        Path to read 2
  --out-prefix OUT_PREFIX, -o OUT_PREFIX
                        Prefix of the output fastq files (default:
                        WFA_OUT_S1_L001_)
  --processes PROCESSES, -p PROCESSES
                        Number of processes to spawn (default: 2)
  --min-bc MIN_BC, -m MIN_BC
                        Minumum barcode multiplicity to include it
  ```
