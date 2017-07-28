# Running the VkJk-Seq pipeline

The data associated with the VkJk publication will soon be available from GEO under the accession number GSE101606. Here we provide two test files of the WT_pre-B_VkJk-seq sample (the first 10,000 sequences) to enable reproduction of how the pipeline was run. The files are called: 

 * Read 1: `VkJk_test_V-region.fastq.gz`
 * Read 2: `VkJk_test_J-region.fastq.gz`
 
 
### Step 1: Quality/adapter trimming and 5' end clipping

In the current implementation of the VkJk-Seq pipeline we assume that the basecall qualities of the first 4 base pairs of Read 2 (the J-reads) are quite poor. We thus clip R2 by 4 bp from the 5' end. In addition, Trim Galore detects and removes read-through adapter contamination and poor quality basecalls from the 5' ends. The trimming command is:

```
trim_galore --clip_r2 4 --paired VkJk_test_V-region.fastq.gz VkJk_test_J-region.fastq.gz
```
The output files are called: `VkJk_test_V-region_val_1.fq.gz` and `VkJk_test_J-region_val_2.fq.gz`.


### Step 2: Mispriming correction

The mispriming correction can be run on the quality/adapter trimmed J-region file(s) like so:

```
./mispriming_correction_mouse_IgK.pl VkJk_test_J-region_val_2.fq.gz
```

This produces a new FastQ file for Read 2 (ending in `.mispriming_corrected.fq.gz`), and prints a mispriming correction report such as this one:

```
Pre-processing report for bait file VkJk_test_J-region_val_2.fq.gz (Read 2)
========================================================
Total number of sequences analysed:     9899
Cases detected as mispriming in total: 520 (5.3%)

Sequence detected       Replaced with   Count
ATTTCCACCTTGGTGCCTCCACCGAACGTC  ATTTCCAGCTTGGTGCCTCCACCGAACGTC  13
ATTTCCAGCTTGGTCACAGC    TCCAGCTTGGTCCCAGC       0
ATTTCCAGCTTGGTCCCAGC    TCCAGCTTGGTCCCAGC       8
ATTTCCAGCTTGGTCCCCCCT   TTTCCAGCTTGGTCCCCCCT    9
ATTTCCAGCTTGGTGCCAGC    TCCAGCTTGGTCCCAGC       2
ATTTCCAGCTTGGTGCCCCCT   TTTCCAGCTTGGTCCCCCCT    7
ATTTCCAGCTTGGTGCCTCCT   TTTCCAGCTTGGTCCCCCCT    132
TCCAGCTTGGTCCCAGCCG     TTATTTCCAACTTTGTCCCCGAGCCG      27
TCCAGCTTGGTCCCAGCTCC    TTTCCAGCTTGGTCCCCCCTCC  36
TCCAGCTTGGTCCCCCCTCC    TTTCCAGCTTGGTCCCCCCTCC  111
TCCAGCTTGGTCCCCGAGCCG   TTATTTCCAACTTTGTCCCCGAGCCG      6
TCCAGCTTGGTGCCTCCA      ATTTCCAGCTTGGTGCCTCCA   22
TTATTTCCAACTTTGTCCCAGCACCGAACGTG        TCCAGCTTGGTCCCAGCACCGAACGTG     0
TTATTTCCAACTTTGTCCCCCCT TTTCCAGCTTGGTCCCCCCT    5
TTATTTCCAACTTTGTCCCCGAACCGAACGTG        TCCAGCTTGGTCCCAGCACCGAACGTG     1
TTATTTCCAACTTTGTCCCCGAACGTC     ATTTCCAGCTTGGTGCCTCCACCGAACGTC  7
TTATTTCCAACTTTGTCCCCGAACGTG     TCCAGCTTGGTCCCAGCACCGAACGTG     1
TTATTTCCAACTTTGTCCCCGCACCGAACGTG        TCCAGCTTGGTCCCAGCACCGAACGTG     0
TTTCCAGCTTGGTCCAAGCACC  TCCAGCTTGGTCCCAGCACC    0
TTTCCAGCTTGGTCCCAGCACC  TCCAGCTTGGTCCCAGCACC    67
TTTCCAGCTTGGTCCCAGCC    TTATTTCCAACTTTGTCCCCGAGCC       2
TTTCCAGCTTGGTCCCCGAGCC  TTATTTCCAACTTTGTCCCCGAGCC       26
TTTCCAGCTTGGTGCCAGCACC  TCCAGCTTGGTCCCAGCACC    0
TTTCCAGCTTGGTGCCTCCA    ATTTCCAGCTTGGTGCCTCCA   38
```


### Step 3: Running BabrahamLinkON_VkJk

Once the mispriming correction has been applied we are ready to run the J-bait detection with this command:

```
BabrahamLinkON_VkJk-Seq VkJk_test_V-region_val_1.fq.gz VkJk_test_J-region_val_2.mispriming_corrected.fq.gz
```

#### This process:

- identifies the J-bait at the start of Read 2 (Jk1, Jk2, Jk4 or Jk5), and extracts 20bp downstream of the bait
- identifies the anchor sequence in Read 1 (either GACTCGT or CTGCTCCT). Sequences without anchor sequence are discarded
- identifies a 6N barcode + either G or C from the anchor sequence in Read 1, and write this out together with the 20bp downstream of the Jk-bait, into a text file (e.g. `.Jk1.bait_side.txt`) like this: 

```
@HWI-1KL136:275:C4HL3ACXX:1:1101:7055:1987 4:N:0:CAGATC
AGATAGG:CGTCGGGTGGGGAACGATGA
+
JJIIIIJFHIJIIFHGFFFF
```
This file will be required for deduplication of V-gene hits in Step 5.

- Will extract the V-gene sequence from Read 1, starting from base 15 (V-reads look like this:
`NNNNNN AAAAAAA(A) T VVVVVVVVVVVV...`,
whereby `N` means random nucleotide (barcode), `A` is the anchor sequence (see above), `T` is a T to bind to A-tailed fragments, and `V` is true V-gene sequence we are interested in. The V-gene sequences are output in FastQ format, e.g. `Jk1.prey_side.fastq`:

```
@HWI-1KL136:275:C4HL3ACXX:1:1101:7055:1987 1:N:0:CAGATC
AGAGGACAAATTGTTCTCACCCAGTCTCCAGCAATCATGTCTGCATCTCTAGGGGAACGGGTCACCATGACCTGCACTGCCAGCTC
+
JIJJJJJJJJJJJHIJJJIGIIJIIIIJJJJJIIJJJJJIJJJJIIJJJJIGGGJHGFFFA=>BDDDCDDDDDDDDDDDDDDDDDD
```

### Step 4: Aligning the V-end to the genome
V-gene reads were aligned very stringently to the mouse genome (build NCBIM37) using Bowtie. Multimapping reads were discarded. Here is an example command:

```
bowtie -t -p 6 --chunkmbs 2048 -S --strata -m 1 --best NCBIM37/Mus_musculus.NCBIM37 VkJk_test_V-region_val_1.Jk1.prey_side.fastq | samtools view -bS - > VkJk_test_V-region_val_1.Jk1.prey_side.bam
```

### Step 5: Deduplicating the output results

This script is simply run by typing:

```
./deduplicate_VkJk-Seq
```

in the folder containing all the results. It automatically picks up files called `*bam` and `*bait_side.txt`, so the above steps 1-4 need to have been run for every single Jk bait file in the folder.

This script:
* stores all starting positions of V-side alignments
* uses 20 bp downstream of the J-region bait fragment
* uses the 7bp barcode sequence (see Step 3 for more details)
    
to deduplicate V-alignments. The J-read sequences in the `*bait_side.txt` files look like this:
    
> TAGACTG:CGTCCACGGGAATGTGTAAA (7N Barcode : 20bp J downstream sequence)

The output of the deduplication procedure is written to a new BAM file ending in `.unique_V.bam`.

