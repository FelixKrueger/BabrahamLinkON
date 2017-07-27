# Running the VkJk-Seq pipeline

The data associated with the VkJk publication will soon be available from GEO under the accession number GSE101606. Here we provide two test files of the WT_pre-B_VkJk-seq sample (the first 10,000 sequences) to enable reproduction of how the pipeline was run. The files are called: 

 * Read 1: `VkJk_test_V-region.fastq.gz`
 * Read 2: `VkJk_test_J-region.fastq.gz`
 
 
### Step 1: Quality/adapter trimming and 5' end clipping

In the current implementation of the VkJk-Seq pipeline we assume that the basecall qualities of the first 4 base pairs of Read 2 (the J-reads) are quite poor. We thus clip R2 by 4 bp from the 5' end. In addition, Trim Galore detects and removes read-through adapter contamination and poor quality basecalls from the 5' ends. The trimmin command is:

```
trim_galore --clip_r2 4 --paired VkJk_test_V-region.fastq.gz VkJk_test_J-region.fastq.gz
```
The output files are called: `VkJk_test_V-region_val_1.fq.gz` and `VkJk_test_J-region_val_2.fq.gz`.


### Step 2: Mispriming correction



### Step 3: Running BabrahamLinkON_VkJK


### Step 4: Aligning the V-end to the genome

### Step 5: Deduplicating the output results


