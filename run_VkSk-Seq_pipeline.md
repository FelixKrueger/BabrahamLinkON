# Running the VkJk-Seq pipeline

### Step 1: Quality/adapter trimming and 5' end clipping

In the current implementation of the VkJk-Seq pipeline we assume that the basecall qualities of the first 4 base pairs of Read 2 (the J-reads) are quite poor (

```
trim_galore --clip_r2 4 --paired VkJk_test_V-region.fastq.gz VkJk_test_J-region.fastq.gz
```



### Step 2: Mispriming correction



### Step 3: Running BabrahamLinkON_VkJK


### Step 4: Aligning the V-end to the genome

### Step 5: Deduplicating the output results


