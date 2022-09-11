# step04_featurecounts

利用基因注释文件，获取每个基因比对的序列数目（read counts）。

首先，在项目目录内创建以下目录

```shell
# In project directory
mkdir -p results/featurecounts
```



## 执行命令

```shell
# In project directory
## M_L3_rep1
featureCounts -t exon -g gene_id -T 16 -a ~/reference/fly/annotation/dmel-all-r6.36.gtf -o results/featurecounts/M_L3_rep1_counts.txt results/align/M_L3_rep1_align/M_L3_rep1.sorted.bam
## M_L3_rep2
featureCounts -t exon -g gene_id -T 16 -a ~/reference/fly/annotation/dmel-all-r6.36.gtf -o results/featurecounts/M_L3_rep2_counts.txt results/align/M_L3_rep2_align/M_L3_rep2.sorted.bam
## M_L3_rep3
featureCounts -t exon -g gene_id -T 16 -a ~/reference/fly/annotation/dmel-all-r6.36.gtf -o results/featurecounts/M_L3_rep3_counts.txt results/align/M_L3_rep3_align/M_L3_rep3.sorted.bam
## M_WP_rep1
featureCounts -t exon -g gene_id -T 16 -a ~/reference/fly/annotation/dmel-all-r6.36.gtf -o results/featurecounts/M_WP_rep1_counts.txt results/align/M_WP_rep1_align/M_WP_rep1.sorted.bam
## M_WP_rep2
featureCounts -t exon -g gene_id -T 16 -a ~/reference/fly/annotation/dmel-all-r6.36.gtf -o results/featurecounts/M_WP_rep2_counts.txt results/align/M_WP_rep2_align/M_WP_rep2.sorted.bam
## M_WP_rep3
featureCounts -t exon -g gene_id -T 16 -a ~/reference/fly/annotation/dmel-all-r6.36.gtf -o results/featurecounts/M_WP_rep3_counts.txt results/align/M_WP_rep3_align/M_WP_rep3.sorted.bam
```

### 参数

`-T`: 使用的线程数目

`-a`: 基因注释文件

`-o`: 输出文件

## 输出结果

```shell
$ ls -lh featurecounts/
total 39M
-rw-rw-r-- 1 bioinfo bioinfo 4.3M Sep 11 16:25 M_BP_rep1_counts.txt
-rw-rw-r-- 1 bioinfo bioinfo  427 Sep 11 16:25 M_BP_rep1_counts.txt.summary
-rw-rw-r-- 1 bioinfo bioinfo 4.3M Sep 11 16:26 M_BP_rep2_counts.txt
-rw-rw-r-- 1 bioinfo bioinfo  427 Sep 11 16:26 M_BP_rep2_counts.txt.summary
-rw-rw-r-- 1 bioinfo bioinfo 4.3M Sep 11 16:28 M_BP_rep3_counts.txt
-rw-rw-r-- 1 bioinfo bioinfo  427 Sep 11 16:28 M_BP_rep3_counts.txt.summary
-rw-rw-r-- 1 bioinfo bioinfo 4.3M Sep 11 16:30 M_L3_rep1_counts.txt
-rw-rw-r-- 1 bioinfo bioinfo  427 Sep 11 16:30 M_L3_rep1_counts.txt.summary
-rw-rw-r-- 1 bioinfo bioinfo 4.3M Sep 11 16:32 M_L3_rep2_counts.txt
-rw-rw-r-- 1 bioinfo bioinfo  427 Sep 11 16:32 M_L3_rep2_counts.txt.summary
-rw-rw-r-- 1 bioinfo bioinfo 4.3M Sep 11 16:34 M_L3_rep3_counts.txt
-rw-rw-r-- 1 bioinfo bioinfo  427 Sep 11 16:34 M_L3_rep3_counts.txt.summary
-rw-rw-r-- 1 bioinfo bioinfo 4.3M Sep 11 16:37 M_WP_rep1_counts.txt
-rw-rw-r-- 1 bioinfo bioinfo  428 Sep 11 16:37 M_WP_rep1_counts.txt.summary
-rw-rw-r-- 1 bioinfo bioinfo 4.3M Sep 11 16:39 M_WP_rep2_counts.txt
-rw-rw-r-- 1 bioinfo bioinfo  427 Sep 11 16:39 M_WP_rep2_counts.txt.summary
-rw-rw-r-- 1 bioinfo bioinfo 4.3M Sep 11 16:40 M_WP_rep3_counts.txt
-rw-rw-r-- 1 bioinfo bioinfo  427 Sep 11 16:40 M_WP_rep3_counts.txt.summary
```

