# step02_Mapping

将质控后的序列与参考基因组进行比对。



首先，在项目目录内创建以下目录

```shell
# In project directory
mkdir -p results/align
```

为了方便结果的输出，对每一个样本创建一个比对结果的输出目录。

```shell
# In project directory
# use for loop to create directories
for id in `ls data/fastq/`;do 
	prefix=`basename ${id} .fastq.gz`
	mkdir -p results/test/${prefix}_align
done
```

## 执行命令

使用STAR进行序列比对

```shell
# In project directory
## M_L3_rep1
STAR --genomeDir ~/reference/fly/index/star_index/dmel_r6.36/ --readFilesCommand zcat --readFilesIn data/clean/M_L3_rep1_trimmed.fq.gz --runThreadN 16 --outSAMtype SAM --outFileNamePrefix results/align/M_L3_rep1_align/
## M_L3_rep2
STAR --genomeDir ~/reference/fly/index/star_index/dmel_r6.36/ --readFilesCommand zcat --readFilesIn data/clean/M_L3_rep2_trimmed.fq.gz --runThreadN 16 --outSAMtype SAM --outFileNamePrefix results/align/M_L3_rep2_align/
## M_L3_rep3
STAR --genomeDir ~/reference/fly/index/star_index/dmel_r6.36/ --readFilesCommand zcat --readFilesIn data/clean/M_L3_rep3_trimmed.fq.gz --runThreadN 16 --outSAMtype SAM --outFileNamePrefix results/align/M_L3_rep3_align/
## M_WP_rep1
STAR --genomeDir ~/reference/fly/index/star_index/dmel_r6.36/ --readFilesCommand zcat --readFilesIn data/clean/M_WP_rep1_trimmed.fq.gz --runThreadN 16 --outSAMtype SAM --outFileNamePrefix results/align/M_WP_rep1_align/
## M_WP_rep2
STAR --genomeDir ~/reference/fly/index/star_index/dmel_r6.36/ --readFilesCommand zcat --readFilesIn data/clean/M_WP_rep2_trimmed.fq.gz --runThreadN 16 --outSAMtype SAM --outFileNamePrefix results/align/M_WP_rep2_align/
## M_WP_rep3
STAR --genomeDir ~/reference/fly/index/star_index/dmel_r6.36/ --readFilesCommand zcat --readFilesIn data/clean/M_WP_rep3_trimmed.fq.gz --runThreadN 16 --outSAMtype SAM --outFileNamePrefix results/align/M_WP_rep3_align/
```

### 参数

`--genomeDir `

`--readFilesCommand`

`--readFilesIn`

`--runThreadN`

`--outSAMtype`

`-outFileNamePrefix`

## 输出结果

```shell
$ ls -lh M_L3_rep1_align/
total 1.9G
-rw-rw-r-- 1 bioinfo bioinfo 1.7G Sep 11 16:30 Aligned.out.sam
-rw-rw-r-- 1 bioinfo bioinfo 2.0K Sep 11 16:30 Log.final.out
-rw-rw-r-- 1 bioinfo bioinfo  74K Sep 11 16:30 Log.out
-rw-rw-r-- 1 bioinfo bioinfo  246 Sep 11 16:30 Log.progress.out
-rw-rw-r-- 1 bioinfo bioinfo 187M Sep 11 16:30 M_L3_rep1.sorted.bam
-rw-rw-r-- 1 bioinfo bioinfo 161K Sep 11 16:30 M_L3_rep1.sorted.bam.bai
-rw-rw-r-- 1 bioinfo bioinfo 1.4M Sep 11 16:30 SJ.out.tab
```



