# step03_samtools

使用samtools：

- 将比对结果的sam转换为bam，以二进制方式保存以节省空间
- sort bam
- index bam

## 执行命令

```shell
# In project directory
## use for loop to convet sam to bam, sort and index bam
for id in `ls results/align/`; do
	prefix=`basename _align`
	samtools view -@ 16 -q 30 -F 4 -hSb results/align/${prefix}_align/Aligned.out.sam -o results/align/${prefix}_align/${prefix}.bam
samtools sort -@ 16 -o results/align/${prefix}_align/${prefix}.sorted.bam
samtools index results/align/${prefix}_align/${prefix}.sorted.bam
done
```



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



