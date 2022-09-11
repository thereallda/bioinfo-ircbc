# step01_trimming_qc

在项目目录中，创建以下目录

```shell
mkdir -p data/clean # for trimmed fastq
mkdir -p results/QC # for fastq quality control report
```

使用 `trim_galore` 进行reads trimming and quality control report 

```shell
# In project directory
## M_L3_rep1
trim_galore -q 30 --gzip -o data/clean/ --cores 8 --basename M_L3_rep1 --fastqc_args "-o results/QC -t 8" data/fastq/M_L3_rep1.fastq.gz
## M_L3_rep2
trim_galore -q 30 --gzip -o data/clean/ --cores 8 --basename M_L3_rep2 --fastqc_args "-o results/QC -t 8" data/fastq/M_L3_rep2.fastq.gz
## M_L3_rep3
trim_galore -q 30 --gzip -o data/clean/ --cores 8 --basename M_L3_rep3 --fastqc_args "-o results/QC -t 8" data/fastq/M_L3_rep3.fastq.gz
## M_WP_rep1
trim_galore -q 30 --gzip -o data/clean/ --cores 8 --basename M_WP_rep1 --fastqc_args "-o results/QC -t 8" data/fastq/M_WP_rep1.fastq.gz
## M_WP_rep2
trim_galore -q 30 --gzip -o data/clean/ --cores 8 --basename M_WP_rep2 --fastqc_args "-o results/QC -t 8" data/fastq/M_WP_rep2.fastq.gz
## M_WP_rep3
trim_galore -q 30 --gzip -o data/clean/ --cores 8 --basename M_WP_rep3 --fastqc_args "-o results/QC -t 8" data/fastq/M_WP_rep3.fastq.gz
```

其中，

`-q <INT>`: 去除序列中低于输入值（`<INT>`）的碱基。这里我们去除质量低于30的碱基

`--gzip`: 将输出压缩为 GZIP文件

`-o`: 文件输出目录

`--cores`: 运行线程数，程序限制不得超过8。

`--basename <PREFERRED_NAME>`: 输出文件的前缀名。 如果是单端测序文件，输出文件将命名为 PREFERRED_NAME_trimmed.fq(.gz)；如果是双端测序文件，输出文件将命名为 PREFERRED_NAME_val_1.fq(.gz) and PREFERRED_NAME_val_2.fq(.gz) 

`--fastq_args`: trim galore调用fastqc时，所用的参数。 其中 `-O` 指定fastqc结果输出的路径，`-t` 指定fastqc所用的线程数。



执行后，在 `data/clean` 下产生相应的trimmed fastq文件；在 `results/QC` 下生成质控报告



进一步，我们可以用 `multiqc` 对质控报告进行汇总

```shell
# In project directory
multiqc -o results/QC/  results/QC/*zip
```

 `multiqc` 在`results/QC`目录下产生 `multiqc_report.html`，可以用本地电脑打开查看。

