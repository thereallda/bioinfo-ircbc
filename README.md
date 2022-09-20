# bioinfo-ircbc

Code and documents for IRCBC bioinformatic class.

## log in

```shell
# In terminal
ssh bioinfo@192.168.221.161
```

在家目录 `/share/home/bioinfo` 中，包括以下目录

```shell
$ ll
total 0
drwxrwxr-x 5 bioinfo bioinfo 4096 Sep 19 20:26 ChIPseq_demo
drwxrwxr-x 4 bioinfo bioinfo 4096 Sep 11 14:28 data
drwxrwxr-x 3 bioinfo bioinfo 4096 Sep 11 14:25 reference
drwxrwxr-x 5 bioinfo bioinfo 4096 Sep 11 14:42 RNAseq_demo
drwxrwxr-x 3 bioinfo bioinfo 4096 Sep 11 13:26 softwares
drwxrwxr-x 2 bioinfo bioinfo 4096 Sep 19 21:46 src
```

`data`: 包含果蝇发育时期third-stage larvae 3（L3）和white pupa（WP）的RNA-seq和ChIP-seq的两组数据。

```shell
$ tree data/
data/
├── ChIP-seq_fly_dev
│   ├── L3_rep1_Input.fastq.gz
│   ├── L3_rep1_IP.fastq.gz
│   ├── L3_rep2_Input.fastq.gz
│   ├── L3_rep2_IP.fastq.gz
│   ├── WP_rep1_Input.fastq.gz
│   ├── WP_rep1_IP.fastq.gz
│   ├── WP_rep2_Input.fastq.gz
│   └── WP_rep2_IP.fastq.gz
└── RNA-seq_fly_dev
    ├── L3_rep1.fastq.gz
    ├── L3_rep2.fastq.gz
    ├── L3_rep3.fastq.gz
    ├── WP_rep1.fastq.gz
    ├── WP_rep2.fastq.gz
    └── WP_rep3.fastq.gz

2 directories, 14 files
```



`reference`: 果蝇基因注释文件，参考基因组和相应的index文件

```shell
reference/
└── fly
    ├── annotation
    │   └── dmel-all-r6.36.gtf
    ├── genome
    │   └── dmel-all-chromosome-r6.36.fasta
    └── index
        ├── bowtie2_index
        │   └── dm6.36
        │       ├── genome.1.bt2
        │       ├── genome.2.bt2
        │       ├── genome.3.bt2
        │       ├── genome.4.bt2
        │       ├── genome.rev.1.bt2
        │       └── genome.rev.2.bt2
        └── star_index
            └── dmel_r6.36
                ├── chrLength.txt
                ├── chrNameLength.txt
                ├── chrName.txt
                ├── chrStart.txt
                ├── exonGeTrInfo.tab
                ├── exonInfo.tab
                ├── geneInfo.tab
                ├── Genome
                ├── genomeParameters.txt
                ├── Log.out
                ├── nohup.out
                ├── SA
                ├── SAindex
                ├── sjdbInfo.txt
                ├── sjdbList.fromGTF.out.tab
                ├── sjdbList.out.tab
                └── transcriptInfo.tab

8 directories, 25 files
```

`src`: 脚本和相关源文件存放的目录

```shell
src
├── ChIPseq_pipeline.sh
├── ChIPseq_py3.yaml
├── contrast.csv
├── DEStream_demo.R
├── mergeCounts.R
├── metadata.csv
├── r_de.yaml
├── RNAseq_pipeline.sh
└── RNAseq_py3.yaml

0 directories, 9 files
```



`softwares`: 软件安装的目录



## conda environment

```shell
$ conda env list
# conda environments:
#
base                  *  /share/home/bioinfo/softwares/miniconda3
ChIPseq_py3              /share/home/bioinfo/softwares/miniconda3/envs/ChIPseq_py3
RNAseq_py3               /share/home/bioinfo/softwares/miniconda3/envs/RNAseq_py3
r_de                     /share/home/bioinfo/softwares/miniconda3/envs/r_de
```

`ChIPseq_py3`: ChIP-seq分析相关的环境

`RNAseq_py3`: RNA-seq分析相关的环境

`r_de`: 使用R进行差异基因表达分析的环境

