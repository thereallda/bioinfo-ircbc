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

**命令中所有的`$`只是为了提示该命令是在终端中输入，实际输入时，无需包括`$`符号。**

> `ll`为`ls -l`的别名（alias）

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
$ tree reference/
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
$ tree src /
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



**conda 相关操作命令**

```{shell}
## 查看已安装命令
conda list 

## 安装软件到指定环境
conda install -n envname pkgname
## or
conda install --prefix=/path/to/install pkgname

## 安装指定channel的指定版本软件
conda install -c conda-forge R r=4.0

## 卸载包
conda uninstall pkgname

## 更新包
conda update pkgname

## 显示所有的虚拟环境                   
conda info --envs  
## or
conda env list

## 搜索包
conda search pkgname

## 创建环境
conda create --name envname
## 复制环境envname1到envname2
conda create --name envname2 --clone envname1
## 删除环境
conda remove --name envname --all
```



## cheat-sheet

| Command                       | Description                                                  |
| :---------------------------- | :----------------------------------------------------------- |
| `ls`                          | List all files and directories in the current working directory |
| `ls -lah`                     | List all files and directories including, hidden files and other information like permissions, size, and owner in human-readable format. |
| `cd`                          | Change the directory to the home directory (`~`)             |
| `cd ..`                       | Change the directory to one level up                         |
| `cat filename`                | Display the content of the file                              |
| `cat file1 file2 > file3`     | Combine two files named file1 and file2 and store the output in a new file file3 |
| `tail filename`               | Display the last 10 lines of a file                          |
| `head filename`               | Display the first 10 lines of a file                         |
| `mv oldfile newfile`          | Rename a file                                                |
| `rm filename`                 | Delete a file                                                |
| `mkdir dirname`               | Create a directory                                           |
| `mkdir -p parentDir/childDir` | Make parent directories if needed                            |
| `rm -rf dirname`              | Remove a directory                                           |
| `history`                     | Print a history list of all commands                         |
| `clear`                       | Clear the terminal                                           |
| `hostnamectl`                 | Get system information including, operating system, kernel, and release version |
| `date`                        | Display the current system date and time                     |
| `hostname`                    | Display the hostname of the system                           |
| `ifconfig`                    | Display the IP and Mac Address of the system                 |
| `w`                           | Display currently logged in users in the system              |
| `free -m`                     | Display free and used memory in the system                   |
| `top`                         | Display all running processes                                |



> For more please check: 
>
> https://www.pcwdld.com/linux-commands-cheat-sheet