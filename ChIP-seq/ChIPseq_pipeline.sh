#!/usr/bin/bash

help() {
		echo ""
		echo -e "Usage: \n\t bash $0 [options] -d <project_dir> -o <output_dir> -i <data_dir> --ref <reference_genome> "
		echo ""
		echo -e "  -d \n\t\t Project directory."
		echo -e "  -i \n\t\t Data directory for raw and trimmed fastqs. All raw fastqs should be put in ${data_dir}/fastq."
		echo -e "  -o \n\t\t Output directory for results."
		echo -e "  --ref \n\t\t Reference for bowtie2 mapping."
		echo ""
		echo "Options:"
		echo -e "  -t, --thread \n\t\t Numbers of threads, default: 1"
		echo -e "  --pair \n\t\t Flag to turn on pair-end mode."
		echo -e "  -h, --help \n\t\t show help message."
		echo ""
}


if [ $# -eq 0 ]; then
        echo "Run 'bash $0 -h' to see more information."
        exit 0
fi

# NOTE: This requires GNU getopt.  On Mac OS X and FreeBSD, you have to install this
# separately; see below.
TEMP=`getopt -o d:o:t:i:h --long help,ref:,threads:,pair \
             -n 'ChIPseq_pipeline.sh' -- "$@"`

if [ $? != 0 ] ; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

THREADS=1
PAIR=false
FQ1=
FQ2=
# Illumina universial adapters
#FORADAPT=AGATCGGAAGAG
#READAPT=AGATCGGAAGAG

while true; do
  case "$1" in
				-d ) PROJECTDIR="$2"; shift 2;;
				-o ) OUTDIR="$2"; shift 2;;
				-i ) DATADIR="$2"; shift 2;;
				-t | --threads ) THREADS="$2"; shift 2;;
				--pair ) PAIR=true; shift ;;
				--ref ) REFERENCE="$2"; shift 2;;
				-h | --help)
                        help
                        exit 0;;
                -- ) shift; break ;;
                * ) echo "Invalid option: ${optionName}" ;
				   echo "Try 'bash `basename $0` -h' for more information" ; 
				   break ;;
        esac
done

if [[ ! -d ${PROJECTDIR}/data ]] || [[ ! -d ${PROJECTDIR}/results ]]; then
	echo " The 'data' and 'results' folders must be in your project directory"
	exit 1
fi

if [[ ! -d ${OUTDIR} ]]; then
	echo " The output directory needed to be specified"
	exit 1
fi

if [[ ! -d ${DATADIR} ]]; then
	echo " The data directory needed to be specified"
	exit 1
fi

if [[ $REFERENCE = "" ]]; then
	echo " The reference genome is needed "
	exit 1
fi

## mkdir
if [[ ! -d ${DATADIR}/clean ]]; then mkdir -p ${DATADIR}/clean; fi
if [[ ! -d ${OUTDIR}/QC ]]; then mkdir -p ${OUTDIR}/QC; fi
FASTQLOC=${DATADIR}/fastq
CLEANLOC=${DATADIR}/clean
if [[ ! -d ${OUTDIR}/align ]]; then mkdir -p ${OUTDIR}/align ; fi
ALIGNLOC=${OUTDIR}/align
## Printing 
echo ""
echo " Project directory: ${PROJECTDIR} "
echo " Fastq directory: ${FASTQLOC}"
echo " Reference path: ${REFERENCE}"
echo " Numbers of threads: ${THREADS} "
echo ""
## list all fastq files
FQs=(`ls ${FASTQLOC}/*fastq.gz`)

## Processing with for loop
if $PAIR; then
	## Paired-end mode
	for ((j=0; j<=${#FQs[@]}-1; j+=2)); do
		FQ1=`basename ${FQs[$j]}`
		FQ2=`basename ${FQs[$j+1]}`
		echo "Processing " ${FQ1}" "${FQ2}
		## Extract the common prefix of fastq1 and fastq2
		PREFIX=`{ echo "$FQ1"; echo "$FQ2";  } | sed -e 'N;s/^\(.*\).*\n\1.*$/\1\n\1/;D'`

		## trimming
		echo ""
		echo ">Step01 Trimming & QC"
		echo ""
		trim_galore --nextseq 30 \
		--phred33 \
		--gzip \
		-o ${CLEANLOC}/ \
		--cores 8 \
		--basename ${PREFIX} \
		--fastqc_args "-o ${OUTDIR}/QC -t 8" \
		--paired \
		${FASTQLOC}/${FQ1} ${FASTQLOC}/${FQ2}

		## Mapping -- ../align/{PREFIX}/
		echo ""
		echo ">Step02 Mapping" 
		echo ""
		if [[ ! -d ${ALIGNLOC}/${PREFIX}_align ]]; then mkdir -p ${ALIGNLOC}/${PREFIX}_align; fi

		bowtie2 -x ${REFERENCE} -1 ${CLEANLOC}/${PREFIX}_val_1.fq.gz -2 ${CLEANLOC}/${PREFIX}_val_2.fq.gz -p ${THREADS} -S ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sam 2>${ALIGNLOC}/${PREFIX}_align/${PREFIX}_Align.log

		## samtools -- ../align/{PREFIX}/
		echo ""
		echo ">Step03 SAMTOOLS" 
		echo ""
		samtools view -@ ${THREADS} -q 30 -f 2 -hSb ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sam |samtools sort - -@ 32 -o ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		samtools index ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam

		## bamCoverage -- ../align/{PREFIX}/
		echo ""
		echo ">Step04 bam2bw" 
		echo ""
		bamCoverage -b ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam -o ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.bw --normalizeUsing CPM -bs 1 -p ${THREADS}

	done
else
	## Single-end mode
	for ((j=0; j<=${#FQs[@]}-1; j++)); do
		FQ=`basename ${FQs[$j]}`
		echo "Processing " ${FQ}
		PREFIX=`basename ${FQ} .fastq.gz`

		## trimming
		echo ""
		echo ">Step01 Trimming & QC"
		echo ""
		trim_galore --nextseq 30 \
		--phred33 \
		--gzip \
		-o ${CLEANLOC}/ \
		--cores 8 \
		--basename ${PREFIX} \
		--fastqc_args "-o ${OUTDIR}/QC -t 8" \
		${FASTQLOC}/${FQ} 

		## Mapping -- ../align/{PREFIX}/
		echo ""
		echo ">Step02 Mapping" 
		echo ""
		if [[ ! -d ${ALIGNLOC}/${PREFIX}_align ]]; then mkdir -p ${ALIGNLOC}/${PREFIX}_align; fi
		
		bowtie2 -x ${REFERENCE} -U ${CLEANLOC}/${PREFIX}_trimmed.fq.gz -p ${THREADS} -S ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sam 2>${ALIGNLOC}/${PREFIX}_align/${PREFIX}_Align.log

		## samtools -- ../align/{PREFIX}/
		echo ""
		echo ">Step03 SAMTOOLS" 
		echo ""
		samtools view -@ ${THREADS} -q 30 -F 4 -hSb ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sam |samtools sort - -@ 32 -o ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam
		samtools index ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam

		## bamCoverage -- ../align/{PREFIX}/
		echo ""
		echo ">Step04 bam2bw" 
		echo ""
		bamCoverage -b ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.sorted.bam -o ${ALIGNLOC}/${PREFIX}_align/${PREFIX}.bw --normalizeUsing CPM -bs 1 -p ${THREADS}
	done
fi

# post multiqc
## trimmed reads qc 
multiqc -o ${OUTDIR}/QC/  ${OUTDIR}/QC/*zip
## bowtie2 alignment reports
multiqc -o ${OUTDIR}/align/  ${OUTDIR}/align/*/*Align.log

# peak calling by macs2
echo ""
echo ">Step05 Peak calling" 
echo ""
if [[ ! -d ${OUTDIR}/peak ]]; then mkdir -p ${OUTDIR}/peak; fi
PEAKLOC=${OUTDIR}/peak
## L3 rep1
macs2 callpeak -t ${ALIGNLOC}/L3_rep1_IP_align/L3_rep1_IP.sorted.bam -c ${ALIGNLOC}/L3_rep1_Input_align/L3_rep1_Input.sorted.bam -f BAM -g dm -n L3_rep1 --outdir ${PEAKLOC}/ 2>${PEAKLOC}/L3_rep1_macs2.log
## L3 rep2
macs2 callpeak -t ${ALIGNLOC}/L3_rep2_IP_align/L3_rep2_IP.sorted.bam -c ${ALIGNLOC}/L3_rep2_Input_align/L3_rep2_Input.sorted.bam -f BAM -g dm -n L3_rep2 --outdir ${PEAKLOC}/ 2>${PEAKLOC}/L3_rep2_macs2.log
## WP rep1
macs2 callpeak -t ${ALIGNLOC}/WP_rep1_IP_align/WP_rep1_IP.sorted.bam -c ${ALIGNLOC}/WP_rep1_Input_align/WP_rep1_Input.sorted.bam -f BAM -g dm -n WP_rep1 --outdir ${PEAKLOC}/ 2>${PEAKLOC}/WP_rep1_macs2.log
## WP rep2
macs2 callpeak -t ${ALIGNLOC}/WP_rep2_IP_align/WP_rep2_IP.sorted.bam -c ${ALIGNLOC}/WP_rep2_Input_align/WP_rep2_Input.sorted.bam -f BAM -g dm -n WP_rep2 --outdir ${PEAKLOC}/ 2>${PEAKLOC}/WP_rep2_macs2.log

# extract overlappping peaks between replicates
## L3
bedtools intersect -wo -f 0.4 -r -a ${PEAKLOC}/L3_rep1_peaks.narrowPeak -b ${PEAKLOC}/L3_rep2_peaks.narrowPeak > ${PEAKLOC}/L3_peaks_final.bed
## WP
bedtools intersect -wo -f 0.4 -r -a ${PEAKLOC}/WP_rep1_peaks.narrowPeak -b ${PEAKLOC}/WP_rep2_peaks.narrowPeak > ${PEAKLOC}/WP_peaks_final.bed
# visualization
echo ""
echo ">Step06 Visualization" 
echo ""
if [[ ! -d ${OUTDIR}/vis ]]; then mkdir -p ${OUTDIR}/vis; fi
VISLOC=${OUTDIR}/vis
## L3
computeMatrix reference-point -S ${ALIGNLOC}/L3_rep1_IP_align/L3_rep1_IP.bw -R ${PEAKLOC}/L3_peaks_final.bed --referencePoint TSS -a 1000 -b 1000 -out ${VISLOC}/L3_K4_rep1.tab.gz -p 16
computeMatrix reference-point -S ${ALIGNLOC}/L3_rep2_IP_align/L3_rep2_IP.bw -R ${PEAKLOC}/L3_peaks_final.bed --referencePoint TSS -a 1000 -b 1000 -out ${VISLOC}/L3_K4_rep2.tab.gz -p 16
plotHeatmap -m ${VISLOC}/L3_K4_rep1.tab.gz -out ${VISLOC}/L3_K4_rep1_Heatmap.png --whatToShow "plot, heatmap and colorbar"
plotHeatmap -m ${VISLOC}/L3_K4_rep2.tab.gz -out ${VISLOC}/L3_K4_rep2_Heatmap.png --whatToShow "plot, heatmap and colorbar"
## WP
computeMatrix reference-point -S ${ALIGNLOC}/WP_rep1_IP_align/WP_rep1_IP.bw -R ${PEAKLOC}/WP_peaks_final.bed --referencePoint TSS -a 1000 -b 1000 -out ${VISLOC}/WP_K4_rep1.tab.gz -p 16
computeMatrix reference-point -S ${ALIGNLOC}/WP_rep2_IP_align/WP_rep2_IP.bw -R ${PEAKLOC}/WP_peaks_final.bed --referencePoint TSS -a 1000 -b 1000 -out ${VISLOC}/WP_K4_rep2.tab.gz -p 16
plotHeatmap -m ${VISLOC}/WP_K4_rep1.tab.gz -out ${VISLOC}/WP_K4_rep1_Heatmap.png --whatToShow "plot, heatmap and colorbar"
plotHeatmap -m ${VISLOC}/WP_K4_rep2.tab.gz -out ${VISLOC}/WP_K4_rep2_Heatmap.png --whatToShow "plot, heatmap and colorbar"