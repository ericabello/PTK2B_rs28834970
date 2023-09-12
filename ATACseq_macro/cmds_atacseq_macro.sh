REFERENCE=/lustre/scratch117/sciops/team233/Erica/reference/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
DATA_DIR=/lustre/scratch117/sciops/team233/Erica/ATAC_seq_PTK2B_2019
SAMPLE=${2}
SAMPLENUMBER=${SAMPLE}
FRACTION=${3}
#R1=${DATA_DIR}/00_fastq/${SAMPLE}_R1_001.fastq.gz
#R2=${DATA_DIR}/00_fastq/${SAMPLE}_R2_001.fastq.gz
#R1TRIMMED=${DATA_DIR}/00_fastq/${SAMPLE}_R1_001_val_1.fq.gz
#R2TRIMMED=${DATA_DIR}/00_fastq/${SAMPLE}_R2_001_val_2.fq.gz
#FILE_STEM=output_example
#if [[ $4 == "ecoli" ]]
#then
#  SAMPLE="${SAMPLE}.ecolik12"
#fi
#if [[ $3 == "trimmed" ]]
#then
#  SAMPLE="${SAMPLE}.trimmed"
#if
case $1 in
  cramtobam)
	echo "converting ${SAMPLE} cram into bam file"
	samtools view -b ${DATA_DIR}/crams/${SAMPLE}.cram > ${DATA_DIR}/bams/${SAMPLE}.bam &
	;;
  unmapped)
    	echo "Removing unmapped and low quality reads from sample ${SAMPLE} and sorting"
  	samtools view -@ 4 -b -h -f 2 -F 4 -F 256 -F 2048 -q 30 ${DATA_DIR}/bams/${SAMPLE}.bam | samtools sort -o ${DATA_DIR}/bams/${SAMPLE}_sorted.bam &
	;;
  removeDuplicates)
	echo "removing duplicates from ${SAMPLE}"
	picard MarkDuplicates I=${DATA_DIR}/bams/${SAMPLE}_sorted.bam O=${DATA_DIR}/bams/${SAMPLE}_sorted_nodupl.bam REMOVE_DUPLICATES=true METRICS_FILE=${SAMPLE}_dupl.txt &
  	;;
  inserts)
    	echo "Checking size of inserts for sample ${SAMPLE}"
    	picard CollectInsertSizeMetrics I=${DATA_DIR}/bams/${SAMPLE}_sorted_nodupl.bam O=${DATA_DIR}/inserts/${SAMPLE}_insert_size.txt H=${DATA_DIR}/inserts/${SAMPLE}_insert_size_hist.pdf M=0.5 &
	;;
  readsinbam)
    	echo "Determining number of reads in BAM file for sample ${SAMPLE}"
    	READS=$(samtools flagstat ${DATA_DIR}/bams/${SAMPLE}_sorted_nodupl.bam)
    	echo -e "${SAMPLE}\t${READS}" >>read_counts_atac_macro_ptk2b.tsv
    	;;
  index)
	echo "Indexing sorted bam for sample ${SAMPLE}"
	samtools index -@ 2 ${DATA_DIR}/bams/${SAMPLE}_sorted_nodupl.bam &
	;;
  bigwig)
    	echo "Making bigwig for sample ${SAMPLE}"
    	bamCoverage --bam ${DATA_DIR}/bams/${SAMPLE}_sorted_nodupl.bam -o ${DATA_DIR}/bw/${SAMPLE}.bw --binSize 1 --normalizeUsing CPM --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX &
    	;;
  subsample_to_C11_2)
    	echo "Subsampling bam to ${FRACTION}"
    	samtools view -b -s ${FRACTION} ${DATA_DIR}/bams/${SAMPLE}_sorted_nodupl.bam > ${DATA_DIR}/bams/${SAMPLE}.subsampled.bam &
    	;;
  bamtobed)
	echo "transforming bam to bed for sample ${SAMPLE}"
	bedtools bamtobed -i ${DATA_DIR}/bams/${SAMPLE}.subsampled.bam > ${DATA_DIR}/beds/${SAMPLE}.subsampled.bed &
	;;
  #clean_bed)
	#echo " Cleaning bed for sample ${SAMPLE}"
	#cat ${DATA_DIR}/beds/${SAMPLE}.subsampled.bed | awk '{print $1 "\t" $2 "\t" $3 "\t" $6}' > ${DATA_DIR}/beds/${SAMPLE}.subsampled.clean.bed
#	;;
  macs2)
	echo "Running macs2 on sample ${SAMPLE}"
	macs2 callpeak -B -t ${DATA_DIR}/beds/${SAMPLE}.subsampled.bed --name ${SAMPLE}_macs2 --outdir ${DATA_DIR}/macs2 -g hs -f BED -q 0.01 --nomodel --keep-dup all --shift -100 --extsize 200 --SPMR
	;;
  mergebeds)
	shift
	shift
	BEDSTOMERGE=("$@")
	OUTPUTFILE=${DATA_DIR}/counts/ATAC_macro_narrowpeaks_allMerged.bed
	echo "Merging the following beds: ${BEDSTOMERGE[@]}"
	echo "Saving to ${OUTPUTFILE}"
	cat ${BEDSTOMERGE[@]} | sort -k1,1 -k2,2n | mergeBed -i stdin >${OUTPUTFILE}
    	;;
  sort_mergedbed)
	bedtools sort -faidx /lustre/scratch117/sciops/team233/Erica/reference/chromosome.sizes -i ${DATA_DIR}/counts/ATAC_macro_narrowpeaks_allMerged.bed > ${DATA_DIR}/counts/ATAC_macro_narrowpeaks_allMerged_sorted.bed &
	;;
  generate_counts)
	echo "Getting counts under peaks for sample ${SAMPLE}"
	coverageBed -a ${DATA_DIR}/counts/ATAC_macro_narrowpeaks_allMerged_sorted.bed -b ${DATA_DIR}/bams/${SAMPLE}_sorted_nodupl.bam -g /lustre/scratch117/sciops/team233/Erica/reference/chromosome.sizes -sorted -counts > ${DATA_DIR}/counts/${SAMPLE}_ATAC_macro_PTK2B.counts &
	;;
   *)
    	echo "Unknown argument for sample ${SAMPLE}"
	;;
esac
