REFERENCE=/lustre/scratch117/sciops/team233/Erica/reference/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
ECOLIREF=/lustre/scratch117/sciops/team233/Erica/reference/E.coliK12/ecoli
DATA_DIR=/lustre/scratch117/sciops/team233/Erica/CR_PTK2B_micro_110822/fastq_CR
SAMPLE=${2}
SAMPLENUMBER=${SAMPLE}
FRACTION=${5}
R1=${DATA_DIR}/00_fastq/${SAMPLE}_R1_001.fastq.gz
R2=${DATA_DIR}/00_fastq/${SAMPLE}_R2_001.fastq.gz
R1TRIMMED=${DATA_DIR}/00_fastq/${SAMPLE}_R1_001_val_1.fq.gz
R2TRIMMED=${DATA_DIR}/00_fastq/${SAMPLE}_R2_001_val_2.fq.gz
FILE_STEM=output_example
if [[ $4 == "ecoli" ]]
then
  SAMPLE="${SAMPLE}.ecolik12"
fi
if [[ $3 == "trimmed" ]]
then
  SAMPLE="${SAMPLE}.trimmed"
fi
case $1 in
  trim)
    echo "Trimming fastqs for sample ${SAMPLE}"
    trim_galore -q 20 -o ${DATA_DIR}/00_fastq/ --paired -a GATCGGAAGAGCACACGTCTGAACTCCAGTCA -a2 GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT ${R1} ${R2} &
    ;;
  align)
    echo "Aligning sample ${SAMPLE} with bowtie2"
    bowtie2 -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x ${REFERENCE} -1 ${R1} -2 ${R2} -S ${DATA_DIR}/sams/${SAMPLE}.sam &
    ;;
  aligntrimmed)
    echo "Aligning trimmed sample ${SAMPLE} with bowtie2"
    bowtie2 -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x ${REFERENCE} -1 ${R1TRIMMED} -2 ${R2TRIMMED} -S ${DATA_DIR}/sams/${SAMPLE}.trimmed.sam &
    ;;
  aligntrimmedecoli)
    echo "Aligning trimmed sample ${SAMPLE} with bowtie2 to e coli K12 genome"
    bowtie2 -q --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x ${ECOLIREF} -1 ${R1TRIMMED} -2 ${R2TRIMMED} -S ${DATA_DIR}/sams/${SAMPLE}.ecolik12.trimmed.sam &
    ;;
  bam)
    echo "Coverting sam to bam for sample ${SAMPLE}"
    samtools view -@ 2 -b -o ${DATA_DIR}/bams/${SAMPLE}.bam ${DATA_DIR}/sams/${SAMPLE}.sam &
    ;;
  sort)
    echo "Sorting bam for sample ${SAMPLE}"
    samtools sort -@ 2 -o ${DATA_DIR}/bams/${SAMPLE}.sorted.bam ${DATA_DIR}/bams/${SAMPLE}.bam &
    ;;
  index)
    echo "Indexing sorted bam for sample ${SAMPLE}"
    samtools index -@ 2 ${DATA_DIR}/bams/${SAMPLE}.sorted.bam &
    ;;
  inserts)
    echo "Checking size of inserts for sample ${SAMPLE}"
    bamPEFragmentSize -b ${DATA_DIR}/bams/${SAMPLE}.sorted.bam --maxFragmentLength 500 -o ${DATA_DIR}/inserts/${SAMPLE}_remapped_insert_size_hist.png --table ${DATA_DIR}/inserts/${SAMPLE}_remapped_insert_size.txt &
    ;;
  bigwig)
    echo "Making bigwig for sample ${SAMPLE}"
    bamCoverage --bam ${DATA_DIR}/bams/${SAMPLE}.sorted.bam -o ${DATA_DIR}/bw/${SAMPLE}.bw --binSize 1 --normalizeUsing CPM --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX &
    ;;
  unmapped)
    echo "Removing unmapped and low quality reads from sample ${SAMPLE}"
    samtools view -@ 4 -b -h -f 2 -F 4 -F 256 -F 2048 -q 30 -o ${DATA_DIR}/macs/${SAMPLE}_4macs.bam ${DATA_DIR}/bams/${SAMPLE}.sorted.bam &
    ;; 
  namesort)
    echo "Sorting bam by read name for sample ${SAMPLE}"
    samtools sort -n -o ${DATA_DIR}/seacr/${SAMPLE}.name.sorted.bam ${DATA_DIR}/macs/${SAMPLE}_4macs.bam &
    ;;
  bamtobed)
    echo "Converting bam to bed for sample ${SAMPLE}"
    bedtools bamtobed -bedpe -i ${DATA_DIR}/seacr/${SAMPLE}.name.sorted.bam > ${DATA_DIR}/seacr/${SAMPLE}.bed &
    # bedtools bamtobed -i ${DATA_DIR}/seacr/${SAMPLE}.name.sorted.bam > ${DATA_DIR}/seacr/${SAMPLE}.bed &
    ;;
  cleanbed)
    echo "Cleaning bed for sample ${SAMPLE}"
    awk '$1==$4 && $6-$2 < 1000 {print $0}' ${DATA_DIR}/seacr/${SAMPLE}.bed > ${DATA_DIR}/seacr/${SAMPLE}.clean.bed
    cut -f 1,2,6 ${DATA_DIR}/seacr/${SAMPLE}.clean.bed | sort -k1,1 -k2,2n -k3,3n > ${DATA_DIR}/seacr/${SAMPLE}.fragments.bed
    ;;
  addlength)
    echo "Calculating fragment length and adding to BED file for sample ${SAMPLE}"
    cat ${DATA_DIR}/seacr/${SAMPLE}.fragments.bed | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > ${DATA_DIR}/seacr/${SAMPLE}.fragments.length.bed
    ;;
  spikeinnorm)
    echo "Determining number of reads in BAM file for sample ${SAMPLE}"  
    READS=$(samtools view -c /lustre/scratch117/sciops/team233/Erica/CR_PTK2B_micro_110822/fastq_CR/macs/${SAMPLE}_4macs.bam)
    echo "Found ${READS} reads"
    SCALE=$(printf "%.0f" $(echo "scale=5; (1/($READS/1000000))*100000" | bc))
    echo "Normalising to spike-in DNA for sample ${SAMPLE}, with scale factor ${SCALE}"
    echo "Normalising to spike-in DNA for sample ${SAMPLE}"
    /lustre/scratch117/sciops/team233/Erica/CR_PTK2B_micro_110822/fastq_CR/spike_in_calibration.csh ${DATA_DIR}/seacr/${SAMPLE}.fragments.length.bed ${DATA_DIR}/seacr/${SAMPLE}.ecolik12.trimmed.fragments.length.bed ${SCALE} bg /lustre/scratch117/sciops/team233/Erica/reference/chromosome.sizes 1 1000 ${DATA_DIR}/seacr &
    #./spike_in_calibration.csh ${DATA_DIR}/seacr/${SAMPLE}.fragments.length.bed ${DATA_DIR}/seacr/${SAMPLENUMBER}.ecolik12.trimmed.fragments.length.bed 10000 bg /home/jon_kerry/CutAndRun/chromosome.sizes 1 1001 ${DATA_DIR}/seacr &
    ;;
  bgbigwig)
    echo "Make bigWig from bedGraph for sample ${SAMPLE}" 
    bedGraphToBigWig ${DATA_DIR}/seacr/${SAMPLE}.fragments.length.1-1000.bedgraph /lustre/scratch117/core/sciops_repository/references/Homo_sapiens/GRCh38_15/all/star/chrNameLength.txt ${DATA_DIR}/bw/${SAMPLE}.spikenorm.bw &
    ;;
  subsample_to_D3)
    echo "Subsampling bam to ${FRACTION}"
    samtools view -b -s ${FRACTION} ${DATA_DIR}/bams/${SAMPLE}.sorted.bam > ${DATA_DIR}/bams/${SAMPLE}.sorted.subsampled.bam &
    ;;
  unmapped_subsampled)
    echo "Removing unmapped and low quality reads from sample ${SAMPLE}"
    samtools view -@ 4 -b -h -f 2 -F 4 -F 256 -F 2048 -q 30 -o ${DATA_DIR}/macs/${SAMPLE}_4macs.subsampled.bam ${DATA_DIR}/bams/${SAMPLE}.sorted.subsampled.bam &
    ;; 
  namesort_subsampled)
    echo "Sorting bam by read name for sample ${SAMPLE}"
    samtools sort -n -o ${DATA_DIR}/seacr/${SAMPLE}.name.sorted.subsampled.bam ${DATA_DIR}/macs/${SAMPLE}_4macs.subsampled.bam &
    ;;
  bamtobed_subsampled)
    echo "Converting bam to bed for sample ${SAMPLE}"
    bedtools bamtobed -bedpe -i ${DATA_DIR}/seacr/${SAMPLE}.name.sorted.subsampled.bam > ${DATA_DIR}/seacr/${SAMPLE}.subsampled.bed &
    # bedtools bamtobed -i ${DATA_DIR}/seacr/${SAMPLE}.name.sorted.bam > ${DATA_DIR}/seacr/${SAMPLE}.bed &
    ;;
  cleanbed_subsampled)
    echo "Cleaning bed for sample ${SAMPLE}"
    awk '$1==$4 && $6-$2 < 1000 {print $0}' ${DATA_DIR}/seacr/${SAMPLE}.subsampled.bed > ${DATA_DIR}/seacr/${SAMPLE}.clean.subsampled.bed
    cut -f 1,2,6 ${DATA_DIR}/seacr/${SAMPLE}.clean.subsampled.bed | sort -k1,1 -k2,2n -k3,3n > ${DATA_DIR}/seacr/${SAMPLE}.fragments.subsampled.bed
    ;;
  addlength_subsampled)
    echo "Calculating fragment length and adding to BED file for sample ${SAMPLE}"
    cat ${DATA_DIR}/seacr/${SAMPLE}.fragments.subsampled.bed | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > ${DATA_DIR}/seacr/${SAMPLE}.fragments.length.subsampled.bed
    ;;
  spikeinnorm_subsampled)
    echo "Determining number of reads in BAM file for sample ${SAMPLE}"  
    READS=$(samtools view -c /lustre/scratch117/sciops/team233/Erica/CR_PTK2B_micro_110822/fastq_CR/macs/${SAMPLE}_4macs.subsampled.bam)
    echo "Found ${READS} reads"
    SCALE=$(printf "%.0f" $(echo "scale=5; (1/($READS/1000000))*100000" | bc))
    echo "Normalising to spike-in DNA for sample ${SAMPLE}, with scale factor ${SCALE}"
    echo "Normalising to spike-in DNA for sample ${SAMPLE}"
    /lustre/scratch117/sciops/team233/Erica/CR_PTK2B_micro_110822/fastq_CR/spike_in_calibration.csh ${DATA_DIR}/seacr/${SAMPLE}.fragments.length.subsampled.bed ${DATA_DIR}/seacr/${SAMPLE}.ecolik12.trimmed.fragments.length.bed ${SCALE} bg /lustre/scratch117/sciops/team233/Erica/reference/chromosome.sizes 1 1000 ${DATA_DIR}/seacr &
    ;;

  #genomecov)
   # echo "Creating bedgraph for sample ${SAMPLE}"
   # bedtools genomecov -bg -i ${DATA_DIR}/seacr/${SAMPLE}.fragments.bed -g /lustre/scratch117/sciops/team233/Erica/reference/chromosome.sizes > ${DATA_DIR}/seacr/${SAMPLE}.fragments.bedgraph &
    #;;
  seacr)
    echo "Running SEACR on sample ${SAMPLE}"
    bash /lustre/scratch117/sciops/team233/Erica/SEACR/SEACR_1.3.sh ${DATA_DIR}/seacr/${SAMPLE}.fragments.length.subsampled.1-1000.bedgraph ${DATA_DIR}/seacr/IgG_CR.fragments.length.1-1000.bedgraph norm stringent ${DATA_DIR}/seacr/${SAMPLE}.seacr.norm.stringent &
    ;;
  clean_bed_seacr)
    echo " Cleaning bedgraph for sample ${SAMPLE}"
    cat ${DATA_DIR}/seacr/${SAMPLE}.seacr.norm.stringent.stringent.bed | awk '{print $1 "\t" $2 "\t" $3}' > ${DATA_DIR}/seacr/${SAMPLE}.seacr.clean.bed
    ;;
  mergebeds)
    shift
    shift
    BEDSTOMERGE=("$@")
    OUTPUTFILE=${DATA_DIR}/counts/CR_CEBPb.norm.stringent_allMerged.bed
    echo "Merging the following beds: ${BEDSTOMERGE[@]}"
    echo "Saving to ${OUTPUTFILE}"
    cat ${BEDSTOMERGE[@]} | sort -k1,1 -k2,2n | mergeBed -i stdin >${OUTPUTFILE}
    ;;
  sort_mergedbed)
    bedtools sort -faidx /lustre/scratch117/sciops/team233/Erica/reference/chromosome.sizes -i ${DATA_DIR}/counts/CR_CEBPb.norm.stringent_allMerged.bed > ${DATA_DIR}/counts/CR_CEBPb.norm.stringent_allMerged.sorted.bed 
    ;;
  counts)
    echo "Getting counts under peaks for sample ${SAMPLE}"
    coverageBed -a ${DATA_DIR}/counts/CR_CEBPb.norm.stringent_allMerged.sorted.bed -b ${DATA_DIR}/bams/${SAMPLE}.sorted.bam -g /lustre/scratch117/sciops/team233/Erica/reference/chromosome.sizes -sorted -counts >${DATA_DIR}/counts/${SAMPLE}_CR_CEBPb_micro_PTK2B.counts &
    ;;
  *)
    echo "Unknown argument for sample ${SAMPLE}"
    ;;
  modbed)
    echo "Modifying sixth column of bed file for sample ${SAMPLE} for bigBed conversion" 
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t+"}' ${DATA_DIR}/seacr/${SAMPLE}.relaxed.bed >${DATA_DIR}/seacr/${SAMPLE}.relaxed.mod.bed
    ;;
  bigbed)
    echo "Creating bigBed for sample ${SAMPLE}"
    /home/jon_kerry/bedToBigBed ${DATA_DIR}/seacr/${SAMPLE}.relaxed.mod.bed /home/jon_kerry/CutAndRun/chromosome.sizes ${DATA_DIR}/bw/hg38/${SAMPLE}.relaxed.bb &
    ;;
  macs)
    echo "Running macs2 on sample ${SAMPLE}"
    macs2 callpeak -B -t ${DATA_DIR}/macs/${SAMPLE}_4macs.bam --name ${SAMPLE}_macs2 --outdir ${DATA_DIR}/macs -g hs -f BAMPE -q 0.01 --keep-dup all --SPMR &
    ;;
  *)
    echo "Unknown argument for sample ${SAMPLE}"
    ;;
esac
# bedtools sort -faidx ${DATA_DIR}/chromosome.sizes -i ${DATA_DIR}/${FILE_STEM}_macs2_peaks.narrowPeak > ${DATA_DIR}/${FILE_STEM}_macs2_sorted.narrowPeak
# coverageBed -a ${DATA_DIR}/${FILE_STEM}_macs2_sorted.narrowPeak -b ${DATA_DIR}/${FILE_STEM}_4macs.bam -g ${DATA_DIR}/chromosome.sizes -sorted -counts > ${DATA_DIR}/${FILE_STEM}_counts
# bedtools bamtobed -bedpe -i ${DATA_DIR}/${FILE_STEM}.bam > ${DATA_DIR}/${FILE_STEM}.bed
# bedtools genomecov -bg -i ${DATA_DIR}/${FILE_STEM}.fragments.bed -g ${DATA_DIR}/chromosome.sizes > ${DATA_DIR}/${FILE_STEM}.fragments.bedgraph
# plotHeatmap -m matrix_k27_ntc.mat.gz -out ExampleHeatmap11.png --startLabel "" --endLabel "" -y "Signal" -x "Peak" -z "K4me3Peaks" --legendLocation none --plotType std --yMin 0 --yMax 30 --samplesLabel "K27ac NTC1 1" "K27ac NTC1 2" "K27ac NTC1 3"
