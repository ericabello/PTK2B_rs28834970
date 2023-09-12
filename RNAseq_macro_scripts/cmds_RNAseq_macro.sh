DATA_DIR=/lustre/scratch117/sciops/team233/Erica/RNA_seq_PTK2B_macro_2019
SAMPLES=( KOLF2_1 KOLF2_2 KOLF2_3 C11_1 C11_2 C11_3 B11_1 B11_2 B11_3 D3_1 D3_2 D3_3 C8_1 C8_2 C8_3 )
for SAMPLE in ${SAMPLES[@]}; do
 #cramtobam)
	echo "converting ${SAMPLE} cram into bam file"
	samtools view -b ${DATA_DIR}/crams/${SAMPLE}.cram > ${DATA_DIR}/bams/${SAMPLE}.bam &
	;;
  #featureCounts)
	echo "performing FeatureCounts for all samples"
	/warehouse/team233_wh01/Erica/subread-1.6.3-Linux-i386/bin/featureCounts -T 5 -a /lustre/scratch115/realdata/mdt3/projects/otcoregen/jeremys/reference/GRCh38/gencode.v27.annotation.gtf -t exon -g gene_id -o counts_ptk2b_macro_271119.txt -p -s 0 B11_1_macro_PTK2B.bam B11_2_macro_PTK2B.bam B11_3_macro_PTK2B.bam KOLF2_1_macro_PTK2B.bam KOLF2_2_macro_PTK2B.bam KOLF2_3_macro_PTK2B.bam C11_1_macro_PTK2B.bam C11_2_macro_PTK2B.bam C11_3_macro_PTK2B.bam D3_1_macro_PTK2B.bam D3_2_macro_PTK2B.bam D3_3_macro_PTK2B.bam C8_1_macro_PTK2B.bam C8_2_macro_PTK2B.bam C8_3_macro_PTK2B.bam &
  #sort)
	echo "Sorting bam"
	samtools view -b -hu -f 2 -F 4 -F 256 -F 2048 -q 30 ${SAMPLE}_macro_PTK2B.bam | samtools sort -o ${SAMPLE}_macro_PTK2B_RNA_sorted.bam &
  #index)
	echo "Indexing bam"
	samtools index ${SAMPLE}_macro_PTK2B_RNA_sorted.bam &
	;;
  #bigwig)
    	echo "Making bigwig for sample ${SAMPLE} CPM normalised"
	bamCoverage --bam ${SAMPLE}_macro_PTK2B_RNA_sorted.bam -o ${SAMPLE}RNA_macro.SeqDepthNorm_CPM.bw --binSize 1  --normalizeUsing CPM --effectiveGenomeSize 2913022398 --ignoreForNormalization chrX &  	
    	;;

  