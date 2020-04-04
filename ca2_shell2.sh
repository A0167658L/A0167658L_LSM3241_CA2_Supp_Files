set -e
cd ~/ca2
mkdir -p fastqc
echo 'running fastqc'
fastqc ~/ca2/rawdata/*.fq* -o ~/ca2/fastqc
cd ~/ca2/fastqc

echo "unzipping"
for filename in ~/ca2/fastqc/*.zip
	do
	unzip $filename
	done

echo "saving summary"
mkdir -p ~/ca2/docs
cat */summary.txt > ~/ca2/docs/fastqc_summaries.txt

echo "trimming"
cd ~/ca2
mkdir -p ~/ca2/cleaned
mkdir -p ~/ca2/cleaned/trimmed
mkdir -p ~/ca2/cleaned/orphaned
TrimmomaticPE rawdata/A0167658L_1.fq rawdata/A0167658L_2.fq \
cleaned/trimmed/A0167658L_1.trimmed.fq cleaned/orphaned/A0167658L_1.orphaned.fq \
cleaned/trimmed/A0167658L_2.trimmed.fq cleaned/orphaned/A0167658L_2.orphaned.fq \
MINLEN:20 SLIDINGWINDOW:4:20

echo "running fastqc"
fastqc ~/ca2/cleaned/trimmed/*.fq* -o ~/ca2/fastqc
cd ~/ca2/fastqc
for filename in ~/ca2/fastqc/*trimmed*.zip
	do
	unzip $filename
	done
cat */summary.txt > ~/ca2/docs/fastqc_summaries.txt

echo "indexing"
cd ~/ca2
mkdir -p ~/ca2/allreferences/all
bowtie2-build rawdata/sacCer3.fa,rawdata/ty5_6p.fa allreferences/all/reference

echo "alignment against all"
cd ~/ca2
mkdir -p ~/ca2/results/sam

export BOWTIE2_INDEXES=$(pwd)/allreferences/all
bowtie2 -x reference \
	--sensitive -p 4 \
	-1 cleaned/trimmed/A0167658L_1* \
	-2 cleaned/trimmed/A0167658L_2* \
	-S results/sam/A0167658L_results.sam

echo "converting to bam"
mkdir -p ~/ca2/results/bam
samtools view -S -b results/sam/A0167658L_results.sam > results/bam/A0167658L_results.bam

echo "sorting"
mkdir results/bam/converttobed
samtools sort results/bam/A0167658L_results.bam -o results/bam/converttobed/A0167658L_sortedresults.bam

echo "indexing for tview"
samtools index results/bam/converttobed/A0167658L_sortedresults.bam

echo "finding unmapped reads"
samtools view -b -F 2 results/bam/converttobed/A0167658L_sortedresults.bam > results/bam/A0167658L_intersect.bam

echo "searching for transposon locations"
samtools view results/bam/A0167658L_intersect.bam | grep -i 'ty' > results/bam/A0167658L_grep.bam
samtools view -H results/bam/A0167658L_intersect.bam > results/sam/header.sam
cat results/sam/header.sam results/bam/A0167658L_grep.bam > results/sam/A0167658L_locations.sam
samtools view -S -b results/sam/A0167658L_locations.sam > results/bam/converttobed/A0167658L_locations.bam

echo "list of transposons inserted in reverse orientation"
samtools view -b -f 20 -F 10 results/bam/converttobed/A0167658L_locations.bam > results/bam/converttobed/A0167658L_orientation.bam

echo "converting to bed"
mkdir ~/ca2/results/bed
for filename in results/bam/converttobed/*.bam; do
	name=$(basename $filename .bam)
	bedtools bamtobed -i results/bam/converttobed/${name}.bam > results/bed/${name}.bed
	done

echo "sorting bedfiles"
for filename in results/bed/*.bed; do
	name=$(basename $filename .bed)
	sort -k1,1 -k2,2n results/bed/${name}.bed > results/bed/${name}_sorted.bed
	done

echo "merging bedfiles"
for filename in results/bed/*_sorted.bed; do
	name=$(basename $filename _sorted.bed)
	bedtools merge -i results/bed/${name}_sorted.bed > results/bed/${name}_merged.bed
	bedtools merge -i results/bed/${name}_sorted.bed > results/${name}.txt
	done

echo "calculating genome coverage"
for filename in rawdata/*.fa; do 
	samtools faidx $filename
	done

awk -v OFS='\t' {'print $1,$2'} rawdata/sacCer3.fa.fai > rawdata/sacCer3_genomefile.txt
awk -v OFS='\t' {'print $1,$2'} rawdata/ty5_6p.fa.fai > rawdata/ty5_6p_genomefile.txt

cat rawdata/*_genomefile.txt > rawdata/genome.txt

bedtools genomecov -max 15 -i results/bed/A0167658L_sortedresults_sorted.bed \
-g rawdata/genome.txt > results/genomecoverage.xls

