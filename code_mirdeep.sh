#Initial Predictions with cleaned data file 
#create config file as follows
#path to file					three letter code 							
ES_1003_S_cutadapt.fastq        CT1 												   
ES_1007_S_cutadapt.fastq        CT2
ES_1594_S_cutadapt.fastq        ES4
ES_1595_S_cutadapt.fastq        ES5

# download mirBase mature and precursor sequences
wget --directory-prefix=/data/references/horse/ncRNA ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
wget --directory-prefix=/data/references/horse/ncRNA ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

#Downloaded files have white spaces that need to be removed (for miRDeep2)
zcat hairpin.fa.gz | cut -f1 -d" ">$RESULTS_DIR/hairpin_noSpace.fa
zcat mature.fa.gz | cut -f1 -d" " >$RESULTS_DIR/mature_noSpace.fa

# Remove non-canonical nucleotides and convert to SL | fastaparse.pl is part of the miRDeep 2 package | fasta_formatter is part of FASTX-Toolkit
fastaparse.pl $RESULTS_DIR/hairpin_noSpace.fa -b | fasta_formatter -w 0  >hairpin_cleaned.fa
fastaparse.pl $RESULTS_DIR/mature_noSpace.fa -b | fasta_formatter -w 0 >mature_cleaned.fa

#grep out the horse miRNAscd ..
grep -A1 "eca" mature_cleaned.fa  | grep -v -- "^--$" >eca_matureMirnas.fa
grep -A1 "eca" hairpin_cleaned.fa | grep -v -- "^--$" >eca_hairpinMirnas.fa

#Begin miRDeep process with mapper.pl (with config file)
mapper.pl config.txt -d -v -u -e -h -o 16 -m -n -q -p /data/references/horse/igenomes/Equus_caballus/Ensembl/EquCab2/Sequence/BowtieIndex/genome -s 	-t ES_S_reads_collapsed_vs_genome.arf

# Identify conserved miRNAs and predict novel miRNAs with miRDeep2
miRDeep2.pl ES_S_reads_collapsed.fa /data/references/horse/igenomes/Equus_caballus/Ensembl/EquCab2/Sequence/BowtieIndex/genome.fa ES_S_reads_collapsed_vs_genome.arf /data/references/horse/ncRNA/mature_eca_dna.fa /data/references/horse/ncRNA/mature_hsa_dna.fa /data/references/horse/ncRNA/hairpin_eca_dna.fa -t eca -z .initialPred -P 2>ES_S_report.log

# filter by mirdeep2 score cutoff ;in our case a score of 4 was applied
name_for_split=$(cat result_22_11_2016_t_14_17_20.bed|awk '{if ( $5<4)print}'|head -n 1|cut -f4)
csplit result_22_11_2016_t_14_17_20.bed “/${name_for_split}/“
mv xx00 novel_filtered.bed

## Remove predicted novel miRNAs that have a hit to Rfam (all of them showed non significant randfold p-values, so no additional step required)
cat result_22_11_2016_t_14_17_20.csv|grep RNA|grep '^[0-9]'|cut -f1>rfam_names  #rfam names has now just the names of the mirna
for line in $(cat rfam_names);do grep -v $line novel_filtered.bed>tmp;cat tmp>novel_filtered.bed;done

#merge filtered novel and mirbase (known) mirnas
cat novel_filtered.bed known.bed >known_novel_mirna.bed

# Remove overlapping miRNA predictions
# use bedtools to sort miRNA predictions by score
# the above step actually involves two steps 
# first overlap predicted miRNAs then overlap the eca miRNAs with the merged predictions.

# intially test that there are not any overlaps between known and novel
# miRNAs
# check the bedtools version, since the latest version does not have the -nms option
known_novel_mirna.bed | bedtools sort -i - |bedtools merge -d 1 -nms | grep known | grep novel  

# Now merge only the novel miRNAs 
known_novel_mirna.bed |cat ov >all_miRNAoverlap.bed

#Run quantifier with the generated combined set of known and novel miRNAs:
quantifier.pl -p precursors_Horse.fa -m Eca_mature.fa -r ALL_collapsed.fa -t Horse -k -j

#the resulting count file was used for differential expression analysis with DESeq2 and EdgeR (code provided in S5 Fil)

asympt vs c
eca-miR-10_78
eca-miR-14_7387
eca-miR-15_9066
eca-miR-20_16663
eca-miR-6_33076
eca-miR-7_34136
eca-miR-9_38236
eca-miR-142-3p
eca-miR-193a-3p
eca-miR-363
eca-miR-505

sympt vs n
eca-miR-142-3p
eca-miR-212
eca-miR-26a
eca-miR-31
eca-miR-502-3p

asy vs sy
eca-miR-130b
eca-miR-146b-5p
eca-miR-155
eca-miR-212
eca-miR-31
