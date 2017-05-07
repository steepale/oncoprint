# Determine orthologues of mutated genes

cd /home/proj/MDW_genomics/steepale/oncoprint

# Files with annotated variants
./data/somatic_snvs_final.txt
./data/somatic_indels_final.txt
# Annotated CNVs and LOH
./data/scnas_annotation.bed
# Chicken-human orthologs from ensembl
./data/biomart_chicken-human_orthologs_2017_01_04.txt 



# Prepare data in proper format for oncoprint, essentially in a matrix

# Prepare SNVs

# Remove rows with Gene names of "NA"
head -n 1 ./data/somatic_snvs_final.txt > ./data/somatic_snvs_ohne_NA.tsv 
for gene in `grep -v "^#" ./data/somatic_snvs_final.txt | cut -f7 | sort | uniq | grep -v "NA"`
do 
grep -w "$gene" ./data/somatic_snvs_final.txt >> ./data/somatic_snvs_ohne_NA.tsv
done

# Generate the header for the matrix (sample names)
# Format sought ,s1,s2,s3
grep -v "^#" /home/proj/MDW_genomics/steepale/pathway_analysis/layout_uniq.txt | cut -f14 | sort | uniq | tr '\n' '\t' | sed 's/^/GENE_ID\t/' | sed 's/\t$/$/' > test.txt
# Output: ,738-1_S1,741-1_S2,756-3_S3,766-1_S4,777-3_S14,787-2_S15,788-1_S16,794-1_S17,798-1_S5,833-1_S6,834-2_2_S12,834-2_S7,835-1_S18,841-3_S19,842-2_2_S25,842-2_S20,855-1_S8,863-1_S9,884-2_S21,901-2_2_S26,901-2_S22,906-1_S23,911-1_2_S13,911-1_S24,918-3_S10,927-2_S11

# Create list of sample names in the same order
grep -v "^#" /home/proj/MDW_genomics/steepale/pathway_analysis/layout_uniq.txt | cut -f14 | sort | uniq > ./data/samples_column_order.txt

C1: 738-1_S1
C2: 741-1_S2
C3: 756-3_S3
C4: 766-1_S4
C5: 777-3_S14
C6: 787-2_S15
C7: 788-1_S16
C8: 794-1_S17
C9: 798-1_S5
C10: 833-1_S6
C11: 834-2_2_S12
C12: 834-2_S7
C13: 835-1_S18
C14: 841-3_S19
C15: 842-2_2_S25
C16: 842-2_S20
C17: 855-1_S8
C18: 863-1_S9
C19: 884-2_S21
C20: 901-2_2_S26
C21: 901-2_S22
C22: 906-1_S23
C23: 911-1_2_S13
C24: 911-1_S24
C25: 918-3_S10
C26: 927-2_S11

declare -a gene_list

gene_list[1]="ADAMTS13"
gene_list[2]="ADRA1B"
gene_list[3]="AKAP9"
gene_list[4]="AKAP9"
gene_list[5]="ANGPTL2"
gene_list[6]="ATP8B3"
gene_list[7]="ATVR1"
gene_list[8]="BCORL1"
gene_list[9]="CHST1"
gene_list[10]="CLASP2"
gene_list[11]="DHX35"
gene_list[12]="DMXL1"
gene_list[13]="DPH6"
gene_list[14]="GOLGB1"
gene_list[15]="IKZF1"
gene_list[16]="JAK2"
gene_list[17]="KCNMA1"
gene_list[18]="LAS1L"
gene_list[19]="LRP4"
gene_list[20]="LZTR1"
gene_list[21]="MED12"
gene_list[22]="MYOD1"
gene_list[23]="N4BP1"
gene_list[24]="NES"
gene_list[25]="NETO2"
gene_list[26]="NFKB1"
gene_list[27]="NOTCH2"
gene_list[28]="NSD1"
gene_list[29]="NTRK1"
gene_list[30]="PHC1"
gene_list[31]="PHKB"
gene_list[32]="PKP4"
gene_list[33]="POLK"
gene_list[34]="PRDM12"
gene_list[35]="PRDM16"
gene_list[36]="PRKAA1"
gene_list[37]="PRKAR1A"
gene_list[38]="PROM2"
gene_list[39]="RABEP1"
gene_list[40]="RASSF9"
gene_list[41]="RBM15"
gene_list[42]="RHOBTB2"
gene_list[43]="RUNX1T1"
gene_list[44]="SLC6A17"
gene_list[45]="SLX4"
gene_list[46]="SPTAN1"
gene_list[47]="SUFU"
gene_list[48]="TCEA1"
gene_list[49]="THAP1"
gene_list[50]="THSD7B"
gene_list[51]="TOR3A"
gene_list[52]="UCK1"
gene_list[53]="VCAN"
gene_list[54]="VPS18"
gene_list[55]="WDR17"
gene_list[56]="ZNF143"
gene_list[57]="ZNF367"
gene_list[58]="ZNF82"

# Check which of these genes is in CGC
for qgene in ${gene_list[@]}
do
#grep "^$qgene" ./data/cosmic_CGC_gene_list_2017_01_03.tsv
grep "$qgene" ./data/annotated_snps_indels/all_nonsyn_snps.txt
done

# Generate a tab seperated matrix to be loaded into R
# Header
grep -v "^#" /home/proj/MDW_genomics/steepale/pathway_analysis/layout_uniq.txt | cut -f14 | sort | uniq | tr '\n' '\t' | sed 's/^/GENE_ID\t/' | sed 's/\t$//' > test.txt
printf "\n" >> test.txt
# Body
for qgene in ${gene_list[@]}
do
if (grep -w -q "$qgene" ./data/somatic_snvs_ohne_NA.tsv); then
printf "$qgene\t" >> test.txt
line=`grep -w "$qgene" ./data/somatic_snvs_ohne_NA.tsv`
presample=`echo "$line" | cut -f11`
snum=`echo "$presample" | grep -o '|' | wc -l`
((snum++))
unset 'zample[*]'
declare -a zample

echo "snum: $snum"

while [ $snum -gt 0 ]
do
zample[$snum]=`echo "$presample" | cut -d '|' -f "$snum"`
((snum--))
done

echo "samples: ${zample[@]}"

if [[ "${zample[@]}" =~ "738-1_S1" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "741-1_S2" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "756-3_S3" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "766-1_S4" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "777-3_S14" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "787-2_S15" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "788-1_S16" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "794-1_S17" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "798-1_S5" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "833-1_S6" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_2_S12" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_S7" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "835-1_S18" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "841-3_S19" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_2_S25" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_S20" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "855-1_S8" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "863-1_S9" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "884-2_S21" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_2_S26" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_S22" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "906-1_S23" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_2_S13" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_S24" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "918-3_S10" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "927-2_S11" ]]; then
printf "SNV\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi

printf "\n" >> test.txt

else
printf "$qgene\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\n" >> test.txt
fi
done <./data/somatic_snvs_ohne_NA.tsv

(grep "GENE_ID" test.txt; grep -v "GENE_ID" test.txt | sort | uniq) > test2.txt

grep -v "^SNV" test2.txt | grep -v "^NULL" | sed "s/\t$//g" > ./data/somatic_snvs_oncoprint.mat


# Prepare INDELs

# Remove rows with Gene names of "NA"
head -n 1 ./data/somatic_indels_final.txt > ./data/somatic_indels_ohne_NA.tsv 
for gene in `grep -v "^#" ./data/somatic_indels_final.txt | cut -f7 | sort | uniq | grep -v "NA"`
do 
grep -w "$gene" ./data/somatic_indels_final.txt >> ./data/somatic_indels_ohne_NA.tsv
done

# Generate a tab seperated matrix to be loaded into R
# Header
grep -v "^#" /home/proj/MDW_genomics/steepale/pathway_analysis/layout_uniq.txt | cut -f14 | sort | uniq | tr '\n' '\t' | sed 's/^/GENE_ID\t/' | sed 's/\t$//' > test.txt
printf "\n" >> test.txt
# Body
for qgene in ${gene_list[@]}
do
if (grep -w -q "$qgene" ./data/somatic_indels_ohne_NA.tsv); then
printf "$qgene\t" >> test.txt
line=`grep -w "$qgene" ./data/somatic_indels_ohne_NA.tsv`
presample=`echo "$line" | cut -f11`
snum=`echo "$presample" | grep -o '|' | wc -l`
((snum++))
unset 'zample[*]'
declare -a zample

echo "snum: $snum"

while [ $snum -gt 0 ]
do
zample[$snum]=`echo "$presample" | cut -d '|' -f "$snum"`
((snum--))
done

echo "samples: ${zample[@]}"

if [[ "${zample[@]}" =~ "738-1_S1" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "741-1_S2" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "756-3_S3" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "766-1_S4" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "777-3_S14" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "787-2_S15" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "788-1_S16" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "794-1_S17" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "798-1_S5" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "833-1_S6" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_2_S12" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_S7" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "835-1_S18" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "841-3_S19" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_2_S25" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_S20" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "855-1_S8" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "863-1_S9" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "884-2_S21" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_2_S26" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_S22" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "906-1_S23" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_2_S13" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_S24" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "918-3_S10" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "927-2_S11" ]]; then
printf "INDEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi

printf "\n" >> test.txt

else
printf "$qgene\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\n" >> test.txt
fi
done <./data/somatic_indels_ohne_NA.tsv

(grep "GENE_ID" test.txt; grep -v "GENE_ID" test.txt | sort | uniq) > test2.txt

grep -v "^INDEL" test2.txt | grep -v "^NULL" | sed "s/\t$//g" > ./data/somatic_indels_oncoprint.mat


# Prepare AMPs

# Adjust the sample names accordingly
# Replace all of the samples ID's using an array and sed to replace the word only strings in file
declare -A birds
birds["S1"]="738-1_S1"
birds["S2"]="741-1_S2"
birds["S3"]="756-3_S3"
birds["S4"]="766-1_S4"
birds["S5"]="798-1_S5"
birds["S6"]="833-1_S6"
birds["S7"]="834-2_S7"
birds["S8"]="855-1_S8"
birds["S9"]="863-1_S9"
birds["S10"]="918-3_S10"
birds["S11"]="927-2_S11"
birds["S12"]="834-2_2_S12"
birds["S13"]="911-1_2_S13"
birds["S14"]="777-3_S14"
birds["S15"]="787-2_S15"
birds["S16"]="788-1_S16"
birds["S17"]="794-1_S17"
birds["S18"]="835-1_S18"
birds["S19"]="841-3_S19"
birds["S20"]="842-2_S20"
birds["S21"]="884-2_S21"
birds["S22"]="901-2_S22"
birds["S23"]="906-1_S23"
birds["S24"]="911-1_S24"
birds["S25"]="842-2_2_S25"
birds["S26"]="901-2_2_S26"

for n in "${!birds[@]}"
do
echo "$n - ${birds[$n]}"
sed -i "s/\b$n\b/${birds[$n]}/g" ./data/scnas_annotation.bed
done



# Reduce the file to only uniq values

# Create a file with just the amplifications
cut -f4,5,11 ./data/scnas_annotation.bed | sort | uniq | grep "^gain" > ./data/amps.txt

# get uniq genes
# Header
#printf "#GENE_SYMBOL\tALT\tSAMPLES\n" > ./data/amps_uniq.txt
rm ./data/amps_uniq.txt
for gene in `cut -f3 ./data/amps.txt | sort | uniq | grep -w -v '\-1'`
do
sample=`grep -w "$gene" ./data/amps.txt | cut -f2 | awk '!seen[$0]++'`
#echo "gene: $gene"
#echo "sample: $sample"
samples=`echo -e "$sample" | tr '\n' '|' | sed 's/|$//'`
#echo "samples: $samples"
printf "$gene\tAMP\t$samples\n" >> ./data/amps_uniq.txt
done


declare -a gene_list

gene_list[1]="ADAMTS13"
gene_list[2]="ADRA1B"
gene_list[3]="AKAP9"
gene_list[4]="AKAP9"
gene_list[5]="ANGPTL2"
gene_list[6]="ATP8B3"
gene_list[7]="ATVR1"
gene_list[8]="BCORL1"
gene_list[9]="CHST1"
gene_list[10]="CLASP2"
gene_list[11]="DHX35"
gene_list[12]="DMXL1"
gene_list[13]="DPH6"
gene_list[14]="GOLGB1"
gene_list[15]="IKZF1"
gene_list[16]="JAK2"
gene_list[17]="KCNMA1"
gene_list[18]="LAS1L"
gene_list[19]="LRP4"
gene_list[20]="LZTR1"
gene_list[21]="MED12"
gene_list[22]="MYOD1"
gene_list[23]="N4BP1"
gene_list[24]="NES"
gene_list[25]="NETO2"
gene_list[26]="NFKB1"
gene_list[27]="NOTCH2"
gene_list[28]="NSD1"
gene_list[29]="NTRK1"
gene_list[30]="PHC1"
gene_list[31]="PHKB"
gene_list[32]="PKP4"
gene_list[33]="POLK"
gene_list[34]="PRDM12"
gene_list[35]="PRDM16"
gene_list[36]="PRKAA1"
gene_list[37]="PRKAR1A"
gene_list[38]="PROM2"
gene_list[39]="RABEP1"
gene_list[40]="RASSF9"
gene_list[41]="RBM15"
gene_list[42]="RHOBTB2"
gene_list[43]="RUNX1T1"
gene_list[44]="SLC6A17"
gene_list[45]="SLX4"
gene_list[46]="SPTAN1"
gene_list[47]="SUFU"
gene_list[48]="TCEA1"
gene_list[49]="THAP1"
gene_list[50]="THSD7B"
gene_list[51]="TOR3A"
gene_list[52]="UCK1"
gene_list[53]="VCAN"
gene_list[54]="VPS18"
gene_list[55]="WDR17"
gene_list[56]="ZNF143"
gene_list[57]="ZNF367"
gene_list[58]="ZNF82"

# Generate a tab seperated matrix to be loaded into R
# Header
grep -v "^#" /home/proj/MDW_genomics/steepale/pathway_analysis/layout_uniq.txt | cut -f14 | sort | uniq | tr '\n' '\t' | sed 's/^/GENE_ID\t/' | sed 's/\t$//' > test.txt
printf "\n" >> test.txt
# Body

for qgene in ${gene_list[@]}
do
if (grep -w -q "$qgene" ./data/amps_uniq.txt); then
printf "$qgene\t" >> test.txt
line=`grep -w "$qgene" ./data/amps_uniq.txt`
presample=`echo "$line" | cut -f3`
snum=`echo "$presample" | grep -o '|' | wc -l`
((snum++))
unset 'zample[*]'
declare -a zample

echo "snum: $snum"

while [ $snum -gt 0 ]
do
zample[$snum]=`echo "$presample" | cut -d '|' -f "$snum"`
((snum--))
done

echo "samples: ${zample[@]}"

if [[ "${zample[@]}" =~ "738-1_S1" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "741-1_S2" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "756-3_S3" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "766-1_S4" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "777-3_S14" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "787-2_S15" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "788-1_S16" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "794-1_S17" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "798-1_S5" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "833-1_S6" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_2_S12" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_S7" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "835-1_S18" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "841-3_S19" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_2_S25" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_S20" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "855-1_S8" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "863-1_S9" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "884-2_S21" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_2_S26" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_S22" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "906-1_S23" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_2_S13" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_S24" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "918-3_S10" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "927-2_S11" ]]; then
printf "AMP\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi

printf "\n" >> test.txt

else
printf "$qgene\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\n" >> test.txt
fi
done <./data/amps_uniq.txt

(grep "GENE_ID" test.txt; grep -v "GENE_ID" test.txt | sort | uniq) > test2.txt

grep -v "^AMP" test2.txt | grep -v "^NULL" | sed "s/\t$//g" > ./data/somatic_amps_oncoprint.mat


# Prepare DELs

# Reduce the file to only uniq values

# Create a file with just the deletions
cut -f4,5,11 ./data/scnas_annotation.bed | sort | uniq | grep "^loss" > ./data/dels.txt

# get uniq genes
# Header
#printf "#GENE_SYMBOL\tALT\tSAMPLES\n" > ./data/amps_uniq.txt
rm ./data/dels_uniq.txt
for gene in `cut -f3 ./data/dels.txt | sort | uniq | grep -w -v '\-1'`
do
sample=`grep -w "$gene" ./data/dels.txt | cut -f2 | awk '!seen[$0]++'`
#echo "gene: $gene"
#echo "sample: $sample"
samples=`echo -e "$sample" | tr '\n' '|' | sed 's/|$//'`
#echo "samples: $samples"
printf "$gene\tDEL\t$samples\n" >> ./data/dels_uniq.txt
done

declare -a gene_list

gene_list[1]="ADAMTS13"
gene_list[2]="ADRA1B"
gene_list[3]="AKAP9"
gene_list[4]="AKAP9"
gene_list[5]="ANGPTL2"
gene_list[6]="ATP8B3"
gene_list[7]="ATVR1"
gene_list[8]="BCORL1"
gene_list[9]="CHST1"
gene_list[10]="CLASP2"
gene_list[11]="DHX35"
gene_list[12]="DMXL1"
gene_list[13]="DPH6"
gene_list[14]="GOLGB1"
gene_list[15]="IKZF1"
gene_list[16]="JAK2"
gene_list[17]="KCNMA1"
gene_list[18]="LAS1L"
gene_list[19]="LRP4"
gene_list[20]="LZTR1"
gene_list[21]="MED12"
gene_list[22]="MYOD1"
gene_list[23]="N4BP1"
gene_list[24]="NES"
gene_list[25]="NETO2"
gene_list[26]="NFKB1"
gene_list[27]="NOTCH2"
gene_list[28]="NSD1"
gene_list[29]="NTRK1"
gene_list[30]="PHC1"
gene_list[31]="PHKB"
gene_list[32]="PKP4"
gene_list[33]="POLK"
gene_list[34]="PRDM12"
gene_list[35]="PRDM16"
gene_list[36]="PRKAA1"
gene_list[37]="PRKAR1A"
gene_list[38]="PROM2"
gene_list[39]="RABEP1"
gene_list[40]="RASSF9"
gene_list[41]="RBM15"
gene_list[42]="RHOBTB2"
gene_list[43]="RUNX1T1"
gene_list[44]="SLC6A17"
gene_list[45]="SLX4"
gene_list[46]="SPTAN1"
gene_list[47]="SUFU"
gene_list[48]="TCEA1"
gene_list[49]="THAP1"
gene_list[50]="THSD7B"
gene_list[51]="TOR3A"
gene_list[52]="UCK1"
gene_list[53]="VCAN"
gene_list[54]="VPS18"
gene_list[55]="WDR17"
gene_list[56]="ZNF143"
gene_list[57]="ZNF367"
gene_list[58]="ZNF82"

# Generate a tab seperated matrix to be loaded into R
# Header
grep -v "^#" /home/proj/MDW_genomics/steepale/pathway_analysis/layout_uniq.txt | cut -f14 | sort | uniq | tr '\n' '\t' | sed 's/^/GENE_ID\t/' | sed 's/\t$//' > test.txt
printf "\n" >> test.txt
# Body

for qgene in ${gene_list[@]}
do
if (grep -w -q "$qgene" ./data/dels_uniq.txt); then
printf "$qgene\t" >> test.txt
line=`grep -w "$qgene" ./data/dels_uniq.txt`
presample=`echo "$line" | cut -f3`
snum=`echo "$presample" | grep -o '|' | wc -l`
((snum++))
unset 'zample[*]'
declare -a zample

echo "snum: $snum"

while [ $snum -gt 0 ]
do
zample[$snum]=`echo "$presample" | cut -d '|' -f "$snum"`
((snum--))
done

echo "samples: ${zample[@]}"

if [[ "${zample[@]}" =~ "738-1_S1" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "741-1_S2" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "756-3_S3" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "766-1_S4" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "777-3_S14" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "787-2_S15" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "788-1_S16" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "794-1_S17" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "798-1_S5" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "833-1_S6" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_2_S12" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_S7" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "835-1_S18" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "841-3_S19" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_2_S25" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_S20" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "855-1_S8" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "863-1_S9" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "884-2_S21" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_2_S26" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_S22" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "906-1_S23" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_2_S13" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_S24" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "918-3_S10" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "927-2_S11" ]]; then
printf "DEL\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi

printf "\n" >> test.txt

else
printf "$qgene\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\n" >> test.txt
fi
done <./data/dels_uniq.txt

(grep "GENE_ID" test.txt; grep -v "GENE_ID" test.txt | sort | uniq) > test2.txt

grep -v "^DEL" test2.txt | grep -v "^NULL" | sed "s/\t$//g" > ./data/somatic_dels_oncoprint.mat


# Prepare LOHs

# Reduce the file to only uniq values

# Create a file with just the LOHs
cut -f4,5,11 ./data/scnas_annotation.bed | sort | uniq | grep "^neutral" > ./data/loh.txt

# get uniq genes
# Header
#printf "#GENE_SYMBOL\tALT\tSAMPLES\n" > ./data/amps_uniq.txt
rm ./data/loh_uniq.txt
for gene in `cut -f3 ./data/loh.txt | sort | uniq | grep -w -v '\-1'`
do
sample=`grep -w "$gene" ./data/loh.txt | cut -f2 | awk '!seen[$0]++'`
#echo "gene: $gene"
#echo "sample: $sample"
samples=`echo -e "$sample" | tr '\n' '|' | sed 's/|$//'`
#echo "samples: $samples"
printf "$gene\tLOH\t$samples\n" >> ./data/loh_uniq.txt
done

declare -a gene_list

gene_list[1]="ADAMTS13"
gene_list[2]="ADRA1B"
gene_list[3]="AKAP9"
gene_list[4]="AKAP9"
gene_list[5]="ANGPTL2"
gene_list[6]="ATP8B3"
gene_list[7]="ATVR1"
gene_list[8]="BCORL1"
gene_list[9]="CHST1"
gene_list[10]="CLASP2"
gene_list[11]="DHX35"
gene_list[12]="DMXL1"
gene_list[13]="DPH6"
gene_list[14]="GOLGB1"
gene_list[15]="IKZF1"
gene_list[16]="JAK2"
gene_list[17]="KCNMA1"
gene_list[18]="LAS1L"
gene_list[19]="LRP4"
gene_list[20]="LZTR1"
gene_list[21]="MED12"
gene_list[22]="MYOD1"
gene_list[23]="N4BP1"
gene_list[24]="NES"
gene_list[25]="NETO2"
gene_list[26]="NFKB1"
gene_list[27]="NOTCH2"
gene_list[28]="NSD1"
gene_list[29]="NTRK1"
gene_list[30]="PHC1"
gene_list[31]="PHKB"
gene_list[32]="PKP4"
gene_list[33]="POLK"
gene_list[34]="PRDM12"
gene_list[35]="PRDM16"
gene_list[36]="PRKAA1"
gene_list[37]="PRKAR1A"
gene_list[38]="PROM2"
gene_list[39]="RABEP1"
gene_list[40]="RASSF9"
gene_list[41]="RBM15"
gene_list[42]="RHOBTB2"
gene_list[43]="RUNX1T1"
gene_list[44]="SLC6A17"
gene_list[45]="SLX4"
gene_list[46]="SPTAN1"
gene_list[47]="SUFU"
gene_list[48]="TCEA1"
gene_list[49]="THAP1"
gene_list[50]="THSD7B"
gene_list[51]="TOR3A"
gene_list[52]="UCK1"
gene_list[53]="VCAN"
gene_list[54]="VPS18"
gene_list[55]="WDR17"
gene_list[56]="ZNF143"
gene_list[57]="ZNF367"
gene_list[58]="ZNF82"

# Generate a tab seperated matrix to be loaded into R
# Header
grep -v "^#" /home/proj/MDW_genomics/steepale/pathway_analysis/layout_uniq.txt | cut -f14 | sort | uniq | tr '\n' '\t' | sed 's/^/GENE_ID\t/' | sed 's/\t$//' > test.txt
printf "\n" >> test.txt
# Body
for qgene in ${gene_list[@]}
do
if (grep -w -q "$qgene" ./data/loh_uniq.txt); then
printf "$qgene\t" >> test.txt
line=`grep -w "$qgene" ./data/loh_uniq.txt`
presample=`echo "$line" | cut -f3`
snum=`echo "$presample" | grep -o '|' | wc -l`
((snum++))
unset 'zample[*]'
declare -a zample

echo "snum: $snum"

while [ $snum -gt 0 ]
do
zample[$snum]=`echo "$presample" | cut -d '|' -f "$snum"`
((snum--))
done

echo "samples: ${zample[@]}"

if [[ "${zample[@]}" =~ "738-1_S1" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "741-1_S2" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "756-3_S3" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "766-1_S4" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "777-3_S14" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "787-2_S15" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "788-1_S16" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "794-1_S17" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "798-1_S5" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "833-1_S6" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_2_S12" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "834-2_S7" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "835-1_S18" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "841-3_S19" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_2_S25" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "842-2_S20" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "855-1_S8" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "863-1_S9" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "884-2_S21" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_2_S26" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "901-2_S22" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "906-1_S23" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_2_S13" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "911-1_S24" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "918-3_S10" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi
if [[ "${zample[@]}" =~ "927-2_S11" ]]; then
printf "LOH\t" >> test.txt
else
printf "NULL\t" >> test.txt
fi

printf "\n" >> test.txt

else
printf "$qgene\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\n" >> test.txt
fi
done <./data/loh_uniq.txt

(grep "^GENE_ID" test.txt; grep -v "^GENE_ID" test.txt | sort | uniq) > test2.txt

grep -v "^LOH" test2.txt | grep -v "^NULL" | sed "s/\t$//g" > ./data/somatic_loh_oncoprint.mat


























head -n 1 ./data/biomart_chicken-human_orthologs_2017_01_04.txt > ./results/human_orthologs_mutated_genes.tsv
# Combine all of the annotated ensemble gene ID's from somatic SNVs and Indels and search for human ortholgues
for chick_ensebml in `(grep -v "^#" ./data/somatic_indels_final.txt; grep -v "^#" ./data/somatic_snvs_final.txt) | cut -f8 | sort | uniq`
do
grep -w "$chick_ensebml" ./data/biomart_chicken-human_orthologs_2017_01_04.txt >> ./results/human_orthologs_mutated_genes.tsv
done

# Check for overlap between SNVs, indels, CNVs, and LOH
rm ./results/shared_snvs_indels_cnvs_loh.bed
for chick_ensebml in `(grep -v "^#" ./data/somatic_indels_final.txt; grep -v "^#" ./data/somatic_snvs_final.txt) | cut -f7 | sort | uniq | grep -v "^NA"`
do
grep -w "$chick_ensebml" ./data/scnas_annotation.bed >> ./results/shared_snvs_indels_cnvs_loh.bed
done

head -n 1 ./data/cosmic_CGC_gene_list_2017_01_03.tsv | cut -f1,8,18 > ./results/human_orthologs_cosmic.tsv
for humo_name in `grep -v "^Human associated gene name" ./results/human_orthologs_mutated_genes.tsv | cut -f3 | sort | uniq`
do
grep "^$humo_name" ./data/cosmic_CGC_gene_list_2017_01_03.tsv | cut -f1,8,18 >> ./results/human_orthologs_cosmic.tsv
done

for gene in `grep -v "^Gene Symbol" ./results/human_orthologs_cosmic.tsv | cut -f1`
do
grep "$gene" ./data/biomart_chicken-human_orthologs_2017_01_04.txt
done