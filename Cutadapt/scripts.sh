
# ----- Aug 2, 2024 at 17:48:46 -----

# The phyllosphere microbiome shifts toward combating melanose pathogen
# EN_16S Forward Primer
cutadapt -a GTGYCAGCMGCCGCGGTAA -o /dev/null {SRR12131106,SRR12131107,SRR12131108,SRR12131109,SRR12131110,SRR12131111,SRR12131112,SRR12131113,SRR12131114,SRR12131154,SRR12131155,SRR12131156,SRR12131157,SRR12131158,SRR12131159,SRR12131160,SRR12131161,SRR12131162,SRR12131163,SRR12131165}_1.fastq.gz
# → 複数fileを一度には不可能

# EN_16S → アダプター配列を含むReadsそのものが除去されてしまう(-aは3'adapterに対して作用)
cutadapt -a GTGYCAGCMGCCGCGGTAA -A GACTACHVGGGTATCTAATCC -o out_SRR12131106_1.fastq.gz -p out_SRR12131106_2.fastq.gz SRR12131106_1.fastq.gz SRR12131106_2.fastq.gz > score.txt

cutadapt --revcomp -a GTGYCAGCMGCCGCGGTAA -A GACTACHVGGGTATCTAATCC -o out_SRR12131106_1.fastq.gz -p out_SRR12131106_2.fastq.gz SRR12131106_1.fastq.gz SRR12131106_2.fastq.gz > score.txt

cutadapt -a GTGYCAGCMGCCGCGGTAA -A GGACTACHVGGGTWTCTAAT -o out_SRR12131106_1.fastq.gz -p out_SRR12131106_2.fastq.gz SRR12131106_1.fastq.gz SRR12131106_2.fastq.gz > score_v1.txt

# ----- Aug 7, 2024 at 10:08:23 -----
# (-gは5'adapterに対して作用)
cutadapt -g GTGYCAGCMGCCGCGGTAA -G GACTACHVGGGTATCTAATCC -o out_SRR12131106_1.fastq.gz -p out_SRR12131106_2.fastq.gz SRR12131106_1.fastq.gz SRR12131106_2.fastq.gz --json=SRR12131106.cutadapt.json 

# /Users/saitoikuto/Documents/015.RStudio/NGS_Practice/The phyllosphere microbiome shifts toward combating melanose pathogen/EN_16S
for ID in $(cat ../SRR_Acc_List_EN.txt); do \
cutadapt -g GTGYCAGCMGCCGCGGTAA -G GACTACHVGGGTATCTAATCC --discard-untrimmed --quality-cutoff 10 --quiet \
-o P_Trim_${ID}_1.fastq.gz -p P_Trim_${ID}_2.fastq.gz \
${ID}_1.fastq.gz ${ID}_2.fastq.gz \
; done 

# FastQC CUI 
# /Users/saitoikuto/Documents/015.RStudio/NGS_Practice/The phyllosphere microbiome shifts toward combating melanose pathogen/EN_16S
for ID in $(cat ../EN_P_Trim_fastq.txt); do \
fastqc --outdir ../EN_16S_FastQC_report ${ID} \
; done


# /Users/saitoikuto/Documents/015.RStudio/NGS_Practice/The phyllosphere microbiome shifts toward combating melanose pathogen/EP_16S
for ID in $(cat ../SRR_Acc_List_EP.txt); do \
cutadapt -g GTGCCAGCMGCCGCGGTAA -G GGACTACHVGGGTWTCTAAT --discard-untrimmed --quality-cutoff 10 --quiet \
-o P_Trim_${ID}_1.fastq.gz -p P_Trim_${ID}_2.fastq.gz \
${ID}_1.fastq.gz ${ID}_2.fastq.gz \
; done 

# FastQC CUI 
# /Users/saitoikuto/Documents/015.RStudio/NGS_Practice/The phyllosphere microbiome shifts toward combating melanose pathogen/EP_16S
for ID in $(cat ../EP_P_Trim_fastq.txt); do \
fastqc --outdir ../EP_16S_FastQC_report ${ID} \
; done

# P_Trim mv
# /Users/saitoikuto/Documents/015.RStudio/NGS_Practice/The phyllosphere microbiome shifts toward combating melanose pathogen/EP_16S
mv P_Trim* ../EP_16S_P_Trim
mv P_Trim* ../EN_16S_P_Trim
rename 's/_1.fastq.gz/_R1.fastq.gz/' P_Trim_*.fastq.gz \
;rename 's/_2.fastq.gz/_R2.fastq.gz/' P_Trim_*.fastq.gz







