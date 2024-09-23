
#!/bin/zsh

cat SRR_Acc_List.txt | 
while read SRR_ID
do
	cmd='fastq-dump --split-files ${SRR_ID};
	gzip ${SRR_ID}*fastq';
	eval ${cmd}
done