source activate /mnt/data/tools
cd /mnt/data/data
s3bucket="s3://lucymli/bangladesh_chikv/"
aws s3 sync $s3bucket .
muscle -in inseq.fasta -out ncbi_chrf_aln.fasta
aws s3 cp ncbi_chrf_aln.fasta $s3bucket
# curl "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta&id=296124573&conwithfeat=on&withparts=on" -o "outgroup_HM045817.fasta" -L
# sed -i "1s/.*/>outgroup_HM045817_Senegal_2005/" outgroup_HM045817.fasta
# muscle -profile -in1 ncbi_chrf_aln.fasta -in2 outgroup_HM045817.fasta -out ncbi_chrf_aln_outgroup.fasta
Rscript cut_alignment.R
modeltest-ng --force -p 2 -q partition_modeltest.txt -i ncbi_chrf_aln.fasta
mv ncbi_chrf_aln.fasta.part.aicc partition_raxml.txt
mv ncbi_chrf_aln.fasta.part.bic partition_beast.txt
for x in `ls partition*.txt`; do aws s3 cp $x $s3bucket; done
aws s3 cp ncbi_chrf_aln.fasta.log $s3bucket
# use HM045817 as an outgroup as it is probably african strain and only has 85% sequence identity to the Bangladesh sequences. Collected Nov-2005
# curl "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?tool=portal&save=file&log$=seqview&db=nuccore&report=fasta_cds_na&id=296124573&conwithfeat=on&withparts=on" -o "outgroup_HM045817.fasta" -L
# head -n 108 outgroup_HM045817.fasta > outgroup_gp1_HM045817.fasta
# sed -n -e '109,164p' outgroup_HM045817.fasta > outgroup_gp2_HM045817.fasta
# sed -i "1s/.*/>outgroup_gp1_HM045817_Senegal_2005/" outgroup_gp1_HM045817.fasta
# sed -i "1s/.*/>outgroup_gp2_HM045817_Senegal_2005/" outgroup_gp2_HM045817.fasta
../raxml-ng/raxml-ng --all --msa ncbi_chrf_aln.fasta --model partition_raxml.txt --bs-trees 200
for x in `ls *raxml*`; do aws s3 cp $x $s3bucket; done