###########
### Script for counting cell-level read counts from 10x run (BAMTagHistogram)
### D.Croucher
### September 8, 2020
### https://github.com/broadinstitute/Drop-seq
############

#Run create_bamtaghistogram.sh to set up sample-specific scripts:

#samplenames=`cat sample_names.txt`

#for i in $samplenames
#do
#    sed s/XXXX/"$i"/g templates_bamtaghistogram.sh > scripts/bamtaghistogram/$i.sh
#done

# Run bamtaghistogram

module load java/7

DROPSEQ_ROOT='/cluster/home/dcrouche/pughlab/bin/Drop-seq_tools-1-2.12'
TOOL='BAMTagHistogram'

$DROPSEQ_ROOT'/'$TOOL \
I=../crcount/XXXX_CRouts/outs/possorted_genome_bam.bam \
O=../../BAMTagHistogram/XXXX.bam.readcounts.txt.gz \
TAG=CB 
