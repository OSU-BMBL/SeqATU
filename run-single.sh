#/bin/nash/

##Download and install bwa

# Download and install bwa
cd data
wget https://nchc.dl.sourceforge.net/project/bio-bwa/bwa-0.7.17.tar.bz2           
tar -xjvf bwa-0.7.17.tar.bz2
rm bwa-0.7.17.tar.bz2
cd bwa-0.7.17
make
cd ..
echo "The installation of BWA was successful!"

# Run bwa
	#1：sequence_ecoli_NCBI.fasta； 2：Condition1_R1.fastq； 3: Condition1.sam
		
bwa-0.7.17/bwa index -a bwtsw $1
bwa-0.7.17/bwa aln $1 $2 > R1.sai
bwa-0.7.17/bwa samse $1 R1.sai $2 > $3 
cd ..
echo "The alignment result (.sam) was acquired successfully!"