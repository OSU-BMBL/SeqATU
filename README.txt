SeqATU
-----------------
Brief Description:

The SeqATU program is the program to predict ATUs (alternative transcription units) of bacterial organisms.
----------



-----------------
Environment:

SeqATU is an integrated C++ package requires a basic UNIX/Linux environment. The gcc compiler with version 4.4.7 or 4.8.5 is required to be prior installed. More details can be found https://gcc.gnu.org/wiki/InstallingGCC. Currently, SeqATU does not support Mac or Windows system.
-----------------



-----------------
Usage:

1. Download SeqATU
	The source code of SeqATU is freely available at https://github.com/OSU-BMBL/SeqATU.
	To install the SeqATU, first, download the zip file manually from github, or use the code below in Unix:
		cd /your/working/path/
		wget https://github.com/OSU-BMBL/SeqATU
	Unzip the file:
		unzip SeqATU.zip && rm -rf SeqATU.zip
	
2. Prepare data
	First create a folder named data, then open it:
		cd SeqATU
		mkdir data
		cd data
	(1) Download and unzip the input data file of SeqATU
		wget https://bmbl.bmi.osumc.edu/downloadFiles/data.zip
		unzip data.zip && rm -rf data.zip
	The input data of SeqATU contain five datasets:
		The pair-end RNA-Seq datasets (Condition1_R1.fastq and Condition1_R2.fastq), genomic sequence data (sequence_ecoli_NCBI.fasta), gene annotation file (GCF_000005845.2_ASM584v2_genomic.gff) and maximal ATU cluster data (Condition1_FinalTUTable_forPlot).
	Note: you can also run the rSeqTU manually from https://github.com/OSU-BMBL/rSeqTU to get the maximal ATU cluster data (Condition1_FinalTUTable_forPlot).
	         The input data of rseqTU are GCF_000005845.2_ASM584v2_genomic.gff, sequence_ecoli_NCBI.fasta, Condition1_R1.fastq, and Condition1_R2.fastq. Choose paired="fr" of the function of qAlign for rSeqTU. 
	(2) Download, install and run bwa, after this step you will acquire the read alignment result (Condition1.sam) of the two RNA-Seq datasets (Condition1_R1.fastq and Condition1_R2.fastq):
	Input:
		The genomic sequence of Escherichia coli str. K-12 substr. MG1655 (sequence_ecoli_NCBI.fasta), and the pair-end RNA-Seq data under Condition1 (Condition1_R1.fastq and Condition1_R2.fastq).
	Output:
		The read alignment result of RNA-Seq data under Condition1 (Condition1.sam).
	Argument:
		cd ..
		sh run.sh sequence_ecoli_NCBI.fasta Condition1_R1.fastq Condition1_R2.fastq Condition1.sam
		
3. Run SeqATU
	First create a folder named out, and open the code folder:
		mkdir out
		cd code
	Then, run the SeqATU:
	Input:
		Three input datasets: the Read alignment results of RNA-Seq data under Condition1 (Condition1.sam), the gene annotation file of Escherichia coli str. K-12 substr. MG1655 (GCF_000005845.2_ASM584v2_genomic.gff), and the maximal ATU cluster data from rseqATU (Condition1_FinalTUTable_forPlot). 
		Two cutoffs: the cutoff of the distance between two maximal ATU clusters (40) and the cutoff of the expression value of maixmal ATU clusters (10).
	Output:
		Predicted ATUs under Condition1 (data_atu_Condition1.txt)
	Argument:	
		g++ SeqATU.cc Array.cc QuadProg++.cc -o SolveQuadProg
		./SolveQuadProg ../data/Condition1.sam ../data/GCF_000005845.2_ASM584v2_genomic.gff ../data/Condition1_FinalTUTable_forPlot ../out/data_atu_Condition1.txt 40 10
	Finally, you can obtain the predicted ATUs by SeqATU.
-----------------



-----------------
Data Interpretation:

1. Maximal ATU cluster data from rseqATU (Condition1_FinalTUTable_forPlot):

e.g.:
     TU name      start         end         strand       confidence          gene names
41:  TU0041      117109   117549       -               1                 ASAP:ABE-0000372,ECOCYC:EG12107,EcoGene:EG12107,GeneID:945817
42:  TU0042      117752   121551       -               1                 ASAP:ABE-0000374,ECOCYC:EG11546,EcoGene:EG11546,GeneID:948869     ASAP:ABE-0000384,ECOCYC:EG10084,EcoGene:EG10084,GeneID:946018
	
2. Predicted ATUs by SeqATU (data_atu_Condition1.txt):

1: A maximal ATU cluster and its strand 
2: ATUs under the corresponding maximal ATU cluster
3: A maximal ATU cluster and its strand
4: ATUs under the corresponding maximal ATU cluster
e.g. :
1: fhuA	fhuC	fhuD	fhuB	forward
2: fhuA	fhuA|fhuC 	fhuA|fhuC|fhuD|fhuB	
3: clcA	yadW	erpA	forward
4: clcA|yadW|erpA	yadW|erpA	erpA	
-----------------	