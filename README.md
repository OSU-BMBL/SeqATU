# SeqATU #
## Brief Description: ##
The SeqATU program is the program to predict ATUs (alternative transcription units) of bacterial organisms.
## Environment: ##
SeqATU is an integrated C++ package requires a basic **UNIX/Linux environment**. The gcc compiler with **version 4.4.7 or 4.8.5** is required to be prior installed. More details can be found [here](https://gcc.gnu.org/wiki/InstallingGCC). Currently, SeqATU does not support Mac or Windows system.
# Usage: #
## 1. Download SeqATU ##
The source code of SeqATU is freely available at <https://github.com/OSU-BMBL/SeqATU>.  
To install the SeqATU, first, download the zip file manually from github, or use the code below in Unix:  

    cd /your/working/path/
	wget https://github.com/OSU-BMBL/SeqATU/archive/master.zip

Unzip the file:

	unzip master.zip && rm -rf master.zip

## 2. Prepare data ##
First open SeqATU-master:

	cd SeqATU-master

1.Download and unzip the input data file of SeqATU

	wget https://bmbl.bmi.osumc.edu/downloadFiles/data.zip
	unzip data.zip && rm -rf data.zip
	

The input data of SeqATU contains five datasets:
|Data|Description|
| :---: | :--- |
| ./SeqATU-master/data/**M9Enrich_R1.fastq** | *pair-end RNA-Seq data* |
| ./SeqATU-master/data/**M9Enrich_R2.fastq** | *pair-end RNA-Seq data* |
| ./SeqATU-master/data/**sequence_ecoli_NCBI.fasta** | *genomic sequence data* |
| ./SeqATU-master/data/**GCF_000005845.2_ASM584v2_genomic.gff** | *gene annotation file* |
| ./SeqATU-master/data/**M9Enrich_FinalTUTable_forPlot** | *maximal ATU cluster data* |


**Note 1:** You can also run the rSeqTU manually from <https://github.com/OSU-BMBL/rSeqTU> to get the maximal ATU cluster data (M9Enrich_FinalTUTable_forPlot).  
**Note 2:** The input data of rSeqTU are GCF_000005845.2_ASM584v2_genomic.gff, sequence_ecoli_NCBI.fasta, M9Enrich_R1.fastq, and M9Enrich_R2.fastq. Choose paired="fr" of the function of qAlign for rSeqTU.  

2.Download, install and run bwa, after this step you will acquire the read alignment result (M9Enrich.sam) for the two RNA-Seq datasets (M9Enrich_R1.fastq and M9Enrich_R2.fastq): 
 
***Input:***  
The genomic sequence of Escherichia coli str. K-12 substr. MG1655 (sequence_ecoli_NCBI.fasta), and the pair-end RNA-Seq data under M9Enrich (M9Enrich_R1.fastq and M9Enrich_R2.fastq).  
***Output:***  
The read alignment result for RNA-Seq data under M9Enrich (M9Enrich.sam).  
***Argument:***  
	
	sh run.sh sequence_ecoli_NCBI.fasta M9Enrich_R1.fastq M9Enrich_R2.fastq M9Enrich.sam
	
## 3. Run SeqATU ##
First create a folder named out, then open the code folder:

	mkdir out && cd code

Then, run the SeqATU:

***Input:***  
Three input datasets: Read alignment results for RNA-Seq data under M9Enrich (M9Enrich.sam), the gene annotation file of Escherichia coli str. K-12 substr. MG1655 (GCF_000005845.2_ASM584v2_genomic.gff), and the maximal ATU cluster data from rseqATU (M9Enrich_FinalTUTable_forPlot).  
Two cutoffs: the cutoff of the distance between two maximal ATU clusters (40) and the cutoff of the expression value of maixmal ATU clusters (10).  
***Output:***  
Predicted ATUs under M9Enrich (data_atu_M9Enrich.txt)  
***Argument:***  

	g++ SeqATU.cc Array.cc QuadProg++.cc -o SolveQuadProg
	./SolveQuadProg ../data/M9Enrich.sam ../data/GCF_000005845.2_ASM584v2_genomic.gff ../data/M9Enrich_FinalTUTable_forPlot ../out/data_atu_M9Enrich.txt 40 10
	
Finally, you can obtain the predicted ATUs by SeqATU.

**Note :** You can also obtain the predicted ATUs by SeqATU under condition RiEnrich by using the three datasets RiEnrich_R1.fastq, RiEnrich_R2.fastq, and RiEnrich_FinalTUTable_forPlot.
# Data Interpretation: #
## 1. Maximal ATU cluster data from rseqATU (M9Enrich_FinalTUTable_forPlot): ##

| TU name | Start | End | Strand | Confidence | Gene names |
| :---: | :---: | :---: | :---: | :---: | :---: |
| TU0041 | 117109 | 117549 | - | 1 | ASAP:ABE-0000372,ECOCYC:EG12107,EcoGene:EG12107,GeneID:945817 |
| TU0042 | 117752 | 121551 | - | 1 |ASAP:ABE-0000374,ECOCYC:EG11546,EcoGene:EG11546,GeneID:948869<br>ASAP:ABE-0000384,ECOCYC:EG10084,EcoGene:EG10084,GeneID:946018 |
## 2.Predicted ATUs by SeqATU (data_atu_M9Enrich.txt): ##
| Description | Data |
| :--- | :--- |
| A maximal ATU cluster and its strand | fhuA fhuC fhuD fhuB forward |
| ATUs | fhuA fhuA&#124;fhuC fhuA&#124;fhuC&#124;fhuD&#124;fhuB	|
| A maximal ATU cluster and its strand | clcA yadW erpA forward |
| ATUs | clcA&#124;yadW&#124;erpA yadW&#124;erpA erpA |
# Contact #
Any questions, problems, bugs are welcome and should be dumped to [Qin Ma](Qin.Ma@osumc.edu).
