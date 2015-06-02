
all: fakeFASTQ.fastq enrichFiles forward_paired.fq reverse_paired.fq forward_unpaired.fq reverse_unpaired.fq


clean: 
	-rm -f forward_paired.fq reverse_paired.fq forward_unpaired.fq reverse_unpaired.fq fakeFASTQ.fastq out.* water.txt read.fasta

forward_paired.fq reverse_paired.fq forward_unpaired.fq reverse_unpaired.fq:
	java -jar /Users/jimstapleton/scientificPrograms/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 20 test1.fastq test2.fastq forward_paired.fq forward_unpaired.fq reverse_paired.fq reverse_unpaired.fq ILLUMINACLIP:/Users/jimstapleton/scientificPrograms/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:6:1 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20

fakeFASTQ.fastq:
	python /Users/jimstapleton/Dropbox/Whitehead_lab/Maynard/mutationCounter.py forward_paired.fq forward_unpaired.fq reverse_paired.fq reverse_unpaired.fq
	rm -f out.hist out.histogram water.txt read.fasta

enrichFiles:
	cp fakeFASTQ.fastq ./Enrich/M1H5-P5/data/raw/sel_example_F
	cp fakeFASTQ.fastq ./Enrich/M1H5-P5/data/raw/unsel_example_F
	python Enrich/enrich.py --mode=run_all --config_file=/Users/jimstapleton/Dropbox/Whitehead_lab/Maynard/Enrich/M1H5-P5/input/example_local_config	
