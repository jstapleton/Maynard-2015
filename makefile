
all: forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq reverse_unpaired1.fq selFASTQ.fastq forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq reverse_unpaired2.fq unselFASTQ.fastq enrichFiles


clean:
	-rm -f forward_paired.fq reverse_paired.fq forward_unpaired.fq reverse_unpaired.fq fakeFASTQ.fastq out.* water.txt read.fasta

# replace seleceted_F.fastq and selected_R.fastq with the names of the selected forward and reverse paired-end fastq files to be converted
forward_paired1.fq reverse_paired1.fq forward_unpaired1.fq reverse_unpaired1.fq:
	java -jar Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 20 selected_F.fastq selected_R.fastq forward_paired1.fq forward_unpaired1.fq reverse_paired1.fq reverse_unpaired1.fq ILLUMINACLIP:/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:6:1 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20

selFASTQ.fastq:
	python mutationCounter.py forward_paired1.fq forward_unpaired1.fq reverse_paired1.fq reverse_unpaired1.fq
	mv fakeFASTQ.fastq selFASTQ.fastq
	cp selFASTQ.fastq ./data/raw/sel_example_F
	rm -f out.hist out.histogram water.txt read.fasta

# replace unselected_F.fastq and unselected_R.fastq with the names of the unselected forward and reverse paired-end fastq files to be converted
forward_paired2.fq reverse_paired2.fq forward_unpaired2.fq reverse_unpaired2.fq:
	java -jar Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 20 unselected_F.fastq unselected_R.fastq forward_paired2.fq forward_unpaired2.fq reverse_paired2.fq reverse_unpaired2.fq ILLUMINACLIP:/Trimmomatic-0.32/adapters/TruSeq3-PE.fa:2:30:6:1 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:20

unselFASTQ.fastq:
	python mutationCounter.py forward_paired2.fq forward_unpaired2.fq reverse_paired2.fq reverse_unpaired2.fq
	mv fakeFASTQ.fastq unselFASTQ.fastq
	cp unselFASTQ.fastq ./data/raw/unsel_example_F
	rm -f out.hist out.histogram water.txt read.fasta

enrichFiles:
	python enrich.py --mode=run_all --config_file=/input/example_local_config
