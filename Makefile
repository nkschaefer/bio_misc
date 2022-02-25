SHELL=/bin/bash
COMP=g++

all: bin/bcf2eigenstrat bin/bcf2treemix bin/eig_upgma bin/eig_dstat bin/filter_pairs bin/sort_huge_bed bin/bam_fq_pairs bin/bam_split_snps bin/bam_subs_pipe bin/split_read_file bin/atac_fq_preprocess
MAXHAPS ?= 500

bin/bcf2eigenstrat: src/bcf2eigenstrat.cpp
	$(COMP) -std=c++11 --std=gnu++11  src/bcf2eigenstrat.cpp -o bin/bcf2eigenstrat -lhts -lz

bin/bcf2treemix: src/bcf2treemix.cpp
	$(COMP) -std=c++11 --std=gnu++11 src/bcf2treemix.cpp -o bin/bcf2treemix -lhts -lz

bin/eig_upgma: src/eig_upgma.cpp treeNode.o
	$(COMP) -D MAXHAPS=$(MAXHAPS) -std=c++11 --std=gnu++11  src/eig_upgma.cpp -o bin/eig_upgma treeNode.o -lz

bin/eig_dstat: src/eig_dstat.cpp
	$(COMP) -std=c++11 --std=gnu++11 src/eig_dstat.cpp -o bin/eig_dstat -lz

bin/filter_pairs: src/filter_pairs.cpp
	$(COMP) -std=c++11 --std=gnu++11 src/filter_pairs.cpp -o bin/filter_pairs -lhts -lz

bin/split_read_file: src/split_read_file.cpp
	$(COMP) -std=c++11 --std=gnu++11 src/split_read_file.cpp -o bin/split_read_file -lhts -lz

bin/sort_huge_bed: src/sort_huge_bed.cpp
	$(COMP) -std=c++11 --std=gnu++11 src/sort_huge_bed.cpp -o bin/sort_huge_bed -lz

bin/bam_fq_pairs: src/bam_fq_pairs.cpp src/bam.h bam.o
	$(COMP) -std=c++11 --std=gnu++11 src/bam_fq_pairs.cpp -o bin/bam_fq_pairs bam.o -lz -lhts

bin/bam_split_snps: src/bam_split_snps.cpp src/bam.h bam.o
	$(COMP) -std=c++11 --std=gnu++11 src/bam_split_snps.cpp -o bin/bam_split_snps bam.o -lz -lhts

bin/bam_subs_pipe: src/bam_subs_pipe.cpp src/bam.h bam.o
	$(COMP) -std=c++11 --std=gnu++11 src/bam_subs_pipe.cpp -o bin/bam_subs_pipe bam.o -lz -lhts

bam.o: src/bam.cpp src/bam.h
	$(COMP) -std=c++11 -c src/bam.cpp -lhts

treeNode.o: src/treeNode.cpp
	$(COMP) -D MAXHAPS=$(MAXHAPS) -std=c++11 --std=gnu++11 $(OPTS) -c src/treeNode.cpp

clean:
	rm *.o
