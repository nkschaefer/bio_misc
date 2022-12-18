SHELL=/bin/bash
COMP=g++
FLAGS=-std=c++11 --std=gnu++11
LDFLAGS=-lz -lhts
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CELLAR=$(shell brew info argp-standalone | grep Cellar | cut -d' ' -f1)
	FLAGS += -I$(CELLAR)/include/
	LDFLAGS += -L$(CELLAR)/lib/ -largp
endif

all: bin/bcf2eigenstrat bin/bcf2treemix bin/bam_dummy_rg bin/eig_upgma bin/eig_dstat bin/filter_pairs bin/sort_huge_bed bin/bam_fq_pairs bin/split_read_file
MAXHAPS ?= 500

bin/bcf2eigenstrat: src/bcf2eigenstrat.cpp
	$(COMP) $(FLAGS)  src/bcf2eigenstrat.cpp -o bin/bcf2eigenstrat $(LDFLAGS)

bin/bcf2treemix: src/bcf2treemix.cpp
	$(COMP) $(FLAGS) src/bcf2treemix.cpp -o bin/bcf2treemix $(LDFLAGS)

bin/bam_dummy_rg: src/bam_dummy_rg.cpp src/bam.h bam.o
	$(COMP) $(FLAGS) src/bam_dummy_rg.cpp -o bin/bam_dummy_rg bam.o $(LDFLAGS)

bin/eig_upgma: src/eig_upgma.cpp treeNode.o
	$(COMP) -D MAXHAPS=$(MAXHAPS) $(FLAGS)  src/eig_upgma.cpp -o bin/eig_upgma treeNode.o $(LDFLAGS)

bin/eig_dstat: src/eig_dstat.cpp
	$(COMP) $(FLAGS) src/eig_dstat.cpp -o bin/eig_dstat $(LDFLAGS)

bin/filter_pairs: src/filter_pairs.cpp
	$(COMP) $(FLAGS) src/filter_pairs.cpp -o bin/filter_pairs $(LDFLAGS)

bin/split_read_file: src/split_read_file.cpp
	$(COMP) $(FLAGS) src/split_read_file.cpp -o bin/split_read_file $(LDFLAGS)

bin/sort_huge_bed: src/sort_huge_bed.cpp
	$(COMP) $(FLAGS) src/sort_huge_bed.cpp -o bin/sort_huge_bed $(LDFLAGS)

bin/bam_fq_pairs: src/bam_fq_pairs.cpp src/bam.h bam.o
	$(COMP) $(FLAGS) src/bam_fq_pairs.cpp -o bin/bam_fq_pairs bam.o $(LDFLAGS)

bam.o: src/bam.cpp src/bam.h
	$(COMP) $(FLAGS) -c src/bam.cpp $(LDFLAGS)

treeNode.o: src/treeNode.cpp
	$(COMP) -D MAXHAPS=$(MAXHAPS) $(FLAGS) $(OPTS) -c src/treeNode.cpp

clean:
	rm *.o
