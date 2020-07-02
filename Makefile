SHELL=/bin/bash
COMP=g++

all: bin/bcf2eigenstrat bin/bcf2treemix bin/eig_upgma bin/eig_dstat bin/filter_pairs bin/sort_huge_bed
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

bin/sort_huge_bed: src/sort_huge_bed.cpp
	$(COMP) -std=c++11 --std=gnu++11 src/sort_huge_bed.cpp -o bin/sort_huge_bed -lz
	
treeNode.o: src/treeNode.cpp
	$(COMP) -D MAXHAPS=$(MAXHAPS) -std=c++11 --std=gnu++11 $(OPTS) -c src/treeNode.cpp
	
clean:
	rm *.o
