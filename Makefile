all: refRevComp refBismap mksary minUniqueKmer meanKmerCoverage kmerUniqueMapping

CC 	= g++
CFLAGS	= -pthread -std=c++11

refRevComp: refRevComp.cpp
	-mkdir -p bin
	$(CC) $(CFLAGS) -o bin/$@ $^

refBismap: refBismap.cpp
	-mkdir -p bin
	$(CC) $(CFLAGS) -o bin/$@ $^

minUniqueKmer: minUniqueKmer.cpp CTPL/ctpl_stl.h
	-mkdir -p bin
	$(CC) $(CFLAGS) -o bin/$@ $< -ICTPL

meanKmerCoverage: meanKmerCoverage.cpp CTPL/ctpl_stl.h
	-mkdir -p bin
	$(CC) $(CFLAGS) -o bin/$@ $< -ICTPL

kmerUniqueMapping: kmerUniqueMapping.cpp
	-mkdir -p bin
	$(CC) $(CFLAGS) -o bin/$@ $<

CTPL/ctpl_stl.h:
	git submodule update --init --recursive

mksary: bin/mksary
	touch $<

bin/mksary:
	-mkdir -p bin
	git submodule update --init --recursive
	cd libdivsufsort && \
	cmake -S . -B build
	sed -i -e 's/int32_t/int64_t/g' libdivsufsort/build/include/divsufsort.h
	cd libdivsufsort/build && \
	make
	cp libdivsufsort/build/examples/mksary $@

clean:
	rm -r bin