#!/bin/bash
BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`
ref_untrimmed=$1
echo 'BINDIR '$BINDIR

refbis=$ref_untrimmed'.refbismap'
if [[ ! -x $BINDIR/refBismap ]]
then
  g++ $BINDIR/refBismap.cpp -o $BINDIR/refBismap --std=c++11
fi

if [[ ! -r $refbis ]]
then
  $BINDIR/refBismap $ref_untrimmed $refbis
fi

echo 'Bisulfite-treated reference: ' $refbis
if [[ ! -x $BINDIR/libdivsufsort/build/examples/mksary ]]
then
	cd $BINDIR
	git submodule update --init --recursive
	cd libdivsufsort
	mkdir -p build
	cd build
	cmake -DCMAKE_BUILD_TYPE="Release" \
	-DCMAKE_INSTALL_PREFIX="/usr/local" ..
	sed -i 's/int32_t/int64_t/g' include/divsufsort.h
	make
	cd $WORKINGDIR
fi

sa=$ref_untrimmed'.sa'
if [[ ! -r $sa || $sa -ot $refbis ]]
then
	$BINDIR/libdivsufsort/build/examples/mksary $refbis $sa
fi

mul=$ref_untrimmed'.bis.mul.wig'
mur=$ref_untrimmed'.bis.mur.wig'
if [[ ! -x $BINDIR/minUniqueBismapKmer ]]
then
	g++ $BINDIR/minUniqueBismapKmer.cpp -o $BINDIR/minUniqueBismapKmer --std=c++11
fi

echo 'MUL: '$mul
echo 'MUR: '$mur
if [[ ! -r $mul || ! -r $mur || $mul -ot $sa || $mur -ot $sa ]]
then
  $BINDIR/minUniqueBismapKmer $ref_untrimmed $sa $mul $mur
fi
