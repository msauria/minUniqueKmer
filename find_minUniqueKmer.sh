#!/bin/bash
BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`
ref_untrimmed=$1
echo 'BINDIR '$BINDIR

refrev=$ref_untrimmed'.refrev'
if [[ ! -x $BINDIR/refRevComp ]]
then
  g++ $BINDIR/refRevComp.cpp -o $BINDIR/refRevComp --std=c++11
fi

if [[ ! -r $refrev ]]
then
  $BINDIR/refRevComp $ref_untrimmed $refrev
fi

echo 'Reverse-complimented reference: ' $refrev
if [[ ! -x $BINDIR/libdivsufsort/build/examples/mksary ]]
then
	cd $BINDIR
	git submodule update --init --recursive
	cd libdivsufsort
	mkdir -f build
	cd build
	cmake -DCMAKE_BUILD_TYPE="Release" \
	-DCMAKE_INSTALL_PREFIX="/usr/local" ..
	sed -i 's/int32_t/int64_t/g' include/divsufsort.h
	make
	cd $WORKINGDIR
fi

sa=$ref_untrimmed'.sa'
if [[ ! -r $sa || $sa -ot $refrev ]]
then
	$BINDIR/libdivsufsort/build/examples/mksary $refrev $sa
fi

mul=$ref_untrimmed'.mul.bg'
mur=$ref_untrimmed'.mur.bg'
if [[ ! -x $BINDIR/minUniqueKmer.cpp ]]
then
	g++ $BINDIR/minUniqueKmer.cpp -o $BINDIR/minUniqueKmer --std=c++11
fi

echo 'MUL: '$mul
echo 'MUR: '$mur
if [[ ! -r $mul || ! -r $mur || $mul -ot $sa || $mur -ot $sa ]]
then
  $BINDIR/minUniqueKmer $ref_untrimmed $sa $mul $mur
fi
