#!/bin/bash

BINDIR=`dirname $(readlink -f "$0")/bin`
WORKINGDIR=`pwd`
ref_untrimmed=$1
threads=${2 : 1}
echo 'BINDIR '$BINDIR

refrev=$ref_untrimmed'.refrev'

if [[ ! -r $refrev ]]
then
  $BINDIR/refRevComp $ref_untrimmed $refrev
fi

echo 'Reverse-complimented reference: ' $refrev

sa=$ref_untrimmed'.sa'
if [[ ! -r $sa || $sa -ot $refrev ]]
then
	$BINDIR/mksary $refrev $sa
fi

mul=$ref_untrimmed'.mul.wig'
mur=$ref_untrimmed'.mur.wig'
echo 'MUL: '$mul
echo 'MUR: '$mur
if [[ ! -r $mul || ! -r $mur || $mul -ot $sa || $mur -ot $sa ]]
then
  $BINDIR/minUniqueKmer -t$threads -p$ref_untrimmed $ref_untrimmed $sa
fi
