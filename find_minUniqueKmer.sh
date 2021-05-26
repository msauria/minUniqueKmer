BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`
echo 'BINDIR '$BINDIR

cd $BINDIR
git submodule update --init --recursive
cd libdivsufsort
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE="Release" \
-DCMAKE_INSTALL_PREFIX="/usr/local" ..
sed -i 's/int32_t/int64_t/g' include/divsufsort.h
make

cd $WORKINGDIR

ref_untrimmed=$1
refrev=$ref_untrimmed'.refrev'
if [ ! -r refrev ]
then
  g++ $BINDIR/refRevComp.cpp -o $BINDIR/refRevComp
  $BINDIR/refRevComp $ref_untrimmed $refrev
fi
sa=$ref_untrimmed'.sa'

echo 'Reverse-complimented reference: ' $refrev
$BINDIR/libdivsufsort/build/examples/mksary $refrev $sa
mul=$ref_untrimmed'.mul.bg'
mur=$ref_untrimmed'.mur.bg'

g++ $BINDIR/minUniqueKmer.cpp -o $BINDIR/minUniqueKmer
echo 'MUL: '$mul
echo 'MUR: '$mur
if [ ! -r $mul || -r $mur ]
then
  $BINDIR/minUniqueKmer $ref_untrimmed $sa $mul $mur
fi
