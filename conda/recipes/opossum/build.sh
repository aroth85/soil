#!/bin/bash
set -eu -o pipefail

outdir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM
mkdir -p $outdir
mkdir -p $PREFIX/bin
cp -r * $outdir
cp $RECIPE_DIR/opossum.sh $outdir/opossum
ln -s $outdir/opossum $PREFIX/bin
