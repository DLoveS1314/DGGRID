#!/bin/bash

dggridExe=../../build/src/apps/dggrid/dggrid
# f="aigenerate"
f="gridgenPureKML"
fmeta="gridgenPureKML"
# fmeta="transformseq2geo"
##f和=之间不要有空格 sh会自动识别
echo "#### running example out=>${fmeta}.txt" 
cd $f
$dggridExe ${fmeta}.meta >& outputfiles/${fmeta}.txt


