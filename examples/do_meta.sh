#!/bin/bash

dggridExe=../../build/src/apps/dggrid/dggrid
# f="aigenerate"
f="genChildren" #文件夹所在路径
fmeta="genchildren" #meta名字
# fmeta="gridgenDiamond" #meta名字

 
# fmeta="transformseq2geo"
##f和=之间不要有空格 sh会自动识别
echo "#### running example out=>${fmeta}.txt" 
cd $f
$dggridExe ${fmeta}.meta >& outputfiles/${fmeta}.txt


