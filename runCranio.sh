#!/bin/bash

cd /shared/metanetworkSynapse/

qsub -v s3="s3://metanetworks/CRANIO/",dataFile="/shared/CRANIO/cranioRNAseq.csv",pathv="/shared/metanetworkSynapse/",sparrowZ=1,sparrow2Z=1,lassoCV1se=1,lassoCVmin=1,lassoAIC=1,lassoBIC=1,ridgeCV1se=1,ridgeCVmin=1,ridgeAIC=1,ridgeBIC=1,genie3=1,tigressRootN=1,elasticNetAIC=1,elasticNetBIC=1,elasticNetCVmin=1,elasticNetCV1se=1,numberCore=319,outputpath="/shared/CRANIO/" -pe orte 319 -S /bin/bash -V -cwd -N cranio -e /shared/CRANIO/error.txt -o /shared/CRANIO/out.txt buildNet.sh
