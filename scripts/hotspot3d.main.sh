#!/bin/bash
# Author: Adam D Scott (adamscott@wustl.edu)
# Wrapper for all Main steps of HotSpot3D except for visual module.

function run {
	maf=$1
	dataDir=$2
	prefix=$3

	pairwiseFile="${prefix}.pairwise"
	intraFile="${prefix}.pairwise.singleprotein.collapsed"
	interFile="${prefix}.pairwise.complex.collapsed"
	clustersIntraFile="${intraFile}.clusters"
	clustersInterFile="${interFile}.clusters"
	summaryIntraFile="${clustersIntraFile}.summary"
	summaryInterFile="${clustersInterFile}.summary"

	echo "hotspot3d search --maf-file ${maf} --data-dir ${dataDir} --output-prefix ${prefix} 1> ${prefix}.search.out 2> ${prefix}.search.err"
	hotspot3d search --maf-file ${maf} --data-dir ${dataDir} --output-prefix ${prefix} 1> ${prefix}.search.out 2> ${prefix}.search.err

	echo "hotspot3d post --maf-file ${maf} --input-prefix ${prefix} 1> ${prefix}.post.out 2>${prefix}.post.err"
	hotspot3d post --maf-file ${maf} --input-prefix ${prefix} 1> ${prefix}.post.out 2>${prefix}.post.err
	
	echo "hotspot3d cluster --collapsed-pairs-file ${intraFile} --pairwise-file ${pairwiseFile} --output-file ${clustersIntraFile} --maf-file ${maf} 1> ${prefix}.cluster.out 2> ${prefix}.cluster.err"
	hotspot3d cluster --collapsed-pairs-file ${intraFile} --pairwise-file ${pairwiseFile} --output-file ${clustersIntraFile} --maf-file ${maf} 1> ${prefix}.cluster.out 2> ${prefix}.cluster.err
	
	echo "hotspot3d summary --clusters-file ${clustersIntraFile} --output-file ${summaryIntraFile} 1> ${prefix}.summary.out 2> ${prefix}.summary.err"
	hotspot3d summary --clusters-file ${clustersIntraFile} --output-file ${summaryIntraFile} 1> ${prefix}.summary.out 2> ${prefix}.summary.err
}

defaultOut="hotspot3d.results"
defaultDir="./"
if [ ! -z $1 ]; then
	if [ ! -z $2 ]; then
		if [ ! -z $3 ]; then
			run $1 $2 $3
		else
			run $1 $2 ${defaultOut}
		fi
	else
		run $1 ${defaultDir} ${defaultOut}
	fi
else
	echo "bash $0 <maf-file> /preprocessing-dir/ \"output-prefix\""
fi
