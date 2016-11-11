Density Scripts
===========

These scripts could be used to do density based clustering and visualization.

Usage
-----

Usage: perl DensityAll.pl pairwise_file_name_in_./Test/ Epsilon MinPts PDB_ID

Results will be copied to ./Results/ directory.

The DensityAll.pl script runs three different scripts in the following order:

1) OpticsWithR.pl : OPTICS clustering - outputs an ordered list of variants with corresponding reachability distances(Output file name: RD.$Epsilon.$MinPts.pairwise_file_name).

2) SuperClustersID.pl : Performs clustering. Whenever an important event like merging two clusters or appearing new cluster happens, records clusters with IDs.
 
(Output file name1: RD.$Epsilon.$MinPts.pairwise_file_name.SuperClustersID.clusters, 
Output file name2: RD.$Epsilon.$MinPts.pairwise_file_name.SuperClustersID.plot, 
Output file name3: SuperClustersID.RD.$Epsilon.$MinPts.pairwise_file_name.pdf ,
Output file name4: RD.$Epsilon.$MinPts.pairwise_file_name.clusters.shiny.R). 

3) DensityVisual.pl : Writes a pymol script for visualization

4) ClusterProbability.pl : Determines the membership probability of each variant in clusters.

Additional Softwares
-----------------------------

Install RStudio and the R package "shiny" for better visualization of the reachability plot.


