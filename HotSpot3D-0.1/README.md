HotSpot3D
===========

This 3D proximity tool can be used to identify the mutation hotspots in the linear 1D sequence and correlates these hotspots with known or potential interacting domains based on both known intermolecular interactions and calculated proximity for potential intramolecular interactions.

Usage
-----

        Program:     HotSpot3D - 3D mutation proximity analysis program.
        Version:     V0.2
         Author:     Beifang Niu && John Wallis

        Usage:  hotspot3d <command> [options]

Key commands:

        search    --  3D mutation proximity searching
        visual    --  Visulization of 3D proximity
        cluster   --  Determine mutation clusters from HotSpot3D inter, intra, and druggable data 

        drugport  --  Parse drugport database 
        uppro     --  Update proximity files
        calpro    --  Calculate proximity file for one UniProt ID
        calroi    --  Generate region of interest (ROI) information
        statis    --  Calculate p_values for pairs of mutations
        anno      --  Add region of interest (ROI) annotation
        trans     --  Add transcript annotation 
        homo      --  Add homology PDB structures 
        cosmic    --  Add COSMIC annotation to proximity file
        prior     --  Prioritization
        help      --  this message

SUPPORT
For user support please mail bniu@genome.wustl.edu


Install (Ubuntu 14.04.01)
-------

Prerequisites:

In order to install HotSpot3D package, we need CPANM program
(cpanm - get, unpack build and install modules from CPANM)


        sudo apt-get install cpanminus

Intall LWP::Simple module

        sudo apt-get install libwww-perl

Install HotSpot3D package: 
        
        git clone https://github.com/ding-lab/hotspot3d
        cd hotspot3d
        cpanm --force HotSpot3D-0.1.tar.gz


example
-------





