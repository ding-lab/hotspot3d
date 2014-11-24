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

Intall Perl5 local lib

        cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)

Intall LWP::Simple module

        sudo apt-get install libwww-perl

Intall Test::Most module
        
        wget http://search.cpan.org/CPAN/authors/id/O/OV/OVID/Test-Most-0.34.tar.gz
        cpanm Test-Most-0.34.tar.gz

Install HotSpot3D package: 
        
        git clone https://github.com/ding-lab/hotspot3d
        cd hotspot3d
        cpanm HotSpot3D-0.1.tar.gz


example
-------

Preprocessing procedure


1. Run drugport module to parse Drugport data to generate a drugport parsing results flat file :

        hotspot3d drugport --pdb-file-dir=pdb_files_dir --output-file=drugport_parsing_results_file


2. Run 3D proximity calculation ( this step need one output directory to store all of the data from pre-processing procedure
and one directory which contains pdb files downloaded from PDB website. This step will automatically download PDB files if
there is no needed PDB files in pdb files directory) :

        hotspot3d uppro --output-dir=preprocessing_output --pdb-file-dir=pdb_files_dir --drugport-file=drugport_parsing_results_file --max-3d-dis=100 1>hotspot3d.preprocessing.t.err 2>hotspot3d.preprocessing.t.out

3. Calculate protein domain information for each Uniprot ID : 

        hotspot3d calroi --output-dir=preprocessing_output

4. Significance determination calculation :  

        hotspot3d statis --output-dir=preprocessing_output

5. Add protein domain annotation information to 3D proximity information :

        hotspot3d anno --output-dir=preprocessing_output

6. Choose transcripts based on the alignment between Uniprot sequence and human peptides sequences :

        hotspot3d trans --output-dir=preprocessing_output

7. Add cosmic v67 information to 3D proximity results :

        mkdir preprocessing_output/cosmic
        cp cosmic_67_for_HotSpot3D_missense_only.tsv ./preprocessing_output/cosmic/
        hotspot3d cosmic --output-dir=preprocessing_output

8. Prioritization :

        hotspot3d prior --output-dir=preprocessing_output --p-value=0.1 --3d-dis=20 --linear-dis=0.5


3D proximity searching based on prioritization results and visualization

1. Proximity searching :

        hotspot3d search --maf-file=pancan19_input.maf --data-dir=preprocessing_output --output-prefix=pancan19 --skip-silent 1>pancan19.t.out 2>pancan19.t.err

2. Visualization :

        hotspot3d visual --pymol-dir=/usr/bin/pymol --output-dir=pymol_out --pdb-dir=pdb_files_dir

3. 3D proximity results clustering : 

        hotspot3d cluster --inter-intra-proximity-file=interactions_file --data-location-file=location_data --output-file=clustering.out --target-nontarget-file=drug_data_file


