HotSpot3D
===========

This 3D proximity tool can be used to identify mutation hotspots from linear protein sequence and correlate the hotspots with known or potentially interacting domains, mutations, or drugs. Mutation-mutation and mutation-drug clusters can also be identified and viewed.

Usage
-----

        Program:     HotSpot3D - 3D mutation proximity analysis program.
        Version:     V0.4
         Author:     Beifang Niu, John Wallis, Adam D Scott, & Sohini Sengupta

          Usage:     hotspot3d <command> [options]

Key commands:

        search    --  3D mutation proximity searching
		post      --  Post-processing on pairwise data
        cluster   --  Determine clusters from HotSpot3D inter, intra, and druggable pairwise data 
		sigclus   --  Determine significance of clusters
		summary   --  Determine cluster-level measures
        visual    --  Visulization of 3D proximity

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
For user support please email adamscott@wustl.edu


Install (Ubuntu 14.04.01)
-------

Prerequisites:

In order to install HotSpot3D package, first install CPANM
(cpanm - get, unpack build and install modules from CPANM)
NOTE: Some steps may require adding --force to install successfully.

        sudo apt-get install cpanminus

Another way to install cpanminus is to just download it, as per the installer
        
        curl -LO http://xrl.us/cpanm
        chmod +x cpanm

Or by using cpan

		cpan App::cpanminus

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
        cpanm HotSpot3D-#.#.tar.gz

(Installations under some organizations may use an internal perl version.
To make use of the /usr/ perl, edit the first line of ~/perl5/bin/hotspot3d.
from: #!/org/bin/perl
to: #!/usr/bin/perl)

Example
-------

Preprocessing procedure

1. (Optional) Run drugport module to parse Drugport data and generate a drugport parsing results flat file :

        hotspot3d drugport --pdb-file-dir=pdb_files_dir --output-file=drugport_parsing_results_file

2. Run 3D proximity calculation ( this step needs one output directory to store all of the data from the pre-processing procedure
and one directory which contains pdb files downloaded from the PDB website. This step will automatically download PDB files if
the necessary PDB file is not yet in the pdb files directory) :

        hotspot3d uppro --output-dir=preprocessing_output --pdb-file-dir=pdb_files_dir --drugport-file=drugport_parsing_results_file --max-3d-dis=100 1>hotspot3d.preprocessing.t.err 2>hotspot3d.preprocessing.t.out

3. Calculate protein domain information for each UniProt ID : 

        hotspot3d calroi --output-dir=preprocessing_output

4. Significance determination calculation :  

        hotspot3d statis --output-dir=preprocessing_output

5. Add protein domain annotation information to 3D proximity information :

        hotspot3d anno --output-dir=preprocessing_output

6. Choose transcripts based on the alignment between Uniprot sequence and human peptides sequences :

        hotspot3d trans --output-dir=preprocessing_output

7. Add cosmic v67 information to 3D proximity results :

        mkdir preprocessing_output/cosmic
        cp COSMIc/cosmic_67_for_HotSpot3D_missense_only.tsv.bz2 ./preprocessing_output/cosmic/
        cd ./preprocessing_output/cosmic/ 
        bzip2 -d cosmic_67_for_HotSpot3D_missense_only.tsv.bz2
        hotspot3d cosmic --output-dir=preprocessing_output

8. Prioritization :

        hotspot3d prior --output-dir=preprocessing_output --p-value=0.1 --3d-dis=20 --linear-dis=0.5


3D proximity searching based on prioritization results and visualization

1. Proximity searching (acquire proximity information for input mutations):

        hotspot3d search --maf-file=pancan19_input.maf --data-dir=preprocessing_output --output-prefix=pancan19 --skip-silent 1>pancan19.t.out 2>pancan19.t.err

2. Post-processing of pairwise data (required for cluster step):

		hotspot3d post --maf-file=pancan19_input.maf --input-prefix=pancan19

3. Cluster pairwise data:

        hotspot3d cluster --inter-intra-proximity-file=interactions_file --data-location-file=location_data --output-file=clustering.out --target-nontarget-file=drug_data_file

4. Cluster significance calculation:

        hotspot3d sigclus --prep-dir=preprocessing_output --pairwise=pairwise_file --clusters=clusters_file --output=output_file

5. Clustering Summary:

        hotspot3d summary --clusters-file=cluster_file --output-file=output_summary

6. Visualization:

        hotspot3d visual --pymol-dir=/usr/bin/pymol --output-dir=pymol_out --pdb-dir=pdb_files_dir

Tips
----

Current Annotation Support
Transcript ID - Ensembl coding transcript ID's (ENST)
Gene name - HUGO symbol
Mutation file - Standard .maf with custom coding and protein annotations (ENST00000275493 and p.L858R)

Clustering with different pairs data:
For intra you need to include the singleprotein pairs without DrugPort results/pairs.
For inter you need complex pairs without DrugPort pairs.
For DrugPort only, do not include singleprotein or complex pairs; include only DrugPort pairs.
For intra+inter you can concatenate the singleprotein and complex pairs without DrugPort pairs.
For intra+DrugPort include singleprotein pairs and DrugPort pairs.
For inter+DrugPort include complex pairs and DrugPort pairs.
For intra+inter+DrugPort include a concatenated singleprotein and complex pairs file with the DrugPort pairs.

Note that if concatenating pairs files, you should take care with removing the second header that will appear in the middle of the file. The .pairwise file contains both intra and inter pairs, so it can be used when involving intra or inter pairs.
