HotSpot3D for AlphaFold DB
===========
This branch is under development by members of the Ding Lab at WashU. Code is being modified to handle structures from the AlphaFold Database instead of PDB.

This 3D proximity tool can be used to identify mutation hotspots from linear protein sequence and correlate the hotspots with known or potentially interacting domains, mutations, or drugs. Mutation-mutation and mutation-drug clusters can also be identified and viewed.

Usage
-----

                     Program:     HotSpot3D for AlphaFold DB - 3D mutation proximity analysis program.

                      Stable:     v1.8.2 

                        Beta:     alphafold_implementation v0.1 (this branch)
        
            Original Authors:     Beifang Niu, John Wallis, Adam D Scott, Sohini Sengupta, & Amila Weerasinghe
	    
	    AlphaFold DB version authors:	  Fernanda Martins Rodrigues

  Usage: hotspot3d <command> [options]

           Preprocessing
             drugport  --  0) Parse drugport database (OPTIONAL)
             uppro     --  1) Update proximity files
             prep      --  2) Run preprocessing steps 2a-2f
                 calroi    --  2a) Generate region of interest (ROI) information
                 statis    --  2b) Calculate p_values for pairs of mutations
                 anno      --  2c) Add region of interest (ROI) annotation
                 trans     --  2d) Add transcript annotation
                 cosmic    --  2e) Add COSMIC annotation to proximity file
                 prior     --  2f) Prioritization

           Analysis
		     main      --  Run analysis steps a-f (beta)
                 search    --  a) 3D mutation proximity searching
                 cluster   --  b) Determine mutation-mutation and mutation-drug clusters
                 sigclus   --  c) Determine significance of clusters (BETA/OPTIONAL)
                 summary   --  d) Summarize clusters (OPTIONAL)
                 visual    --  e) Visulization of 3D proximity (OPTIONAL)

Support
-------

For user support please email ruiyang.liu@wustl.edu & bniu@sccas.cn

For HotSpot3D for AlphaFold DB version please email Fernanda Martins Rodrigues (mrodrigues.fernanda@gmail.com; fernanda@wustl.edu)


Update
------

To reinstall code of the same version (in some cases, may need --sudo):

	cpanm --reinstall HotSpot3D-#.tar.gz


Install (Ubuntu 14.04.01)
-------

Make sure that you have cpanm:

	cpan App::cpanminus

For configuration, we recommend using local::lib:

	cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)

Dependencies include the modules: LWP::Simple, Test::Most, List::Util, List::MoreUtils, Parallel::ForkManager

	cpanm LWP::Simple

	cpanm Test::Most

	cpanm List::Util

	cpanm List::MoreUtils

	cpanm Parallel::ForkManager

Install HotSpot3D package: 
        
	git clone https://github.com/ding-lab/hotspot3d

	cd hotspot3d

For the latest stable version:

	git checkout v0.6.0
    
	cpanm HotSpot3D-0.6.0.tar.gz

For the latest beta version:

	git pull origin v1.8.2
    
	cpanm HotSpot3D-1.8.2.tar.gz

For HotSpot3D for AlphaFold DB version:
	
	git pull origin alphaFold_implementation

Final note: Installations under some organizations may use an internal perl version. 
To make use of the /usr/ perl, edit the first line of ~/perl5/bin/hotspot3d.
    
	from: #!/org/bin/perl
    
	to: #!/usr/bin/perl


We have developed webservice for HotSpot3D, and you can submit your mutation file to do
analysis online here: http://niulab.scgrid.cn/HotSpot3D/ , once you don't want to install
HotSpot3D locally. This does not apply for HotSpot3D for AlphaFold DB version.

Configure Environment
---------------------

	It is helpful to add your perl5 lib directory, and to add your perl5 bin directory.

	You can add the following lines to your ~/.bash_profile. Then run 'source ~/.bash_profile'.

		export PERL5LIB=~/perl5/lib/perl5/:${PERL5LIB}

		export PERL5BIN=~/perl5/bin/:${PERL5BIN}

		export PATH=~/perl5/bin/:${PATH}

	Add cosmic v67 information to 3D proximity results :

		mkdir preprocessing_dir/cosmic

		cp COSMIC/cosmic_67_for_HotSpot3D_missense_only.tsv.bz2 ./preprocessing_dir/cosmic/

		cd ./preprocessing_dir/cosmic/ 

		bzip2 -d cosmic_67_for_HotSpot3D_missense_only.tsv.bz2


Example - Preprocessing
-----------------------

## Download from Synapse
1. Go to https://www.synapse.org/#!Synapse:syn8699796, and check out the wiki for any updates/details.

2. Select the Files tab, then go into the AverageResidueDistance data directory (syn8717211).

3. The DrugPort processing results are located here (syn9704835) for those interested.

4. Select the reference/assembly version of interest (GRCh37 with Ensembl version 74 (syn9701918), or GRCh38 with Ensembl version 87 (syn9704851)).

5. You will need to download the hugo.uniprot.pdb.transcript.csv (syn9704852).

6. Two download options are available, prioritization.tar.gz (syn9704853) contains all human proteins that have been preprocessed. This is a large file that can take an hour or more depending on internet speeds. Alternatively, you can download the prioritization/ (syn9705109) or any specific protein proximity files within. The proximity files are compressed for faster/more targeted downloading. 

NOTE: Proximity data only contains pairs within 20Angstroms between mutations. This should be sufficient for many HotSpot3D applications.

## Generate on your own
1. (Optional) Run drugport module to parse Drugport data and generate a drugport parsing results flat file :

		hotspot3d drugport --pdb-file-dir=pdb_files_dir

2. Run 3D proximity calculation that also updates any existing preprocessed data (default launches LSF jobs) :

		hotspot3d uppro --output-dir=preprocessing_dir --pdb-file-dir=pdb_files_dir --drugport-file=drugport_parsing_results_file 1>hotspot3d.uppro.err 2>hotspot3d.uppro.out

3. Run automated preprocessing for other measurments and annotations (can alternatively run steps 2a-2f individually) :

		hotspot3d prep --output-dir=preprocessing_dir


Example - Analysis
------------------

3D proximity searching based on prioritization results and visualization

1. Proximity searching (acquire proximity information for input mutations):

		hotspot3d search --maf-file=your.maf --prep-dir=preprocessing_dir

2. Cluster pairwise data:

		hotspot3d cluster --pairwise-file=3D_Proximity.pairwise --maf-file=your.maf

3. Cluster significance calculation:

		hotspot3d sigclus --prep-dir=preprocessing_dir --pairwise-file=3D_Proximity.pairwise --clusters-file=3D_Proximity.pairwise.singleprotein.collapsed.clusters

4. Clustering Summary:

		hotspot3d summary --clusters-file=3D_Proximity.pairwise.singleprotein.collapsed.clusters

5. Visualization (works with PyMol):

		hotspot3d visual --pairwise-file=3D_Proximity.pairwise --clusters-file=3D_Proximity.pairwise.singleprotein.collapsed.clusters --pdb=3XSR

Annotations
-----------

Check out scripts/ for various annotation scripts to add more details to the .clusters file.

HGNC download can be found here: http://www.genenames.org/cgi-bin/genefamilies/.

Information on the Ensembl .gtf can be found here: http://useast.ensembl.org/info/website/upload/gff.html, and downloads can be found at the Ensembl ftp site, ftp://ftp.ensembl.org/pub/.

See the scripts/README.annotations for more details.

Tips
----

Mutation file - Standard .maf with custom coding transcript and protein annotations (ENST00000275493 and p.L858R)

There are only a handful of columns necessary from .maf files. They are:

		Hugo_Symbol
		
		Chromosome
		
		Start_Position
		
		End_Position
		
		Variant_Classification
		
		Reference_Allele
		
		Tumor_Seq_Allele1
		
		Tumor_Seq_Allele2
		
		Tumor_Sample_Barcode

And two non-standard columns:

		a transcript ID column
		
		a protein peptide change column (HGVS p. single letter abbreviations, ie p.T790M)

Current Annotation Support:

		Transcript ID - Ensembl coding transcript ID's (ENST)

		Gene name - HUGO symbol

Clustering with different pairs data:

		For monomers, you need to include the option '--meric-type monomer'

		For homomers, you need to include the option '--meric-type homomer'

		For heteromers, you need to include the option '--meric-type heteromer'

		For both homomers & heteromers simultaneously, you need to include the option '--meric-type multimer'

		For no regard to *mer status, you can include the option 
		'--meric-type unspecified', although this is run by default without the option

		For DrugPort only, do not input the .pairwise file; input only DrugPort pairs file.

		For *mer+DrugPort include the .pairwise file with the DrugPort pairs file, 
		and include the appropriate --meric-type as described above.

Clustering based on different distance measures:

        There are some pairs found on multiple structures. 
		In HotSpot3D versions v0.6.2 and earlier, 
		clustering only used the shortest distance among different structures 
		(shortest structure distance, SSD). 
		In HotSpot3D versions v0.6.3 and later, 
		clustering can be done using the average distance among different structures 
		(average structure distance, ASD), and this is now default.

Citation
--------

If you use HotSpot3D in your research, please cite:

* Protein-structure-guided discovery of functional mutations across 19 cancer types; Niu B, Scott AD, Sengupta S, Bailey MH, Batra P, Ning J, Wyczalkowski MA, Liang WW, Zhang Q, McLellan MD, Sun SQ, Tripathi P, Lou C, Ye K, Mashl RJ, Wallis J, Wendl MC, Chen F, Ding L; Nat Genet 2016 Aug;48(8):827-37
