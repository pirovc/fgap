# FGAP: an automated gap closing tool 

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/fgap/README.html)

Piro, V. C., Faoro, H., Weiss, V. a, Steffens, M. B., Pedrosa, F. O., Souza, E. M., & Raittz, R. T. (2014). FGAP: an automated gap closing tool. BMC Research Notes, 7(1), 371. http://doi.org/10.1186/1756-0500-7-371

Install and run (bioconda):
---------------------------

	conda install -c bioconda fgap
	FGAP "-d sample_data/DRAFT_ecoli_hiseq454.fasta -a 'sample_data/DATASET_ecoli_hiseq.fasta,sample_data/DATASET_ecoli_hiseq.fasta' -o sample_data_results -t 2"

Install and run (source):
-------------------------

	# On the same folder of the fgap.m script
	octave --no-gui --eval 'fgap -d sample_data/DRAFT_ecoli_hiseq454.fasta -a sample_data/DATASET_ecoli_hiseq.fasta'

Install and run (compiled version):
-----------------------------------

	# Download MCR
	https://sourceforge.net/projects/fgap/files/MCR_LINUX64b.tar.gz

	# Download compiled FGAP
	https://sourceforge.net/projects/fgap/files/FGAP_1_8_1_LINUX64b.tar.gz
	
	# Install MCR
	tar xf MCR_LINUX64b.tar.gz
	cd MCR_LINUX64b/
	./installMCR.sh /home/user/MCR/
	
	# Run FGAP
	tar xf FGAP_1_8_1_LINUX64b.tar.gz
	cd FGAP_1_8_1_LINUX64b/
	./run_fgap.sh /home/user/MCR/v717/ -d sample_data/DRAFT_ecoli_hiseq454.fasta -a sample_data/DATASET_ecoli_hiseq.fasta
	
Parameters:
-----------

	------------------------------------------
		        FGAP v1.8.1
	------------------------------------------


	Usage in command-line mode (compiled): ./run_fgap.sh <MCR installation folder> -d <draft file> -a "<dataset(s) file(s)>" [parameters]
	Usage in Matlab/Octave (source): fgap -d <draft file> -a '<dataset(s) file(s)>' [parameters]

	-d /--draft-file        Draft genome file [fasta format - Ex: 'draft.fasta']
	-a /--datasets-files    List of datasets files to close gaps [fasta format - Ex: 'dataset1.fasta,dataset2.fasta']

	-s /--min-score         Min Score (raw) to return results from BLAST (integer) - Default: 25
	-e /--max-evalue        Max E-Value to return results from BLAST (float) - Default: 1e-7
	-i /--min-identity      Min identity (%) to return results from BLAST (integer [0-100]) - Default: 70

	-C /--contig-end-length Length (bp) of contig ends to perform BLAST alignment (integer) - Default: 300
	-T /--edge-trim-length  Length of ignored bases (bp) upstream and downstrem of the gap (integer) - Default: 0
	-R /--max-remove-length Max number of bases (bp) that can be removed (integer) - Default: 500
	-I /--max-insert-length Max number of bases (bp) that can be inserted (integer) - Default: 500

	-p /--positive-gap      Enable closing of positive gaps (with insertion) (integer [0-1]) - Default: 1
	-z /--zero-gap          Enable closing of zero gaps (without insert any base) (integer [0-1]) - Default: 0
	-g /--negative-gap      Enable closing of negative gaps (overlapping contig ends) (integer [0-1]) - Default: 0

	-c /--gap-char                          Base that represents the gap (char) - Default: 'N'
	-b /--blast-path                        Blast+ package path (only makeblastdb and blastn are needed, version 2.2.28+ or higher) - Default: ''
	-l /--blast-alignment-parameters        BLAST alignment parameters (opengap,extendgap,match,mismatch,wordsize) - Default: '1,1,1,-3,15'
	-r /--blast-max-results                 Max results from BLAST for each query (integer) - Default: 200
	-t /--threads                           Number of threads (integer) - Default: 1

	-m /--more-output       More output files with gap regions after and before gap closing (integer [0-1]) - Default: 0
	-o /--output-prefix     Output prefix [File or folder - Ex: 'out' or 'out_folder/out' ] - Default: 'output_fgap'
	-h /--help              This help message
