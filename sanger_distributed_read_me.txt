*********************************************************
Supplied code and configuration files
*********************************************************

1) Directory seg contains:

This code takes a feature_intensity file (parsed .cel file). 

* .m files
	* segmentGenome.m -> main script to execute

2) Directory config. 
	*Contains config.mat which has all the coefficients calculated from training.  Needed by code in seg and hmm 
 
3) Directory hmm
	* Contains code to run the HMM. 

4) Directory celConverter
	* Contains the files necessary for the CelFileConverter (.jar file and supporting files).


********************************************************
Execution of PICNIC 
********************************************************
In it's current form, PICNIC can only be run using Matlab.  The execution requires the following software and specifications:

* Matlab: including the optimisation and statistics toolboxes.
* 2G memory 
* perl
* java

The standalone version of the algorithms will be released shortly.

Running PICNIC 
---------------

1. Combine .m files from 'seg' and 'hmm' to a common directory.  This will form the 'root directory' from which the rest of the code is executed.

2.Run the 'runHMM.pl' script. Type the following command:

	perl runHMM.pl 'rootDir' 
	
	where rootDir is the location of the matlab source files from the previous step.

3. Conversion of .cel file to text

	a) Place .cel files to analyse in a separate directory.  It is advisable to place only one .cel file per directory.

	b) Download a 'cdf' file from the Affymetrix website: http://www.affymetrix.com/Auth/support/downloads/library_files/genomewidesnp6_libraryfile.zip.

	c) Run the .cel file converter from directory celConverter (WARNING: you will ned 2G memory for this step!). Navigate to celDir on a command prompt.

	Run the following command:

		java -Xmx2G -jar CelFileConverter.jar -m Snp6FeatureMappings.csv -c 'cdf_file_including_path' -s 'directory name of cel files' -t rootDir/data/cancer/raw
 
	For further instructions on running the cel file converter type 
		java -jar CelFileConverter.jar -h
	
	The parsed .cel file will be put into rootDir/data/cancer/raw.
	
4. Run the normalisation/segmentation part of the algorithm.  The output of this step will go into folder rootDir/output. 
	
	From rootDir, type the command below in matlab.  
		run segmentGenome.m		 

5. Manual inspection of /output.
        open the 'genome_fig_*.fig' file. This displays an initial segmentation of the copy number intensities against genomic position for all chromosomes. Black segments indicate an even copy number. Red segments indicate LOH. Green segments indicate a copy number of at least three. Yellow indicates an unknown copy number. This chart should give an indication of how copy number intensity relates to actual copy number. For example, in 'Reference/genome_fig_sample.jpg' the segments cluster into three groups corresponding to copy number intensities of 0.75, 0.94 and 1.13 corresponding to copy numbers 2,3 and 4. The following parameters need to be selected for step 5:
	init_val: the copy number intensity of homozygous deletions. 
	delta_val: the change in intensity between successive copy numbers.

The parameters for the example given are init_val=0.37 and delta_val=0.19.

6. Run the HMM from rootDir.  The inputs to the function were the values determined in the previous steps.  The output of this step will go into folder rootDir/output2.  
	From rootDir, type the following command at the matlab prompt.
	 run 'HMMRun.m init_val delta_val' 

PICNIC Output
-------------

The /Output2 folder contains the following files.

A.	'fig_n_name.fig' files for chromosomes n=1,2,...,24.
B. 	A 'complexity_name.csv' file that counts state changes per chromsome.
C.	An 'LOH_name' contains the proportion of genome in LOH, split by chromosome.
D.	A 'params_name.mat' file containing the parameters used in the HMM.
E.	A 'genome_fig_name.fig' file used in step 4. above.
F.	A 'genes_name.csv' file with four columns indicating copy number intensity, copy number, probability of deletion and probability of change in copy number state. The rows correspond to the rows in the reference file 'Reference/GeneFootPrint.csv'.
G.	A 'qual_name.csv' file used to estimate the quality of the overal fit.
H.	A params file that contains the output for all the probes. The columns are specifically:

1	SNP Identifier
2	Raw Intensity Ratio
3	Allelic Angle
4	Actual Copy Number
5	Segmented Total Copy Number
6	Segmented Minor Copy Number
7	Middle Fitted Angle Height (above 0.5)
8	Outer Fitted Angle Height (above 0.5)
9	LOH index
10	No. A copies (genotype)
11	No. B copies (genotype)
12	State Change Probability
13	Genotyping Confidence
14	Genotyping Confidence Conditional Upon State Classification
15	Heterozygous Probability
16	Allele A LOH probability
17	Allele B LOH probability

The information relating to rows can be found in 'Reference/ProbeRef.csv'. Columns are respectively SNP identifier, chromosome, position and a SNP flag (1=SNP,0=non-polymorphic probe).


