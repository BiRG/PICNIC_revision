#!usr/bin/perl

use strict;

my $hmmRootDir=$ARGV[0];

if ((-d $hmmRootDir )==1 && (-w $hmmRootDir)==0) {
	print "Directory $hmmRootDir exists but is not writable. Please check permissions\n";
}


if ((-d $hmmRootDir )==0) {
   mkdir($hmmRootDir) ;
} else {
	my $dataDir=$hmmRootDir."/data";
	my $dataDir1=$dataDir."/cancer";
	my $dataDir2=$dataDir1."/normalised";
	my $dataDir3=$dataDir1."/raw";
	my $outputDir=$hmmRootDir."/output";
	my $output2Dir=$hmmRootDir."/output2";
	my $configDir=$hmmRootDir."/config";

	if ((-d $dataDir )==0) {
   		mkdir($dataDir) ;
   		mkdir($dataDir1) ;
   		mkdir($dataDir2) ;
   		mkdir($dataDir3) ;
	}

	if ((-d $outputDir )==0) {
   		mkdir($outputDir) ;
	}

	if ((-d $output2Dir )==0) {
   		mkdir($output2Dir) ;
	}

	if ((-d $configDir )==0) {
   		mkdir($configDir) ;
	}
}

exit;


