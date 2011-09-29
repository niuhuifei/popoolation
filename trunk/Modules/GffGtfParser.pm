#!/usr/bin/perl
package GffGtfParser;

use strict;
use warnings;
use Data::Dumper;
use Pileup;
use VarianceExactCorrection;
use VarianceUncorrected;

#export
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(
	get_characteristics_of_genome_gff_pileup
	print_variance_for_feature
);

our @EXPORT_OK = qw(
	_load_GFF_file_with_chromosome_arms
	_correct_chromosome_name
	_update_chromosome_arms_list
	_check_chromosome_arms
	_load
	_invert_feature_hash
	_calculate_characteristics_GFF
);


our %featHash = (
	"chromosome_arm"=>"*",
	"intron"	=>	"i",
	"exon"	=>	"e",
	"ncRNA"	=>	"1",
	"tRNA"	=>	"2",
	"snoRNA"	=>	"3",
	"snRNA"	=>	"4",
	"rRNA"	=>	"5",
	"CDS"	=>	"C",
	"five_prime_UTR"	=> "F",
	"three_prime_UTR"	=>	"T",
	"enhancer"	=>	"h",
	"miRNA"	=>	"m",
	"regulatory_region"	=>	"r",
	"pseudogene"	=>	"p",
	"transposable_element"	=> "t",
	"pre_miRNA"	=> "a"
);

sub get_characteristics_of_genome_gff_pileup{
	my ($GFF_FILE, $PILEUP_FILE,
	$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,
	$POOL_SIZE, $MEASURE, $UNCORRECTED) = @_;
	
	print "#data generated by script $0\n".
		  "#input gff file: $GFF_FILE\n".
		  "#input pileup file: $PILEUP_FILE\n".
		  "#qual-encoding: $QUAL_ENCODING, ".
		  "min-count: $MIN_COUNT, min-cov: $MIN_COV, max-cov: $MAX_COV, min-qual: $MIN_QUAL, pool-size: $POOL_SIZE, measure: $MEASURE ".
		  "uncorrected: $UNCORRECTED\n";
		  
	my $ptrAnnotation = {};
	$ptrAnnotation = _load_GFF_file_with_chromosome_arms($GFF_FILE);
				
	my $ptrGenomeCharacteristics = _calculate_characteristics_GFF(
			$ptrAnnotation,
			$PILEUP_FILE,
			$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,
			$POOL_SIZE, $MEASURE, $UNCORRECTED);
			
	return $ptrGenomeCharacteristics;
}

sub _load_GFF_file_with_chromosome_arms{
	my ($IN_FILE) = @_;
	
	my $ptrAnnotation={};
	my $ptrGff=[];
	
	($ptrGff) = _load($IN_FILE, \%featHash);
	#$ptrAnnotation, $ptrGff, $ptrFeatures	
	_overwrite_and_add_feat_hash($ptrAnnotation, $ptrGff, undef);
	
	return $ptrAnnotation;
}

sub _correct_chromosome_name{
	my ($chr)=@_;
	if (substr($chr, 0, 3) eq "chr"){
		return substr($chr,3);
	}else{
		return $chr;
	}
}

sub _update_chromosome_arms_list{
	my ($ptrChromosomeArms, $chr, $feat)=@_;

	#testing for chromosome arms
	if (exists($ptrChromosomeArms->{$chr})){
		if($feat eq "chromosome_arm"){
			$ptrChromosomeArms->{$chr}=1;
		}
	}else{
		if($feat eq "chromosome_arm"){
			$ptrChromosomeArms->{$chr}=1;
		}else{
			$ptrChromosomeArms->{$chr}=0;
		}
	}
}

sub _check_chromosome_arms{
	my ($ptrChromosomeArms, $ptrChromosomeArmsMissing)=@_;

	my $all_chromosome_arms_ok=1;
		
	foreach my $chromosome (keys %{$ptrChromosomeArms}){
		next unless ($ptrChromosomeArms->{$chromosome} == 0);
		$all_chromosome_arms_ok=0;
		push @$ptrChromosomeArmsMissing, $chromosome;
	}
		
	return $all_chromosome_arms_ok;	
}


sub _load{
	my($IN_FILE, $ptrFeatures) = @_;
	#$ptrGeneIDs -- list of gene names in the loaded GTF file
	#$ignore_intergenic -- bool that means: skip checking for chromosome_arm definition in GFF/GTF input file
	
	my $ptrGffData={};

	my $ptrWithdrawnInGtf={};
	
	open GFFdatabase, "<", $IN_FILE or die "Could not open gff/gtf file $IN_FILE";
	
	#read data line by line	
	while (my $line = <GFFdatabase>){
		chomp($line);
		
		#ignore comments;
		next if ($line=~m/^#/);
		
		my @tmp = split "\t", $line;
		#check a number of fields in the line;
		#die if the number is lower than 8 -- incorrect lines e.g. DNA sequence of mitochondrial genome in r5.32
		if (scalar(@tmp)<8){die "incorrect line in gff file $IN_FILE";}
		
		#split fields into different variables -- distinguish two possibilities (GFF or GTF input file) 
		my ($chr, $source, $feat, $start, $end, $score, $strand, $offset) = split "\t", $line;
		
		#ignore irrelevant features
		next unless(exists($ptrFeatures->{$feat}));	

		#create temporal record in hash -- it will delete redundant records
		my $record;
		$chr = _correct_chromosome_name($chr);
				
		$record = "$chr\t$source\t$feat\t$start\t$end\t$score\t$strand\t$offset";
		$ptrGffData->{$record} = 1;
		
	}
	close GFFdatabase;
	
	#hash of chromosome arms, needed for checking definition of each chromosome arm
	my $ptrChromosomeArms={};
	
	#from $ptrGffData create a structure with info splitted into different variables
	my $ptrGff=[];
	foreach my $record (keys %$ptrGffData){
		my ($chr,$source,$feat,$start,$end,$score,$strand,$offset) = split "\t", $record;
		
		my $ptrHash =	{
			chromosome=>$chr,
			feat=>$feat,
			start=>$start,
			end=>$end,
			score=>$score,
			strand=>$strand,
			offset=>$offset
		}
		push @$ptrGff, $ptrHash;
		_update_chromosome_arms_list($ptrChromosomeArms, $chr, $feat);		
	}
		
	#check chromosome arms	
	my $ptrChromosomeArmsMissing=[];
	my $all_chromosome_arms_ok = _check_chromosome_arms($ignore_intergenic, $ptrChromosomeArms, $ptrChromosomeArmsMissing);
	
	if (!$all_chromosome_arms_ok){
		die "Chromosome arms misssing for chromosomes @{$ptrChromosomeArmsMissing}. Add them and run again.";
	}	

	return ($ptrGff, $ptrWithdrawnInGtf);
}

sub _invert_feature_hash{
	my $ptrInverse={};
	
	foreach my $feature (keys %featHash){
		my $code = $featHash{$feature};
		if ($feature eq "chromosome_arm"){
			$ptrInverse->{$code} = "intergenic";
		}else{
			$ptrInverse->{$code} = $feature;
		}
	}
	return $rInverse;
}

sub _calculate_characteristics_GFF{
	my ($ptrAnnotation,
		$PILEUP_FILE,
		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,
		$POOL_SIZE, $MEASURE, $UNCORRECTED) = @_;
				
	my $ptrInverseFeatHash = _invert_feature_hash();
	
	my $pileupParser = get_pileup_parser($QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL);
	my $vec; 
	
	if ($UNCORRECTED){
		$vec=VarianceUncorrected->new($POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV);
	}else{
		$vec=VarianceExactCorrection->new($POOL_SIZE,$MIN_COUNT, $MIN_COV, $MAX_COV);
	}
	
	# count features lengths for each feature from annotation hash.
	my $ptrLength=_compute_feature_lengths_GFF($ptrAnnotation, $ptrInverseFeatHash);
	
	my $ptrCoveredLength={};
#	my $rSum_pi={};
#	my $rSum_theta={};
#	my $rSum_D={};
			
	open PILEUPfile, "<", $PILEUP_FILE or die "pileup file problem";
	
	while (my $line = <PILEUPfile>){
		chomp($line);
		my $parsedLine = $pileupParser->($line);
		
		# if this position is covered update couner for each feature of the position
		my $isCovered = $parsedLine->{iscov};
		next unless $isCovered;
		
		my $chromosome = $parsedLine->{chr};
		my $position = $parsedLine->{pos};
		my $featuresString = $ptrAnnotation->{$chromosome}[$position]{feat};
		
		foreach my $code (keys %$rInverseFeatHash){
			next unless defined($featuresString);
			next unless $featuresString =~m/\Q$code\E/;
		
			my $feature = $ptrInverseFeatHash->{$code};
			$ptrCoveredLength->{$feature}+=1;
		}
	
	
	# if this position is a pure SNP then store the SNP information into an array of SNPs for the particular measure that is going to be calculated
		my $isPureSNP = $parsedLine->{ispuresnp};
		next unless $isPureSNP;
		
		my $ptrPiSNPs = {};
		my $ptrThetaSNPs = {};
		my $ptrDSNPs = {};
		
		if ($MEASURE eq "pi"){
			push @{$ptrPiSNPs->{$feature}}, $parsedLine;
		}elsif($MEASURE eq "theta"){
			push @{$ptrThetaSNPs->{$feature}}, $parsedLine;
		}elsif($MEASURE eq "d"){
			push @{$ptrDSNPs->{$feature}}, $parsedLine;
		}else{
			
			push @{$ptrPiSNPs->{$feature}}, $parsedLine;
			push @{$ptrThetaSNPs->{$feature}}, $parsedLine;
			push @{$ptrDSNPs->{$feature}}, $parsedLine;						
			
		}		
	}
	close PILEUPfile;
	
	my $ptrPi = {};
	my $ptrTheta = {};
	my $ptrD = {};

	foreach my $feature (keys %$ptrCoveredLength){
		
		if ($MEASURE eq "pi"){
			$ptrPi->{$feature} = $vec->calculate_measure("pi", $ptrPiSNPs->{$feature}, $ptrCoveredLength->{$feature});
		}elsif($MEASURE eq "theta"){
			$ptrTheta->{$feature} = $vec->calculate_measure("theta", $ptrThetaSNPs->{$feature}, $ptrCoveredLength->{$feature});	
		}elsif($MEASURE eq "d"){
			$ptrD->{$feature} = $vec->calculate_measure("d", $ptrDSNPs->{$feature}, $ptrCoveredLength->{$feature});	
		}else{
			$ptrPi->{$feature} = $vec->calculate_measure("pi", $ptrPiSNPs->{$feature}, $ptrCoveredLength->{$feature});			
			$ptrTheta->{$feature} = $vec->calculate_measure("theta", $ptrThetaSNPs->{$feature}, $ptrCoveredLength->{$feature});	
			$ptrD->{$feature} = $vec->calculate_measure("d", $ptrDSNPs->{$feature}, $ptrCoveredLength->{$feature});				
		}
	}	

	
	my $ptrOut={};
	
	foreach my $code (keys %$ptrInverseFeatHash){
		my $feature = $ptrInverseFeatHash->{$code};

		if (defined($ptrLength->{$feature})){
			$ptrOut->{$feature}{totalLength} = $ptrLength->{$feature};
		}else{
			$ptrOut->{$feature}{totalLength} = 0;
		}

		if (defined($ptrCoveredLength->{$feature})and($ptrCoveredLength->{$feature}!=0)){
			$ptrOut->{$feautre}{coveredLength} = $ptrCoveredLength->{$feature};
			
			if ($MEASURE eq "pi"){
				$ptrOut->{$feature}{pi} = $ptrPi->{$feature};	
			}elsif($MEASURE eq "theta"){
				$ptrOut->{$feature}{theta} = $ptrTheta->{$feature};		
			}elsif($MEASURE eq "d"){
				$ptrOut->{$feature}{d} = $ptrD->{$feature};			
			}else{
				$ptrOut->{$feature}{pi} = $ptrPi->{$feature};	
				$ptrOut->{$feature}{theta} = $ptrTheta->{$feature};					
				$ptrOut->{$feature}{d} = $ptrD->{$feature};							
			}
			
			
		};
	}	
	return $ptrOut;
}

sub print_variance_for_feature{
	my ($ptrGenomeCharacteristics, $outFileHandle, $MEASURE, $scriptName, $paramsFile)=@_;
	
	my $ptrInverseFeatHash = _invert_feature_hash();
	print $outFileHandle "#data generated by a script $scriptName, for more details about parameters see $paramsFile\n";

	if ($MEASURE eq "pi"){
		print $outFileHandle "#feature\ttotal_length\tcovered_length\tpi\n";
	}elsif($MEASURE eq "theta"){
		print $outFileHandle "#feature\ttotal_length\tcovered_length\ttheta\n";
	}elsif($MEASURE eq "d"){
		print $outFileHandle "#feature\ttotal_length\tcovered_length\tTajD\n";
	}else{
		print $outFileHandle "#feature\ttotal_length\tcovered_length\tpi\ttheta\tTajD\n";
	}
	


	if ($MEASURE eq "pi"){

		foreach my $code (keys %$rInverseFeatHash){
			my $feature = $ptrInverseFeatHash->{$code};
			print "$feature\t".
				 "$ptrGenomeCharacteristics->{$feature}{totalLength}\t".
				 "$ptrGenomeCharacteristics->{$feature}{coveredLength}\t".
				 "$ptrGenomeCharacteristics->{$feature}{pi}\n";
		}


	}elsif($MEASURE eq "theta"){
		
		foreach my $code (keys %$rInverseFeatHash){
			my $feature = $ptrInverseFeatHash->{$code};
			print "$feature\t".
				 "$ptrGenomeCharacteristics->{$feature}{totalLength}\t".
				 "$ptrGenomeCharacteristics->{$feature}{coveredLength}\t".
				 "$ptrGenomeCharacteristics->{$feature}{theta}\n";
		}
		
	}elsif($MEASURE eq "d"){

		foreach my $code (keys %$rInverseFeatHash){
			my $feature = $ptrInverseFeatHash->{$code};
			print "$feature\t".
				 "$ptrGenomeCharacteristics->{$feature}{totalLength}\t".
				 "$ptrGenomeCharacteristics->{$feature}{coveredLength}\t".
				 "$ptrGenomeCharacteristics->{$feature}{d}\n";
		}

	}else{

		foreach my $code (keys %$rInverseFeatHash){
			my $feature = $ptrInverseFeatHash->{$code};
			print "$feature\t".
				 "$ptrGenomeCharacteristics->{$feature}{totalLength}\t".
				 "$ptrGenomeCharacteristics->{$feature}{coveredLength}\t".
				 "$ptrGenomeCharacteristics->{$feature}{pi}\t".
				 "$ptrGenomeCharacteristics->{$feature}{theta}\t".
				 "$ptrGenomeCharacteristics->{$feature}{d}\n";
		}
		
	}
}