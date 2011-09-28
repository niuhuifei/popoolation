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
	load_GFF_file_with_chromosome_arms 
	load_GFF_file_without_chromosome_arms
	load_GTF_file_CDS_with_gene_IDs
	load_GTF_file_CDS_with_gene_IDs_skip_withdrawn_genes
	get_characteristics_of_genome_gff_pileup
);

our @EXPORT_OK = qw(
	_correct_chromosome_name
	 
	_check_chromosome_arms 
	_update_chromosome_arms_list 
	
	_load
	_overwrite_feature_whole_annotation 
	
	_add_feature_to_specified_area
	_add_feature_whole_annotation
	_is_in_featuresArray
	_add_features
	
	_overwrite_and_add_feat_hash 
	
	_invert_feature_hash 
	_compute_feature_lengths_GFF 
	_calculate_characteristics_GFF
	
	_add_gene_ID_to_annotation_only_CDS_GTF
	_add_gene_ID_to_hash_only_CDS_GTF
	_add_gene_IDs_only_CDS_GTF
	
	_load_withdrawn_genes_into_hash
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
	$POOL_SIZE,$MEASURE, $UNCORRECTED) = @_;
	

	
	print "#data generated by script 'PopGenTools/misc/variance-for-features.pl'\n".
		  "#input gff file: $GFF_FILE\n".
		  "#input pileup file: $PILEUP_FILE\n".
		  "#qual-encoding: $QUAL_ENCODING, ".
		  "min-count: $MIN_COUNT, min-cov: $MIN_COV, max-cov: $MAX_COV, min-qual: $MIN_QUAL, pool-size: $POOL_SIZE, ".
		  "uncorrected: $UNCORRECTED\n";
		  
	my $rAnnotation = {};
	
	$rAnnotation = load_GFF_file_with_chromosome_arms($GFF_FILE);
				
	my $rGenomeCharacteristics = _calculate_characteristics_GFF(
			$rAnnotation,
			$PILEUP_FILE,
			$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,
			$POOL_SIZE, $UNCORRECTED);
			
	return $rGenomeCharacteristics;
}


###########################					calculating characteristics of a genome
sub _calculate_characteristics_GFF{
	my ($rAnnotation,
		$PILEUP_FILE,
		$QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL,
		$POOL_SIZE, $UNCORRECTED) = @_;
			
	# 		
	my $rInverseFeatHash = _invert_feature_hash();
	my $pileupParser = get_pileup_parser($QUAL_ENCODING, $MIN_COUNT, $MIN_COV, $MAX_COV, $MIN_QUAL);
	
	my $vec; 
	
	if ($UNCORRECTED == 0){
		$vec=VarianceExactCorrection->new($POOL_SIZE,$MIN_COUNT, $MIN_COV, $MAX_COV);
	}else{
		$vec=VarianceUncorrected->new($POOL_SIZE, $MIN_COUNT, $MIN_COV, $MAX_COV);
	}
	
	# count features lengths for each feature from annotation hash.
	my $rLength=_compute_feature_lengths_GFF($rAnnotation, $rInverseFeatHash);
	
	my $rCoveredLength={};
	my $rSum_pi={};
	my $rSum_theta={};
	my $rSum_D={};
			
	open PILEUPfile, "<", $PILEUP_FILE or die "pileup file problem";
	
	while (my $line = <PILEUPfile>){
		chomp($line);
		my $parsedLine = $pileupParser->($line);
		
		# if this position is covered update couner for each feature of the position
		my $isCovered = $parsedLine->{iscov};
		next unless $isCovered;
		
		my $chromosome = $parsedLine->{chr};
		my $position = $parsedLine->{pos};
		my $featuresString = $rAnnotation->{$chromosome}[$position]{feat};
		
		foreach my $code (keys %$rInverseFeatHash){
			next unless defined($featuresString);
			next unless $featuresString =~m/\Q$code\E/;
		
			my $feature = $rInverseFeatHash->{$code};
			$rCoveredLength->{$feature}+=1;
		}
	
	# if this position is a pure SNP compute all three measures for the position and update $sum_measure
		my $isPureSNP = $parsedLine->{ispuresnp};
		next unless $isPureSNP;
		
		my $pi= $vec->calculate_measure("pi",[$parsedLine],1);
		my $theta =$vec->calculate_measure("theta",[$parsedLine],1);
		my $D = $vec->calculate_measure("d",[$parsedLine],1);
		foreach my $code (keys %$rInverseFeatHash){
			next unless defined($featuresString);
			next unless $featuresString =~m/\Q$code\E/;
			my $feature = $rInverseFeatHash->{$code};
			$rSum_pi->{$feature}+=$pi;
			$rSum_theta->{$feature}+=$theta;
			$rSum_D->{$feature}+=$D;
		}
	}
	
	my $rNumbers={};
	
	foreach my $code (keys %$rInverseFeatHash){
		my $feature = $rInverseFeatHash->{$code};
		
		if (defined($rLength->{$feature})){
			$rNumbers->{$feature}{totalLength} = $rLength->{$feature};
		}else{
			$rNumbers->{$feature}{totalLength} = 0;
		}
		
		if (defined($rCoveredLength->{$feature})and($rCoveredLength->{$feature}!=0)){
			$rNumbers->{$feature}{coveredLength} = $rCoveredLength->{$feature};
			$rNumbers->{$feature}{pi}=$rSum_pi->{$feature} / $rCoveredLength->{$feature};
			$rNumbers->{$feature}{theta}=$rSum_theta->{$feature} / $rCoveredLength->{$feature};
			$rNumbers->{$feature}{D}=$rSum_D->{$feature} / $rCoveredLength->{$feature};
		}else{
			$rNumbers->{$feature}{coveredLength} = 0;
			$rNumbers->{$feature}{pi}=0;
			$rNumbers->{$feature}{theta}=0;
			$rNumbers->{$feature}{D}=0;
		}
	}
	close PILEUPfile;
	
	print "#feature\ttotal_length\tcovered_length\tpi\ttheta\tD\n";
	foreach my $code (keys %$rInverseFeatHash){
		my $feature = $rInverseFeatHash->{$code};
		print "$feature\t".
			  "$rNumbers->{$feature}{totalLength}\t".
			  "$rNumbers->{$feature}{coveredLength}\t".
			  "$rNumbers->{$feature}{pi}\t".
			  "$rNumbers->{$feature}{theta}\t".
			  "$rNumbers->{$feature}{D}\n";
	}

	return $rNumbers;
}

sub _compute_feature_lengths_GFF{
	
	my ($rAnnotation, $rInverseFeatHash)=@_;
	
	my $rFeatureLengths = {};
	
	foreach my $code (keys %$rInverseFeatHash){
		my $feature = $rInverseFeatHash->{$code};
	
		$rFeatureLengths->{$feature}=0;
	
		foreach my $chromosome (keys %$rAnnotation){
			my $length_chr=0;
			my $end = @{$rAnnotation->{$chromosome}};
			
			for (my $i=0; $i<$end; $i++){
				my $featuresString = $rAnnotation->{$chromosome}[$i]{feat};
				next unless defined($featuresString);
				next unless $rAnnotation->{$chromosome}[$i]{feat}=~m/\Q$code\E/;
				$length_chr+=1;
			}
			$rFeatureLengths->{$feature}+=$length_chr;	
		}
			
	}
	return $rFeatureLengths;
}

sub _invert_feature_hash{
	my $rInverse={};
	
	foreach my $feature (keys %featHash){
		my $code = $featHash{$feature};
		if ($feature eq "chromosome_arm"){
			$rInverse->{$code} = "intergenic";
		}else{
			$rInverse->{$code} = $feature;
		}
	}
	return $rInverse;
}

################### 			loading GFF file and creating an annotation


#sub load_GTF_file_with_chromosome_arms{}

#sub load_GTF_file_without_chromosome_arms{}

sub load_GTF_file_CDS_with_gene_IDs_skip_withdrawn_genes{
	my ($IN_FILE, $WITHDRAWN_FILE,$ptrGeneIDs)=@_;
	
	my $ptrWithdrawn=_load_withdrawn_genes_into_hash($WITHDRAWN_FILE);
	
	my ($ptrAnnotation, $ptrWithdrawnInGtf) = load_GTF_file_CDS_with_gene_IDs($IN_FILE, $ptrGeneIDs, $ptrWithdrawn);
	
	#print "annotation:\n";
	#print Dumper($ptrAnnotation, $ptrWithdrawnInGtf, $IN_FILE);
	
	
	return ($ptrAnnotation, $ptrWithdrawnInGtf);	
}

sub _load_withdrawn_genes_into_hash{
	my ($WITHDRAWN_FILE)=@_;
	my $ptrWithdrawn={};
	
	open withdrawnHandle, "<", $WITHDRAWN_FILE or die "Could not open withdrawn list $WITHDRAWN_FILE";
	while (my $geneID = <withdrawnHandle>){
		chomp($geneID);	
		$ptrWithdrawn->{$geneID}=1;
	}
	close withdrawnHandle;
	
	return $ptrWithdrawn;
}

sub load_GTF_file_CDS_with_gene_IDs{
	my ($IN_FILE, $ptrGeneIDs, $ptrWithdrawn)=@_;
	
	my $ptrCDS = {"CDS"=>"C"};
	
	
	#$IN_FILE, $ptrFeatures, $ignore_intergenic, $ptrGeneIDs
	my $ignore_intergenic = 1;
	my $is_gtf=1;
	my ($ptrGtf, $ptrWithdrawnInGtf)=_load($IN_FILE,$ptrCDS, $ignore_intergenic, $is_gtf, $ptrWithdrawn);
	my $codon_bool = 1;
	my $ptrAnnotation={};
	$ptrAnnotation = _overwrite_feature_whole_annotation_with_codon_positions($ptrAnnotation, $ptrGtf,"CDS",$ptrCDS, $codon_bool);
	
	_add_gene_IDs_only_CDS_GTF($ptrAnnotation, $ptrGtf, $ptrGeneIDs);
	
	return ($ptrAnnotation, $ptrWithdrawnInGtf);	
}



sub load_GFF_file_without_chromosome_arms{
	my ($IN_FILE) = @_;
	
	my $ptrAnnotation={};
	my $ptrGff=[];
	
	my $ignore_intergenic = 1;
	
	($ptrGff) = _load($IN_FILE, \%featHash, $ignore_intergenic);
	
	_overwrite_and_add_feat_hash($ptrAnnotation, $ptrGff, undef);	
	
	return $ptrAnnotation;
}

sub load_GFF_file_with_chromosome_arms{
	my ($IN_FILE) = @_;
	
	my $ptrAnnotation={};
	my $ptrGff=[];
	
	($ptrGff) = _load($IN_FILE, \%featHash);
	#$ptrAnnotation, $ptrGff, $ptrFeatures	
	_overwrite_and_add_feat_hash($ptrAnnotation, $ptrGff, undef);
	
	return $ptrAnnotation;
}

sub _overwrite_and_add_feat_hash{
	my ($ptrAnnotation, $ptrGff, $ptrFeatures)=@_;

	_overwrite_feature_whole_annotation($ptrAnnotation, $ptrGff, "chromosome_arm",  $ptrFeatures);
	_overwrite_feature_whole_annotation($ptrAnnotation, $ptrGff, "intron",  $ptrFeatures);
	_overwrite_feature_whole_annotation($ptrAnnotation, $ptrGff, "exon", $ptrFeatures);
		
	_overwrite_feature_whole_annotation($ptrAnnotation, $ptrGff, "ncRNA", $ptrFeatures);
	_overwrite_feature_whole_annotation($ptrAnnotation, $ptrGff, "tRNA", $ptrFeatures);
	_overwrite_feature_whole_annotation($ptrAnnotation, $ptrGff, "snoRNA", $ptrFeatures);
	_overwrite_feature_whole_annotation($ptrAnnotation, $ptrGff, "snRNA", $ptrFeatures);
	_overwrite_feature_whole_annotation($ptrAnnotation, $ptrGff, "rRNA", $ptrFeatures);
	
	my $ptrAddFeaturesArray = [
		"CDS", 
		"five_prime_UTR", 
		"three_prime_UTR", 
		"enhancer", 
		"miRNA", 
		"regulatory_region", 
		"pseudogene", 
		"transposable_element", 
		"pre_miRNA"
	];
	_add_features($ptrAnnotation, $ptrGff, $ptrAddFeaturesArray, $ptrFeatures);
#	_add_feature($ptrAnnotation, $ptrGff, "five_prime_UTR", $ptrFeatures);
#	_add_feature($ptrAnnotation, $ptrGff, "three_prime_UTR", $ptrFeatures);
#	_add_feature($ptrAnnotation, $ptrGff, "enhancer", $ptrFeatures);
#	_add_feature($ptrAnnotation, $ptrGff, "miRNA", $ptrFeatures);
#	_add_feature($ptrAnnotation, $ptrGff, "regulatory_region", $ptrFeatures);
#	_add_feature($ptrAnnotation, $ptrGff, "pseudogene", $ptrFeatures);
#	_add_feature($ptrAnnotation, $ptrGff, "transposable_element", $ptrFeatures);
#	_add_feature($ptrAnnotation, $ptrGff, "pre_miRNA", $ptrFeatures);
	
	return $ptrAnnotation;
}

sub _load{
	#loads GFF/GTF file into annotation hash of arrays called $rGff -- reference to a hash 
	# hash->{chromosome}[position_on_chromosome]
	#$IN_FILE,$ptrCDS, $ignore_intergenic, $is_gtf, $ptrWithdrawn
	my($IN_FILE, $ptrFeatures, $ignore_intergenic, $is_gtf, $ptrWithdrawn) = @_;
	#$ptrGeneIDs -- list of gene names in the loaded GTF file
	#$ignore_intergenic -- bool that means: skip checking for chromosome_arm definition in GFF/GTF input file
	
	my $rGffData={};

	my $ptrWithdrawnInGtf={};
	
	open GFFdatabase, "<", $IN_FILE or die "Could not open gff/gtf file $IN_FILE";
	
	#read data line by line	
	while (my $line = <GFFdatabase>){
		chomp($line);
		
		#ignore comments;
		next if ($line=~m/^#/);
		
		#check a number of fields in the line;
		#ignore if the number is lower than 8 -- incorrect lines e.g. DNA sequence of mitochondrial genome in r5.32
		my @tmp = split "\t", $line;
		if (scalar(@tmp)<8){die "incorrect lines in gff/gtf file $IN_FILE";}
		
		#split fields into different variables -- distinguish two possibilities (GFF or GTF input file) 
		my $chr; my $source; my $feat; my $start; my $end; my $score; my $strand; my $offset; my $geneID; my $transcriptID; 

		if ( defined($is_gtf) and ($is_gtf) ){	
			my $IDs;
			($chr,$source,$feat,$start,$end,$score,$strand,$offset, $IDs)= @tmp;
			my ($gene, $transcript)= split ";", $IDs;
			my $g; my $t;

			($g, $geneID)= split / /, $gene;
			$geneID = substr($geneID,1,length($geneID)-2);
			
			$transcript = substr($transcript, 1);
			($t, $transcriptID) = split / /, $transcript;
			$transcriptID = substr($transcriptID, 1, length($transcriptID)-2);
		}else{
			($chr,$source,$feat,$start,$end,$score,$strand,$offset) = @tmp;
		}
		
		#ignore irrelevant features
		next unless(exists($ptrFeatures->{$feat}));	

		#ignore $geneID if is in withdrown hash
		my $geneGr="";	my $geneTr="";
		if ( defined($is_gtf) and ($is_gtf) ){	
			($geneGr,$geneTr)=split /-/, $geneID;
		}
		if (exists($ptrWithdrawn->{$geneGr})){
			$ptrWithdrawnInGtf->{$geneID}=1;
		}else{ 
			#create temporal record in hash -- it will delete redundant records
			#distinguish two possibilities (GFF or GTF input file)
			my $record;
			$chr = _correct_chromosome_name($chr);
			if ( defined($is_gtf) and ($is_gtf) ){
				$record = "$chr\t$source\t$feat\t$start\t$end\t$score\t$strand\t$offset\t$geneID\t$transcriptID";
			}else{
				$record = "$chr\t$source\t$feat\t$start\t$end\t$score\t$strand\t$offset";
			}
			$rGffData->{$record} = 1;
		};
	}
	close GFFdatabase;
	
#	print Dumper($rGffData);
	
	#hash of chromosome arms, needed for checking definition of eac chromosome arm
	my $ptrChromosomeArms={};
	
	#move unique splited records to a list
	my $ptrGff=[];
	foreach my $record (keys %$rGffData){
		my $chr; my $source; my $feat; my $start; my $end; my $score; my $strand; my $offset; my $geneID; my $transcriptID;
		my $ptrHash={};
		if ( defined($is_gtf) and ($is_gtf) ){
			($chr,$source,$feat,$start,$end,$score,$strand,$offset, $geneID, $transcriptID) = split "\t", $record;
			
			$ptrHash =	{
				chromosome=>$chr,
				feat=>$feat,
				start=>$start,
				end=>$end,
				score=>$score,
				strand=>$strand,
				offset=>$offset,
				geneID=>$geneID,
				transcriptID => $transcriptID
			};		
		}else{
			($chr,$source,$feat,$start,$end,$score,$strand,$offset) = split "\t", $record;
		
			$ptrHash =	{
				chromosome=>$chr,
				feat=>$feat,
				start=>$start,
				end=>$end,
				score=>$score,
				strand=>$strand,
				offset=>$offset
			}
		}
		push @$ptrGff, $ptrHash;
		_update_chromosome_arms_list($ptrChromosomeArms, $ignore_intergenic, $chr, $feat);		
	}
	
	#print Dumper($ptrGff);
		
	my $ptrChromosomeArmsMissing=[];
	my $all_chromosome_arms_ok = _check_chromosome_arms($ignore_intergenic, $ptrChromosomeArms, $ptrChromosomeArmsMissing);
	
	if (!$all_chromosome_arms_ok){
		die "Chromosome arms misssing for chromosomes @{$ptrChromosomeArmsMissing}. Add them and run again.";
	}	

	return ($ptrGff, $ptrWithdrawnInGtf);
}

sub _correct_chromosome_name{
	my ($chr)=@_;
	if (substr($chr, 0, 3) eq "chr"){
		return substr($chr,3);
	}else{
		return $chr;
	}
}

sub _check_chromosome_arms{
	my ($ignore_intergenic, $ptrChromosomeArms, $ptrChromosomeArmsMissing)=@_;

	my $all_chromosome_arms_ok=1;
		
	if(!(defined($ignore_intergenic)) or (!$ignore_intergenic)){
		foreach my $chromosome (keys %{$ptrChromosomeArms}){
			next unless ($ptrChromosomeArms->{$chromosome} == 0);
			$all_chromosome_arms_ok=0;
			push @$ptrChromosomeArmsMissing, $chromosome;
			#print Dumper($ptrChromosomeArmsMissing);
		}
		
	}
	return $all_chromosome_arms_ok;	
}

sub _update_chromosome_arms_list{
	my ($ptrChromosomeArms, $ignore_intergenic, $chr, $feat)=@_;

	if (!defined($ignore_intergenic) or (!$ignore_intergenic) ){
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

}

sub _overwrite_feature_whole_annotation{
	my ($rAnnotation,$rGffList,$feature, $ptrFeatures)=@_;
	
	for (my $j=0; $j<scalar @$rGffList; $j++){
		next unless ($feature eq $rGffList->[$j]{feat});
				
		my $feat;
		if (defined($ptrFeatures)){
			$feat=$ptrFeatures->{$feature};
		}else{
			$feat=$featHash{$feature};
		}		
		
		#print Dumper($feat);
		
		my $start=$rGffList->[$j]{start};
		my $end=$rGffList->[$j]{end};
		my $chromosome=$rGffList->[$j]{chromosome};
		my $strand = $rGffList->[$j]{strand}; 
		
		for (my $i=$start; $i<=$end; $i++){
			$rAnnotation->{$chromosome}[$i]{feat}= $feat;
			$rAnnotation->{$chromosome}[$i]{strand}= $strand;
		}				
	}
}



#
# later join with _overwrite_feature_whole_annotation
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sub _overwrite_feature_whole_annotation_with_codon_positions{
	my ($rAnnotation,$rGffList,$feature, $ptrFeatures, $codon_bool)=@_;
	#codon_bool makes sense only for CDS data
	if (!defined($codon_bool)){$codon_bool=0}
	
	for (my $j=0; $j<scalar @$rGffList; $j++){
		next unless ($feature eq $rGffList->[$j]{feat});
				
		my $feat;
		if (defined($ptrFeatures)){
			$feat=$ptrFeatures->{$feature};
		}else{
			$feat=$featHash{$feature};
		}		
		
		#print Dumper($feat);
		
		my $start=$rGffList->[$j]{start};
		my $end=$rGffList->[$j]{end};
		my $chromosome=$rGffList->[$j]{chromosome};
		my $strand = $rGffList->[$j]{strand}; 
		
		my $codon;
		my $offset;
		
		if ($codon_bool){$offset = $rGffList->[$j]{offset}}
		
		
		for (my $i=$start; $i<=$end; $i++){
			$rAnnotation->{$chromosome}[$i]{feat}= $feat;
			$rAnnotation->{$chromosome}[$i]{strand}= $strand;
			next unless ($codon_bool);
			 
			if ($strand eq "+"){
				if ($offset%3 == 0){
					$rAnnotation->{$chromosome}[$i]{codonPosition}= ($i - $start)%3;
				}elsif($offset%3 == 1){
					$rAnnotation->{$chromosome}[$i]{codonPosition}= ($i - $start-1)%3;
				}elsif($offset%3 == 2){
					$rAnnotation->{$chromosome}[$i]{codonPosition}= ($i - $start+1)%3;
				}else{print "no chance ;)\n"}
			}elsif($strand eq "-"){
				if ($offset%3 == 0){
					$rAnnotation->{$chromosome}[$i]{codonPosition}= ($end - $i)%3;
				}elsif($offset%3 == 1){
					$rAnnotation->{$chromosome}[$i]{codonPosition}= ($end - $i-1)%3;
				}elsif($offset%3 == 2){
					$rAnnotation->{$chromosome}[$i]{codonPosition}= ($end - $i+1)%3;
				}
			}	 
		}
						
	}
	return $rAnnotation;
}






sub _add_features{
	#add all features from $ptrFeaturesArray
	#use coding in $ptrFeatures or default
	my ($ptrAnnotation,$ptrGff, $ptrFeaturesArray,$ptrFeatures)=@_;
	
	for (my $j=0; $j<scalar @$ptrGff; $j++){
		my $f = $ptrGff->[$j]{feat};
		next unless ( _is_in_featuresArray($f, $ptrFeaturesArray) );
		
		my $start = $ptrGff->[$j]{start};
		my $end = $ptrGff->[$j]{end};
		my $chromosome = $ptrGff->[$j]{chromosome};
		
		my $featCode;
		if (defined($ptrFeatures)){
			$featCode = $ptrFeatures->{$f};
		}else{
			$featCode = $featHash{$f};
		}

		 _add_feature_to_specified_area($ptrAnnotation, $chromosome, $start, $end, $featCode);
	}
}


sub _add_feature_to_specified_area{
	my ($ptrAnnotation, $chromosome, $start, $end, $featCode)=@_;
	
	for (my $i = $start; $i<=$end; $i++){
		next unless defined($ptrAnnotation->{$chromosome}[$i]);
		next if $ptrAnnotation->{$chromosome}[$i]{feat}=~m/\Q$featCode\E/;			
		$ptrAnnotation->{$chromosome}[$i]{feat} = $ptrAnnotation->{$chromosome}[$i]{feat}.$featCode;
	}
	
	
}

sub _is_in_featuresArray{
	my ($f, $ptrFeaturesArray)=@_;
	
	my $isInArray=0;

	foreach my $feature (@{$ptrFeaturesArray}){
		next unless ($feature eq $f);
		$isInArray=1;
	}

	return $isInArray;
}

sub _add_feature_whole_annotation{
	my ($ptrAnnotation,$rGffList,$feature,$ptrFeatures)=@_;

	for (my $j=0; $j<scalar @$rGffList; $j++){
		next unless ($feature eq $rGffList->[$j]{feat});

		my $feat;
		if (defined($ptrFeatures)){
			$feat=$ptrFeatures->{$feature};
		}else{
			$feat=$featHash{$feature};
		}
		
		my $start=$rGffList->[$j]{start};
		my $end=$rGffList->[$j]{end};
		my $chromosome=$rGffList->[$j]{chromosome};
		
		_add_feature_to_specified_area($ptrAnnotation, $chromosome, $start, $end, $feat);
	}
}

sub _add_gene_ID_to_hash_only_CDS_GTF{
	my ($ptrGeneIDs, $geneID, $chromosome, $start, $end, $strand, $feature)=@_;
	
	my $ptrHash={
		chromosome=>$chromosome, 
		start=>$start,
		end=>$end,
		strand =>$strand,
		feature=>$feature
	};	
	push @{$ptrGeneIDs->{$geneID}}, $ptrHash;
}


sub _add_gene_ID_to_annotation_only_CDS_GTF{
	my ($ptrAnnotation, $geneID, $chromosome, $start, $end)=@_;

	for (my $i=$start; $i<=$end; $i++){	
		push @{$ptrAnnotation->{$chromosome}[$i]{geneID}}, $geneID; 
	}
}

sub _add_gene_IDs_only_CDS_GTF{
	my ($ptrAnnotation, $ptrGtf, $ptrGeneIDs)= @_;	
	
	for (my $j=0; $j<scalar @$ptrGtf; $j++){
		next unless ($ptrGtf->[$j]{feat} eq "CDS");

		my $start=$ptrGtf->[$j]{start};
		my $end=$ptrGtf->[$j]{end};
		my $chromosome=$ptrGtf->[$j]{chromosome};
		my $geneID = $ptrGtf->[$j]{geneID};
		my $strand = $ptrGtf->[$j]{strand};
		my $feature = $ptrGtf->[$j]{feat};
		
		_add_gene_ID_to_annotation_only_CDS_GTF($ptrAnnotation, $geneID, $chromosome, $start, $end);
		_add_gene_ID_to_hash_only_CDS_GTF($ptrGeneIDs, $geneID, $chromosome, $start, $end, $strand, $feature);	
	}
}

1;