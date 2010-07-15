use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use Pileup;
use SynchronizeUtility;
our $verbose=1;

my $pileupfile="";
my $output="";
my $minCount=2;
my $fastqtype="illumina";
my $minQual=20;

my $help=0;
my $test=0;
my $minCoverage=4;
my $maxCoverage=1000000;
my $region="";

GetOptions(
    "input=s"           =>\$pileupfile,
    "output=s"          =>\$output,
    "fastq-type=s"      =>\$fastqtype,
    "min-count=i"       =>\$minCount,
    "min-qual=i"        =>\$minQual,
    "region=s"          =>\$region,
    "min-coverage=i"    =>\$minCoverage,
    "max-coverage=i"    =>\$maxCoverage,
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";

pod2usage(-verbose=>2) if $help;
VarTest::runTests() if $test;
pod2usage(-msg=>"Could not find pileup file",-verbose=>1) unless -e $pileupfile;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;
pod2usage(-msg=>"Min count not provided",-verbose=>1) unless $minCount;
pod2usage(-msg=>"Min quality not valid. Has to be between 0 and 40",-verbose=>1) if $minQual<0 || $minQual > 40;
pod2usage(-msg=>"The minimum coverage hast to be at least two times the minimum count",-verbose=>1) unless $minCoverage >= (2*$minCount);


# qualencoding,mincount,mincov,maxcov,minqual
my $pp=get_pileup_parser($fastqtype,$minCount,$minCoverage,$maxCoverage,$minQual);
open my $ifh, "<",$pileupfile or die "Could not open input file";
open my $ofh,">",$output or die "Could not open output file";

my $reg;
$reg=Utility::parse_region($region) if $region;

while(my $line=<$ifh>)
{
    chomp $line;
    if($reg)
    {
        my($chr,$pos)=split /\t/,$line;
        next unless $chr eq $reg->{chr};
        next if $pos < $reg->{start} or $pos > $reg->{end};
    }
    my $pp=$pp->($line);
    next unless $pp->{ispuresnp};
    print $ofh "$pp->{chr}\t$pp->{pos}\t$pp->{refc}\t".format_parsed_pileup($pp)."\n";
}

{
    package Utility;
    use strict;
    use warnings;
    sub parse_region
    {
        my $reg=shift;
        my ($ref,$start,$end)=$reg=~m/^([^:]+):(\d+)[-:](\d+)/;
        return {
            ref=>$ref,
            start=>$start,
            end=>$end
        };
    }
}

=head1 NAME

perl call-SNPs.pl - A script which identifies SNPs from a pileup file (pooled Population)

=head1 SYNOPSIS

perl call-SNPs.pl --input input.pileup --output output.file --min-count 2 --min-coverage 4 --min-quality 20


=head1 OPTIONS

=over 4

=item B<--input>

The input file in the pileup format. A pooled population sequenced and mapped to the reference genome. Finally the mapping results have been converted to sam output format.
Using the samtools the sam format can be easily converted into the pileup format.  Mandatory.

=item B<--output>

The output file.  Mandatory.

=item B<--fastq-type>
The encoding of the quality characters; Must either be 'sanger' or 'illumina'; 

 Using the notation suggested by Cock et al (2009) the following applies:
 'sanger'   = fastq-sanger: phred encoding; offset of 33
 'solexa'   = fastq-solexa: -> NOT SUPPORTED
 'illumina' = fastq-illumina: phred encoding: offset of 64
 
 See also:
 Cock et al (2009) The Sanger FASTQ file format for sequecnes with quality socres,
 and the Solexa/Illumina FASTQ variants; 

default=illumina
 
=item B<--min-count>

The minimum count of the minor allele. This is important for the identification of SNPs; default=2

=item B<--min-coverage>

The minimum coverage of a site. Sites with a lower coverage will not be considered (for SNP identification and coverage estimation); default=4

=item B<--max-coverage>

the maximum coverage; used for SNP identification, the coverage in ALL populations has to be lower or equal to this threshold, otherwise no SNP will be called. default=1000000

=item B<--min-qual>

The minimum quality; Alleles with a quality lower than this threshold will not be considered (for SNP identification and coverage estimation); default=20

=item B<--region>

A region in the genome for which SNPs should be identified; must be in the format chr:start-end; default=""

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

A pileup file as described here: http://samtools.sourceforge.net/pileup.shtml; example:

 2L	90131	N	11	AaAAAaaAaAA	[aUQ_a`^_\Z
 2L	90132	N	11	AaAAAaaAaAA	_bYQ_^aaT^b
 2L	90133	N	11	A$aAAAaaAaAA	_b[Xaaa__Ua
 2L	90134	N	10	tTTTttTtTT	_`aaa_a[aa
 2L	90135	N	10	a$TAAaaAaAA	aZ^a`ba`\_
 2L	90136	N	9	TTTttTtTT	`aaaaaWaa
 2L	90137	N	9	GGGggGgGG	``aaaaQaa
 2L	90138	N	9	T$TTttTtTT	[U\`a\T^_
 2L	90139	N	8	TTttTtTT	``aaaU_a
 2L	90140	N	9	CCccCcCC^FC	[aaba`aaa

=head2 Output

For example:

 2L      8202    A       75:0:10:1:0:0
 2L      8263    G       20:0:0:80:0:0

 col1: reference sequence
 col2: position in the reference sequence
 col3: reference character
 col4: lightwight representation of the SNP allele frequences in the form A:T:C:G:N:del
 A.. count of allele A (only those with min-quality >=threshold)
 T.. count of allele T
 C.. count of allele C
 G.. count of allele G
 N.. count of Ns
 del.. count of deletions (are deletions spanning this position)

=cut

