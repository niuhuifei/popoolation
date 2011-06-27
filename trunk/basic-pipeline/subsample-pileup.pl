use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/../Modules";
use Pileup;
use Math::Round;

our $verbose=1;

my $input="";
my $output="";
my $help=0;
my $test=0;
my $minqual=0;
my $maxcoverage=0;
my $targetcoverage=0;
my $encoding="illumina";
my $mode="random"; #random/fraction

#--gtf /Users/robertkofler/testfiles/3R-3entries.gtf --input /Users/robertkofler/testfiles/3R.pileup --output /Users/robertkofler/testfiles/output/3R-filtered.pileup

GetOptions(
    "input=s"           =>\$input,
    "output=s"          =>\$output,
    "test"              =>\$test,
    "min-qual=i"        =>\$minqual,
    "max-coverage=i"    =>\$maxcoverage,
    "target-coverage=i" =>\$targetcoverage,
    "mode=s"            =>\$mode,
    "fastq-type"        =>\$encoding,
    "help"              =>\$help
) or die "Invalid arguments";

pod2usage(-verbose=>2) if $help;
VarTest::runTests() if $test;
pod2usage(-msg=>"Could not find pileup file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;
pod2usage(-msg=>"Please provide a valid target coverage (>0)",-verbose=>1) unless $targetcoverage;
pod2usage(-msg=>"Please provide a valid maximum coverage (>0)",-verbose=>1) unless $maxcoverage;


my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using minimum quality\t$minqual\n";
print $pfh "Using maximum coverage\t$maxcoverage\n";
print $pfh "Using subsample modus\t$mode\n";
print $pfh "Using fastq-type\t$encoding\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my $pp=get_pileup_parser($encoding,1,1,$maxcoverage,$minqual);
my $qualc=Utility::get_quality_char($encoding,$minqual);
my $subs=Utility::get_subsampler($mode,$targetcoverage);


open my $ifh,"<",$input or die "Could not open input file";
open my $ofh, ">",$output or die "Could not open output file";


while(my $line=<$ifh>)
{
    chomp $line;
    
    # $p=pos, chr, refc, totcov, eucov, A, T, C, G, del, N, iscov, issnp, ispuresnp
    my $p=$pp->($line);
    
    # check for coverage mincoverage, maxcoverage
    next unless $p->{iscov};
    
    # subsample
    my $newnucs=$subs->($p);
    next unless $newnucs;
    my $newqual=$qualc x length($newnucs);
    
    my @toprint=($p->{chr},$p->{pos},$p->{refc},length($newnucs),$newnucs,$newqual);
    print $ofh join("\t",@toprint)."\n";
}
close $ifh;
close $ofh;


exit;

{
    package Utility;
    use strict;
    use warnings;
    use Math::Round;
    
    sub get_subsampler
    {
        my $mode=shift;
        my $targetcoverage=shift;
    
        my $freqcalculator;
        if($mode eq "random")
        {
            $freqcalculator = \&_get_random_sample;
        }
        elsif($mode eq "fraction")
        {
            $freqcalculator = \&_get_fraction_sample;
        }
        else
        {
            die die "unknown mode $mode";
        }
        
        
        return sub
        {
            my $p=shift;

            my $samplecov=$p->{totcov}<$targetcoverage?$p->{totcov}:$targetcoverage;
            # calculate the new frequencies using the preset frequency calculator (random sampling method)
            my $toret=$freqcalculator->($samplecov,$p->{A},$p->{T},$p->{C},$p->{G},$p->{del},$p->{N});
            return $toret;
        }
        
    }
    
    sub _get_random_sample
    {
        my($samplecov,$a,$t,$c,$g,$del,$n)=@_;
        my $temp="A"x$a."T"x$t."C"x$c."G"x$g."*"x$del."N"x$n;
        my @ar=split //,$temp;
        
        my @novel=();
        for my $i (1..$samplecov)
        {
            my $leng=@ar;
            last unless $leng;
            my $index=int(rand()*$leng);
            my $element=splice(@ar,$index,1);
            push @novel,$element;
        }
        
        my $toret=join("",@novel);
        return $toret;
    }
    
    
    sub _get_fraction_sample
    {
        my($samplecov,$a,$t,$c,$g,$del,$n)=@_;
        my $cov=$a+$t+$c+$g+$del+$n;
        my($af,$tf,$cf,$gf,$delf,$nf)=($a/$cov,$t/$cov,$c/$cov,$g/$cov,$del/$cov,$n/$cov);
        
        my $an=round($af*$samplecov);
        my $tn=round($tf*$samplecov);
        my $cn=round($cf*$samplecov);
        my $gn=round($gf*$samplecov);
        my $deln=round($delf*$samplecov);
        my $nn=round($nf*$samplecov);
        my $str="A"x$an."T"x$tn."C"x$cn."G"x$gn."*"x$deln."N"x$nn;
        
        while(length($str)<$samplecov)
        {
            $str.="N";
        }
        while(length($str)>$samplecov)
        {
            my $index=int(rand() * length($str));
            substr($str,$index,1,"");
        }
        return $str; 
    }


    
    
    sub get_quality_char
    {
        my $encoding=shift;
        my $minqual=shift;
        my $qualhash={
            "illumina"=>    sub{chr(shift(@_)+64)},
            "sanger"        =>sub{chr(shift(@_)+33)}
        };
        
        my $conv=$qualhash->{$encoding};
        die "No quality converter for $encoding" unless $conv;
        return $conv->($minqual);
    }
    
    
}

    


    #"input=s"           =>\$input,
    #"output=s"          =>\$output,
    #"min-qual=i"        =>\$minqual,
    #"max-coverage=i"    =>\$maxcoverage,
    #"target-coverage=i" =>\$targetcoverage,
    #"mode=s"            =>\$mode,
    #"fastq-type"        =>\$encoding,
    #"help"              =>\$help

=head1 NAME

perl subsample-pileup.pl - Reduce the coverage of a pileup file to the given target coverage

=head1 SYNOPSIS

perl subsample-pileup.pl --input input.pileup --output output.pileup --target-coverage 50 --max-coverage 400 --min-qual 20 --mode random --fastq-type sanger

=head1 OPTIONS

=over 4

=item B<--input>

The input file in the pileup format; Mandatory.

For more information use <--help>

=item B<--output>

The output file, will be a pileup file  Mandatory.

=item B<--target-coverage>

Reduce the coverage of the pileup-file to the here provided value; If the coverage is smaller than this the coverage will not be reduced; Mandatory

=item B<--max-coverage>

The maximum coverage; Entries having a larger coverage than this will be ignored (without subsampling); Mandatory

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

=item B<--min-qual>

The minimum quality; Bases in the pileup having a lower quality than this will be ignored.

=item B<--mode>

Specify the method for reducing the coverage;
Either random subsampling without replacement (random) may be used
or the exact fraction of the allele frequency rounded to the next integer (fraciton);
random|fraction; default=random

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input pileup

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

The output will be a coverage reduced pileup.

=head2 Technical Details

The script proceeds in the following ways. First the pileup is filtered for quality while only retaining having a quality higher or equal than the minimum quality.
Second pileup position having a coverage higher than the maximum coverage are entirely discarded. In the third step the coverage is reduced to the targetcoverage either by random sampling without replacement or by calculating the exact fractions.
Note that the base quality of the output will be the minimum-quality for all bases!

  
=cut



