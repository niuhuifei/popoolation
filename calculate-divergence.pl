#!/usr/bin/perl-w
use strict;
use warnings;
use Getopt::Long; # to get options
use Pod::Usage; # to show documentations
use File::Basename; # to get the file path, file name and file extension

# Author: Ram Vinay Pandey

# Modified on 12-04-2010 to implement uniform Window sliding module.

# Define the variables
my $input;
my $output="";
my $help=0;
my $test=0;
my $verbose=1;

my $windowsize=1000;
my $step=100;
my $minCoverageFraction=0.6;


my $usage="perl $0 --input mauve-parsing-putput.txt --output dxy-output.txt --window-size 1000 --step-size 100\n";

GetOptions(
    "input=s"	    =>\$input,
    "output=s"	    =>\$output,
    "window-size=i"  =>\$windowsize,
    "step-size=i"   =>\$step,
    "test"          =>\$test,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
Test::runTests() if $test;

   
pod2usage(-msg=>"A input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"A output file has to be provided\n",-verbose=>1) unless $output;


      

open my $ofh, ">$output" or die "Could not open output file";


my $reader;

$reader=BpSlider->new($input,$windowsize,$step);


print "\n\nCalculating Dxy started at\t". localtime() ." ....\n";



while(my $window=$reader->nextWindow())
{
        my $chr=$window->{chr};
        my $pos=$window->{middle};
        my $win=$window->{window};
        my $above=$window->{count_covered};
        my $data=$window->{data};
	
        next unless @$data;
        
        my $coveredFrac=$above/$win;

        if ($coveredFrac>=$minCoverageFraction) {
	    
	    my($string)=Utility::calculateDxy($data,$pos,$chr);
	    #print "$coveredFrac\t$minCoverageFraction\n";
	    print $ofh "$string\n";
            #print "$string\n";
	    
	}
        else {
            print $ofh "$chr\t$pos\tna\tna\tna\n";
            #print "$chr\t$pos\tna\tna\tna\n";
        }
 
}


print "\n\nCalculating Dxy completed at\t". localtime() ." ....\n";


exit;



{
    use warnings;
    use strict;
    package BpSlider;

    sub new
    {
        my $class=shift;
        my $file=shift;
        my $window=shift;
        my $step=shift;
        
        open my $fh,"<$file" or die "Could not open file handle";
        
        return bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            file=>$file,
            fh=>$fh,
            curwin=>[],
            buffer=>[]
        },__PACKAGE__;
    }
    
    sub count_samples
    {
        my $self=shift;
        my $l=$self->_nextline();
        my $p=Utility::_parseLightwight($l);
        my $c=scalar(@{$p->{samples}});
        $self->_bufferline($l);
        return $c;
    }
    
    sub nextWindow
    {
        my $self=shift;

	
        #get the current window, and the current chromosome
        my $curwin=$self->{curwin};
        
        my $curChr="";
        $curChr=$curwin->[0]{chr} if @$curwin;
        
        my $resetchr=0;
        
        # empty unnecessary entries
        EMPTY: while(@$curwin)
        {
            my $e=shift @$curwin;
            if($e->{pos}>$self->{lower})
            {
                unshift @$curwin, $e;
                last EMPTY;
            }
            
        }
        
        # fill with novel entries
        my $line;
        FILL:while($line=$self->_nextline)
        {
            my $e=Utility::_parseLightwight($line);
            $curChr=$e->{chr} unless $curChr;
            
            
            if($e->{chr} eq $curChr && $e->{pos} <= $self->{upper})
            {
                push @$curwin,$e;
            }
            else
            {
                $resetchr=1 if $e->{chr} ne $curChr;
                $self->_bufferline($line);
                last FILL;
            }
        }
        
        return undef unless $curChr;
        
        
        my $toret=Utility::_annotateWindow($curwin,$curChr,$self->{lower},$self->{upper},$self->{window});
        
        if($resetchr or not defined($line))
        {
            # we transgressed the boundaries to the next chromosome
            # reset the windows and the current buffer
            $self->{lower}=0;
            $self->{upper}=$self->{window};
            $self->{curwin}=[];
        }
        else
        {
            # next time we will still be in the same chromosome
            # increase the upper and lower boundaries by the stepsize and set the current buffer
            $self->{upper}+=$self->{step};
            $self->{lower}+=$self->{step};
            $self->{curwin}=$curwin;
        }

        return $toret;
    }
    
    
    
    sub _nextline
    {
        my $self=shift;
        my $fh=$self->{fh};
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub _bufferline
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }
    
    
}



{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    
    
      
    sub calculateDxy {
	
        my $data=shift;
	my $pos=shift;
	my $chr=shift;
        
        my ($Dxy12,$Dxy13,$Dxy23) = (0,0,0);
        my ($ct12,$ct13,$ct23,$array_size) = (0,0,0,0);
        $array_size = @$data;
        
        foreach my $d (@$data)
        {
            #print "$d->{chr}\t$d->{pos}\t$d->{astate}\t$d->{dstate}\t$d->{sp1}\t$d->{sp2}\t$d->{sp3}\n";
            if ("$d->{sp1}" ne "$d->{sp2}") {
                $ct12++;
            }
            
            if ("$d->{sp1}" ne "$d->{sp3}") {
                $ct13++;
            }
            
            if ("$d->{sp2}" ne "$d->{sp3}") {
                $ct23++;
            }
            
        }
        
        # Calculate Dxy	    
        $Dxy12 = $ct12/$array_size;
        $Dxy13 = $ct13/$array_size;
        $Dxy23 = $ct23/$array_size;
        
        $Dxy12 = sprintf "%.10f",$Dxy12;
        $Dxy13 = sprintf "%.10f",$Dxy13;
        $Dxy23 = sprintf "%.10f",$Dxy23;

        #print $ofh "$chr\t$middle\t$Dxy12\t$Dxy13\t$Dxy23\n";
        
        my $string = "";
        $string = "$chr\t$pos\t$Dxy12\t$Dxy13\t$Dxy23";
        return $string;

    }
    
    sub _annotateWindow
    {
        my $curwin=shift;
        my $chr=shift;
        my $start=shift;
        my $end=shift;
        my $window=shift;

        #my $snps=0;
        my $aboveCoverage=0;
        foreach(@$curwin)
        {
            #$snps++ if $_->{ispuresnp};
            $aboveCoverage++ if $_->{iscov};
        }

        return
        {
            chr=>$chr,
            start=>$start,
            end=>$end,
            middle=>int(($end+1+$start)/2),
            count_covered=>$aboveCoverage,
            window=>$window,
            data=>$curwin      
        };
    }
    
    
    
    sub _parseLightwight
    {
        my $line=shift;
	
        chomp $line;
        my @a=split /\s+/,$line;
        my $chr=shift @a;
        my $pos=shift @a;
        my $astate=shift @a;
        my $dstate=shift @a;
	
	my $en={};
	if (($a[0] eq "-") or ($a[1] eq "-") or ($a[2] eq "-")) {

	    $en={
		chr=>$chr,
		pos=>$pos,
		astate=>$astate,
		dstate=>$dstate,
		iscov=>0,
		sp1=>$a[0],
		sp2=>$a[1],
		sp3=>$a[2]
	    };
        }

	else {
	    
	    $en={
		chr=>$chr,
		pos=>$pos,
		astate=>$astate,
		dstate=>$dstate,
		iscov=>1,
		sp1=>$a[0],
		sp2=>$a[1],
		sp3=>$a[2]
	    };
	    
	}
        
        #$en->{ispuresnp} = ($en->{issnp} and not $taintedsnp)?1:0;
        
        return $en;
    }
    
    
}








    
=head1 NAME

calculate-dxy.pl - Calculates the distance between 2 species within a given window in pairwise comparison.

=head1 SYNOPSIS

 perl calculate-dxy.pl --input mauve-parser-output.txt --output Dxy.txt --window-size 1000 --step-size 100\n";

=head1 OPTIONS

=over 4


=item B<--input>

The input file has to be mauve-parser.pl program output file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--window-size>

the size of the sliding window; default=1000

=item B<--step-size>

the size of the sliding window steps; default=100

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

A mauve-parser.pl program output file; example:

 2L	5783	G	C	C	G	G
 2L	5784	C	na	C	C	C
 2L	5785	na	na	C	T	G

 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: ancestral allelic state
 col 4: derived allelic state
 col 5: allelic state in species 1 (reference species)
 col 6: allelic state in species 2
 col 7: allelic state in species 3 (outgroup species)
 
=head2 Output

An output of this program looks like in the given example:

 2L	5788	0.50	0.50	0.30
 2L	5790	0.60	0.50	0.40
 2L	5791	0.44	0.10	0.33
 2L	5792	0.43	0.12	0.43

 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: D12 (Distance between species1 and species2) within a given window
 col 4: D13 (Distance between species1 and species3) within a given window
 col 5: D23 (Distance between species2 and species3) within a given window
 Note: species 1 refers to reference species; species 3 refers to outgroup species in phylogenetic tree; Example: species 1 is D. melanogaster, species2 is D. simulanes and species 3 is D. yakuba

=head1 AUTHORS

Ram vinay pandey

Robert Kofler

Pablo Orozco terWengel

Christian Schloetterer

=cut
