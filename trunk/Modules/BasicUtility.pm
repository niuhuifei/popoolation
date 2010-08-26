{
    package BasicUtility;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT     =qw(reverse_complement_dna load_fasta get_fasta_reader get_fasta_writer);
    
    
    sub reverse_complement_dna
    {
        my $seq=shift;
        $seq=reverse($seq);
        $seq=~tr/ATCGNatcgn/TAGCNtagcn/;
        return $seq;
    }
    
    sub get_fasta_reader
    {
        my $file=shift;
        open my $ifh,"<",$file or die "Could not open input file";
        my $lastheader=\"";
        
        return sub
        {
            my $header="";
            my $sequence="";
            while(my $l=<$ifh>)
            {
                chomp $l;
                if($l=~m/^>/)
                {
                    my($newheader)=$l=~m/^>(.+)/;
                    if($$lastheader)
                    {
                        $header=$$lastheader;
                        $lastheader=\$newheader;
                        return {
                            head=>$header,
                            seq=>$sequence
                        };
                    }
                    else
                    {
                        $lastheader=\$newheader;
                    }
                    
                }
                else
                {
                    $sequence.=$l;
                }
            }
            if($sequence)
            {
                return
                {
                  head=>$$lastheader,
                  seq=>$sequence
                };
            }
            else
            {
                return undef;
            }
            
        }
        
    }
    
    sub get_fasta_writer
    {
        my $file=shift;
        my $leng=shift || 50;
        open my $ofh, ">",$file or die "Could not open output file";

        
        return sub
        {
            my $fasta=shift;
            my $header=$fasta->{head};
            my $seq=$fasta->{seq};
            print $ofh ">$header\n";
            
            for(my $i=0; $i<length($seq); $i+=$leng)
            {
                my $subseq=substr($seq,$i,$leng);
                print $ofh $subseq."\n";
            }
        }
    }
    
    sub load_fasta
    {
        my $reffile=shift;
        open my $ifh, "<", $reffile or die "Could not open reference file";
        my $refhash={};
        
        my $header="";
        my $seq="";
        while(my $l=<$ifh>)
        {
            chomp $l;
            if($l=~m/^>/)
            {
                if($header)
                {
                    die "Reference $header already exists" if exists($refhash->{$header});
                    $refhash->{$header}=$seq;
                }
                ($header)=$l=~m/^>(\S+)/;
                $seq="";
            }
            else
            {
                $seq.=$l;
            }
        }
        
        if($header)
        {
            die "Reference $header already exists" if exists($refhash->{$header});
            $refhash->{$header}=$seq;
        }
        
        return $refhash;
    }
}
1;