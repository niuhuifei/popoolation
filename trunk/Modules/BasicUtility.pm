{
    package BasicUtility;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT     =qw(reverse_complement_dna load_fasta);
    
    
    sub reverse_complement_dna
    {
        my $seq=shift;
        $seq=reverse($seq);
        $seq=~tr/ATCGNatcgn/TAGCNtagcn/;
        return $seq;
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