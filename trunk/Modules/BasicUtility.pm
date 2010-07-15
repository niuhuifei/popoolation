{
    package BasicUtility;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT     =qw(reverse_complement_dna);
    
    
    sub reverse_complement_dna
    {
        my $seq=shift;
        $seq=reverse($seq);
        $seq=~tr/ATCGNatcgn/TAGCNtagcn/;
        return $seq;
    }
}