{
    package VarianceUncorrected;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use ThetaDivisorBuffer;
    use PiDivisorBuffer;
    
    my $asum;
    
    sub new {
        my $class = shift;
        my $poolSize=shift;
        my $maf=shift;
        my $maxbuffer=shift || 100_000;
        
        $asum=_calculate_a_sum($maxbuffer) unless $asum;
        
        
        my $self = bless {
                          n=>$poolSize,
                          b=>$maf,
                          maxbuf=>$maxbuffer,

                          }, __PACKAGE__;
        return $self;
    }
    
    
    sub _calculate_a_sum
    {
        my $maxbuffer=shift;
        my $tr=[];
        $tr->[0]=0;
        $tr->[1]=1;
        
        for(my $i=2; $i<$maxbuffer; $i++)
        {
            $tr->[$i]= $tr->[$i-1]+(1/$i);
        }
        return $tr;
    }
    
    sub calculate_measure
    {
        my $self=shift;
        my $measure=shift;
        my $snps=shift;
        my $covercount=shift;
        
        die "method not implemented";
        
        
    }


}


1;
