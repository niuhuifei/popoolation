{
    package VarianceExactCorrection;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use VarMath;

    sub new {
        my $class = shift;
        my $poolSize=shift;
        my $maf=shift;

        
        
        # get_theta_calculator($b,$n,$snp)
        # get_pi_calculator($b,$n,$snp)

        
        my $self = bless {
                          n=>$poolSize,
                          b=>$maf,
                          pi=>get_pi_calculator(),
                          theta=>get_theta_calculator(),
                          d=>get_D_calculator()
                          }, __PACKAGE__;
        return $self;
    }
    

    
    
    sub calculate_measure
    {
        my $self=shift;
        my $measure=shift;
        my $snps=shift;
        my $covercount=shift;
        
        $measure=lc($measure);
        
        my $measure_sum=0;
        foreach my $snp(@$snps)
        {
            
            my $meas_snp=$self->{$measure}->($self->{b}, $self->{n}, $snp);
            $measure_sum+=$meas_snp;
        }
        
        my $toret=0;
        $toret=$measure_sum/$covercount if $covercount;
        return $toret;
    }

    
    
    

    


}
1;
