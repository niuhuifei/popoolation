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
        
        if($measure eq "pi")
        {
            return $self->_calculate_pi($snps,$covercount);
        }
        elsif($measure eq "theta")
        {
            return $self->_calculate_theta($snps,$covercount);
        }
        elsif($measure eq "d")
        {
            return $self->_calculate_d($snps,$covercount);
        }
        else
        {
            die "unknown measure to calculate $measure";
        }
    }
    
    
    sub _calculate_pi
    {
        my $self=shift;
        my $snps=shift;
        my $covercount=shift;
    
        my $measurecalculater=$self->{pi};
        my $pi_sum=$measurecalculater->($self->{b},$self->{n},$snps);

        my $toret=0;
        $toret=$pi_sum/$covercount if $covercount;
        return $toret;
        
    }
    
    
    sub _calculate_theta
    {
        my $self=shift;
        my $snps=shift;
        my $covercount=shift;
        my $measurecalculator=$self->{theta};
        my $theta=$measurecalculator->($self->{b},$self->{n},$snps);
        my $toret=0;
        $toret=$theta/$covercount if $covercount;
        return $toret;
    }
    
    
    sub _calculate_d
    {
        my $self=shift;
        my $snps=shift;
        my $covercount=shift;
        my $measurecalculator=$self->{d};
        my $d=$measurecalculator->($self->{b},$self->{n},$snps);
        my $toret=0;
        $toret=$d if $covercount;
        return $toret;
    }
}
1;
