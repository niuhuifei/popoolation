{
    package VarianceExactCorrection;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use ThetaDivisorBuffer;
    use PiDivisorBuffer;

    sub new {
        my $class = shift;
        my $poolSize=shift;
        my $maf=shift;
        my $maxbuffer=shift || 1000_000;
        
        
        my $self = bless {
                          n=>$poolSize,
                          b=>$maf,
                          pidivbuf=>PiDivisorBuffer->new($poolSize,$maf,$maxbuffer),
                          thetadivbuf=>ThetaDivisorBuffer->new($poolSize,$maf,$maxbuffer),
                          maxbuf=>$maxbuffer,

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
            return $self->calculate_Pi($snps,$covercount);
        }
        elsif($measure eq "theta")
        {
            return $self->calculate_Theta($snps,$covercount);
        }
        elsif($measure eq "d")
        {
            return $self->calculate_d($snps,$covercount);
        }
        else
        {
            die "measure not supported"
        }
        
    }

    
    
    sub calculate_Theta
    {
        my $self=shift;
        my $snps=shift; # a simple array of snps
        my $covercount=shift;
    
        my $theta_sum=0;
        foreach my $snp(@$snps)
        {
            
            my $theta_snp=$self->_thetaSNP($snp);
            $theta_sum+=$theta_snp;
        }
        
        my $toret=0;
        $toret=$theta_sum/$covercount if $covercount;
        return $toret;
    }

    sub calculate_d
    {
        my $self=shift;
        my $snps=shift; # a simple array of snps
        my $covercount=shift;
        
        return $self->calculate_Pi($snps,$covercount) - $self->calculate_Theta($snps,$covercount);
    }


    sub calculate_Pi
    {
        my $self=shift;
        my $snps=shift; # a simple array of snps
        my $covercount=shift;
    
        my $pi_sum=0;
        foreach my $snp (@$snps)
        {
            my $pi_snp = $self->_piSNP($snp);
            $pi_sum += $pi_snp;
        }
        my $toret=0;
        $toret=$pi_sum/$covercount if $covercount;
        return $toret;
    }
    
    sub _thetaSNP
    {
        my $self=shift;
        my $snp=shift;
        my $thetadivbuffer=$self->{thetadivbuf};
        my $theta_snp=1 / ($thetadivbuffer->get_divisor($snp->{eucov}));
       return $theta_snp;
    }
    
    sub _piSNP
    {
        my $self=shift;
        my $snp = shift;
        my $pidivbuffer=$self->{pidivbuf};
        
        my $M=$snp->{eucov};
        my $pi_snp=1;
        $pi_snp-=($snp->{A}/$M)**2;
        $pi_snp-=($snp->{T}/$M)**2;
        $pi_snp-=($snp->{C}/$M)**2;
        $pi_snp-=($snp->{G}/$M)**2;
        $pi_snp*=$M/($M-1);
        
        $pi_snp/=$pidivbuffer->get_divisor($M);
        return $pi_snp;
    }

}
1;
