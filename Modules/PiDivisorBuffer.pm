{
package PiDivisorBuffer;
use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use Test;
use VarMath;


    sub new
    {
        my $class=shift;
        my $poolsize=shift;
        my $b=shift;
        my $maxbuffer=shift;
        
        my $toret={
            b=>$b,
            n=>$poolsize,
            maxbuf=>$maxbuffer,
            buffer=>[]
            };
        return bless $toret, __PACKAGE__;
    }
    
    sub get_divisor
    {
        my $self=shift;
        my $M=shift;
        return $self->{buffer}[$M] if $self->{buffer}[$M];
        my $div=$self->_calculate_divisor($M);
        
        # store in buffer if smaller than maxbuffer
        $self->{buffer}[$M]=$div if $M<=$self->{maxbuf};
        return $div;        
    }
    
    
    sub _calculate_divisor
    {
        my $self=shift;
        my $M=shift;
        my $b=$self->{b};
        my $n=$self->{n};
        
        my $divisor=0;
        for my $m ($b..$M-$b)
        {
            my $term1=(2*$m*($M-$m))/($M*($M-1));
            $term1*=a_Mnm($M,$n,$m);
            $divisor+=$term1;
        }
        return $divisor;
    }


}


1;
