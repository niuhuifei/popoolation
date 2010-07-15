package VarMath;
use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT=qw(noverk a_Mnm);
our @EXPORT_OK = qw();

sub noverk
{
    my $n=shift;
    my $k=shift;
    die "n over k; n has to be larger than zero" unless $n>0;
    die "n over k; k has to be larger than zero" unless $k>0;
    die "$k mus not be larger than $n" if $k>$n;
    
    my @above=(($n-$k+1)..$n);
    my @below=(1..$k);
    
    my $val=1;
    while(@above and @below)
    {
        if($val<1)
        {
            $val*=shift @above;
        }
        else
        {
            $val/=shift @below;
        }
    }
    
    foreach(@above)
    {
        $val*=$_;
    }
    foreach(@below)
    {
        $val/=$_;
    }
    return $val
}

sub binomial_term
{
    my $M=shift; # coverage
    my $n=shift; # pool size
    my $m=shift; # running variable for $b..$M-$b
    my $k=shift; # running variable for 1..$n-1
    
    my $val=noverk($M,$m);
    die "$val is zero for M: $M and m: $m\n" unless $val;
    my $t1=($k/$n)**$m;
    my $t2=(($n-$k)/$n)**($M-$m);
    my $toret=$val*$t1*$t2;
    return $toret;
}

sub a_Mnm
{
    my $M=shift;
    my $n=shift;
    my $m=shift;
    
    my $toret=0;
    foreach my $k (1..$n-1)
    {
        my $t1=binomial_term($M,$n,$m,$k);
        $t1*=(1/$k);
        $toret+=$t1;
    }
    return $toret;
}



1;
