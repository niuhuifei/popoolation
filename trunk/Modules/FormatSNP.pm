{
    package FormatSNP;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT =qw(get_WindowSNPFormater get_gtfSNPFormater get_SNPwriter);
    
    
            #pos=>$pos,
            #chr=>$chr,
            #refc=>$rc,
            #consc=>"N",
            #consc_confidence=>0,
            #totcov=>0,
            #eucov=>0,
            #A=>0,
            #T=>0,
            #C=>0,
            #G=>0,
            #del=>0,
            #N=>0,
            #iscov=>0,
            #issnp=>0,
            #ispuresnp=>0
            
    sub get_SNPwriter
    {
        my $file =shift;
        open my $ofh, ">",$file or die "Could not open SNP output file";
        return sub
        {
            my $snp=shift;
            print $ofh _formatSNP($snp);
        }
    }
            
    sub get_gtfSNPFormater
    {
        my $file =shift;
        open my $ofh, ">",$file or die "Could not open SNP output file";
        return sub
        {
            my $geneid=shift;
            my $snps=shift;
            my $snpcount=@$snps;
            print $ofh ">$geneid snps:$snpcount\n";
            foreach my $snp (@$snps)
            {
                print $ofh _formatSNP($snp);
            }
            print $ofh "\n";
            
            
        }
        
    }
    
    sub get_WindowSNPFormater
    {
        
        my $file=shift;
        open my $ofh, ">", $file or die "Could not open SNP output file";
        
        return sub
        {
            my $win=shift;
            my $chr=$win->{chr};
            my $start=$win->{start};
            my $end=$win->{end};
            my $middle=$win->{middle};
            my $window=$win->{window};
            my $data=$win->{data};
            my $snpcount=$win->{count_snp};
                
            my $snps=[];
            foreach(@$data)
            {
                push @$snps,$_ if $_->{ispuresnp};
            }
                

            print $ofh ">$chr:$middle $chr:$start-$end snps:$snpcount\n";
            foreach my $snp (@$snps)
            {
                print $ofh _formatSNP($snp);
            }
            print $ofh "\n";
        }
    }
    
    
    sub _formatSNP
    {
        my $snp=shift;
        return "$snp->{chr}\t$snp->{pos}\t$snp->{refc}\t$snp->{eucov}\t$snp->{A}\t$snp->{T}\t$snp->{C}\t$snp->{G}\t$snp->{N}\n";
    }
    
    


}

1;