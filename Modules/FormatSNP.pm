{
    package FormatSNP;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT =qw(get_WindowSNPFormater get_gtfSNPFormater get_SNPwriter get_syn_nonsyn_SNPFormater);
    
    
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
    
    
    sub get_syn_nonsyn_SNPFormater
    {
        my $outfile=shift;
        open my $ofh, ">",$outfile or die "Could not open output file $outfile";
        
        return sub
        {
            my $triplets=shift;
            my $chr=shift;
            my $start=shift;
            my $end=shift;
            
            my @snptriplets =grep {$_->{count_snps}==1} @$triplets;
            my $snpcount=@snptriplets;
            return unless @snptriplets;
            
            
            print $ofh ">$chr:$start-$end snps: $snpcount\n";
            
            foreach my $snp (@snptriplets)
            {
                my $af=$snp->{affected_snp};
                my $codon_old=$snp->{codon};
                my $codon_novel=$snp->{snp_codon};
                my $strand=$snp->{strand};
                my $aa_old=$snp->{aa_old};
                my $aa_new=$snp->{aa_new};
                my $syn=$snp->{syn};
                $syn=$syn?"syn":"non-syn";
                
                my $codon="$codon_old->$codon_novel";
                my $aa="$aa_old->$aa_new";
                
                my $toprint="$af->{chr}\t$af->{pos}\t$af->{refc}\t$af->{eucov}\t$af->{A}\t$af->{T}\t$af->{C}\t$af->{G}\t$af->{N}\t$syn\t$strand\t$codon\t$aa";
                print $ofh $toprint."\n";
                
            }
            print $ofh "\n";
            
            
            # data, chr, lower, upper, window, count_codons, cound_valid_codons, 
            
        }
    }
    
    


}

1;