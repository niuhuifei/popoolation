{
    package Test::Variance;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use Test;
    use VarianceExactCorrection;
    use VarianceUncorrected;
    use FstMath;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT=qw(run_VarianceTests);
    our @EXPORT_OK = qw();
    
    
    sub run_VarianceTests
    {
        test_exact_Pi();
        test_exact_Theta();
        test_FstPi();
        test_exact_calculate_measure_Pi();        
    }
    
    sub troubleshoot
    {

    }
    
    sub test_exact_calculate_measure_Pi
    {
        my $vc;
        my $pi;
        
        $vc=VarianceExactCorrection->new(100,2);
        $pi=$vc->calculate_measure("pi",[{eucov=>4,A=>2,T=>2,C=>0,G=>0}],10);
        ok(abs($pi-0.20002) < 0.00001,"Testing small pi calculation: small pi is correct");
        
        $pi=$vc->calculate_measure("pi",[{eucov=>4,A=>2,T=>2,C=>0,G=>0},{eucov=>8,A=>2,T=>6,C=>0,G=>0}],10);
        ok(abs($pi-0.2600284) < 0.00001,"Testing small pi calculation: small pi is correct");
    }
    
    
    sub test_exact_Theta
    {
        my $vc;
        my $theta;
        
        # COVERAGE
        $vc=VarianceExactCorrection->new(100,2);
        $theta=$vc->_thetaSNP({eucov=>4,A=>2,T=>2,C=>0,G=>0});
        ok(abs($theta-2.0002) < 0.00001,"Testing exact theta correction: Theta value correct");
        
        $theta=$vc->_thetaSNP({eucov=>10,A=>4,T=>6,C=>0,G=>0});
        ok(abs($theta-0.582248) < 0.00001,"Testing exact theta correction: Theta value correct");        
        
        $theta=$vc->_thetaSNP({eucov=>20,A=>2,T=>18,C=>0,G=>0});
        ok(abs($theta-0.401039) < 0.00001,"Testing exact theta correction: Theta value correct");
        
        $theta=$vc->_thetaSNP({eucov=>100,A=>2,T=>98,C=>0,G=>0});
        ok(abs($theta-0.242309) < 0.00001,"Testing exact theta correction: Theta value correct");
        
        
        # MAF
        $vc=VarianceExactCorrection->new(100,1);
        $theta=$vc->_thetaSNP({eucov=>100,A=>4,T=>96,C=>0,G=>0});
        ok(abs($theta-0.211977) < 0.00001,"Testing exact theta correction: Theta value correct");        
        
        $vc=VarianceExactCorrection->new(100,3);
        $theta=$vc->_thetaSNP({eucov=>100,A=>0,T=>96,C=>4,G=>0});
        ok(abs($theta-0.273512) < 0.00001,"Testing exact theta correction: Theta value correct");

        $vc=VarianceExactCorrection->new(100,4);
        $theta=$vc->_thetaSNP({eucov=>100,A=>0,T=>96,C=>4,G=>0});
        ok(abs($theta-0.301776) < 0.00001,"Testing exact theta correction: Theta value correct");        
        
        # POOLSIZE
        $vc=VarianceExactCorrection->new(10,2);
        $theta=$vc->_thetaSNP({eucov=>100,A=>0,T=>60,C=>40,G=>0});
        ok(abs($theta-0.35353) < 0.00001,"Testing exact theta correction: Theta value correct");
        
        $vc=VarianceExactCorrection->new(50,2);
        $theta=$vc->_thetaSNP({eucov=>100,A=>0,T=>60,C=>40,G=>0});
        ok(abs($theta-0.248991) < 0.00001,"Testing exact theta correction: Theta value correct");        
        
        
        
        $vc=VarianceExactCorrection->new(100,2);
        $theta=$vc->_thetaSNP({eucov=>500,A=>2,T=>498,C=>0,G=>0});
        ok(abs($theta-0.194667) < 0.00001,"Testing exact theta correction using a high coverage of M=500: Theta value correct");        
        
        
        
    }
    
    
    sub test_exact_Pi
    {
        my $vc;
        my $pi;
        
        $vc=VarianceExactCorrection->new(50,2);
        $pi=$vc->_piSNP({eucov=>10,A=>5,T=>5,C=>0,G=>0});
        ok(abs($pi-0.714776)<0.00001,"Testing exact pi correction: Pi value correct");
        
        $pi=$vc->_piSNP({eucov=>100,A=>50,T=>50,C=>0,G=>0});
        ok(abs($pi-0.518725)<0.00001,"Testing exact pi correction: Pi value correct");        

        # pool 100; maf 1
        $vc=VarianceExactCorrection->new(100,1);
        $pi=$vc->_piSNP({eucov=>10,A=>5,T=>5,C=>0,G=>0});
        ok(abs($pi-0.561167)<0.00001,"Testing exact pi correction: Pi value correct");
        
        $pi=$vc->_piSNP({eucov=>10,A=>6,T=>4,C=>0,G=>0});
        ok(abs($pi-0.538721)<0.00001,"Testing exact pi correction: Pi value correct");
        
        $pi=$vc->_piSNP({eucov=>100,A=>70,T=>0,C=>30,G=>0});
        ok(abs($pi-0.428528)<0.00001,"Testing exact pi correction: Pi value correct");
        
        $pi=$vc->_piSNP({eucov=>110,A=>0,T=>20,C=>0,G=>90});
        ok(abs($pi-0.303283)<0.00001,"Testing exact pi correction: Pi value correct");
        
        # Coverage influence
        # pool 100; maf 2
        $vc=VarianceExactCorrection->new(100,2);
        $pi=$vc->_piSNP({eucov=>10,A=>0,T=>0,C=>6,G=>4});
        ok(abs($pi-0.685832) < 0.00001,"Testing exact pi correction: Pi value correct");
        
        $pi=$vc->_piSNP({eucov=>20,A=>0,T=>0,C=>12,G=>8});
        ok(abs($pi-0.564895) < 0.00001,"Testing exact pi correction: Pi value correct");

        $pi=$vc->_piSNP({eucov=>100,A=>0,T=>0,C=>60,G=>40});
        ok(abs($pi-0.495659) < 0.00001,"Testing exact pi correction: Pi value correct");

        # MAF Influence
        $vc=VarianceExactCorrection->new(100,3);
        $pi=$vc->_piSNP({eucov=>100,A=>0,T=>60,C=>0,G=>40});
        ok(abs($pi-0.505289) < 0.00001,"Testing exact pi correction: Pi value correct");

        $vc=VarianceExactCorrection->new(100,4);
        $pi=$vc->_piSNP({eucov=>100,A=>0,T=>0,C=>60,G=>40});
        ok(abs($pi-0.516117) < 0.00001,"Testing exact pi correction: Pi value correct");

        # POOLSIZE
        $vc=VarianceExactCorrection->new(10,2);
        $pi=$vc->_piSNP({eucov=>100,A=>0,T=>60,C=>0,G=>40});
        ok(abs($pi-0.538724) < 0.00001,"Testing exact pi correction: Pi value correct");

        $vc=VarianceExactCorrection->new(50,2);
        $pi=$vc->_piSNP({eucov=>100,A=>0,T=>60,C=>0,G=>40});
        ok(abs($pi-0.497976) < 0.00001,"Testing exact pi correction: Pi value correct");

        $vc=VarianceExactCorrection->new(100,1);
        $pi=$vc->_piSNP({eucov=>500,A=>0,T=>20,C=>0,G=>480});
        ok(abs($pi-0.0777312)<0.00001,"Testing exact pi correction using a coverage of 500: Pi value correct");
       
    }

    
    sub test_FstPi
    {

        my $pi;
    }
}


1;
