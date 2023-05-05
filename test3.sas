
   
 
**SAMPLE Data;
data sample(drop = i j ahbp1-ahbp8 abmi1-abmi8);

    call streaminit(5027);
    
    do i=1 to 1000;
        baseage = int( 35 + 25*rand('uniform'));
        
        array ahbp(8);
        array abmi(8);

        do j=1 to 8;
            ahbp(j) = (0.4>rand('uniform'));
            abmi(j) = log(25+5*(rand('normal')));
            end;
        
        do j=3 to 8 until (dia or censlost or censdead);
            id=i;
            time=j-3;
               
            age = baseage + 2*j;

            hbp     = ahbp(j);
            hbp_l1  = ahbp(j-1);
            hbp_l2  = ahbp(j-2);
            hbp_b   = ahbp(2);

            bmi     = abmi(j);
            bmi_l1  = abmi(j-1);
            bmi_l2  = abmi(j-2);
            bmi_b   = abmi(2);

            dia = ( (j/500) >rand('uniform'));
            censlost  = (0.05>rand('uniform'));
            censdead  = (0.05>rand('uniform'));
                       
            output;

            end;
       end;
   
run;

proc means data=sample;
title 'Means of SAMPLE data';
run;

**INTERV Calls;
%let interv1  =
    intno     = 1, 
    intlabel  = 'BMI Less Than 25 and No HBP',
    nintvar   = 2,
    intvar1   = bmi,
    inttype1  = 2,
    intmax1   = 3.22,
    inttimes1 = 1 2 3 4 5,
    intvar2   = hbp,
    inttype2  = 1,
    intvalue2 = 0,
    inttimes2 = 1 2 3 4 5;

    
**GFORMULA Call;
 
%include "gformula4.1p.sas";
*%include "dist.test3.sas" ;

options notes mprint nomlogic;
options autosignon sascmd="sas" ;
 

 
 

options notes mprintnest ;
options  nomlogic nosymbolgen nonotes nomprint ;

data _hosts ;
input host $ script $ ;
cards ;
host1 sas
host2 sas
;
run;

 data _hosts1 ;
 input host $21. script $13. ;
 cards;
jerzy.sph.harvard.edu tcpunix.scr
jerzy.sph.harvard.edu tcpunix.scr
jerzy.sph.harvard.edu tcpunix2.bat
jerzy.sph.harvard.edu tcpunix2.bat
;
run;

options comamid=tcp ;


%gformula(
    data= sample,
    id=id,
    time=time,
    timepoints = 5, 

    outc=dia,

    outctype = binsurv,
    
    fixedcov = hbp_b bmi_b,

    ncov=2,
    
    cov1 = hbp,    cov1otype = 2, cov1ptype = lag1bin,
    cov2 = bmi,    cov2otype = 3, cov2ptype = lag2cub, 
   
    seed= 9458,
   bootstrap_method = 0 ,
    nsimul=,
    nsamples = 500, 

numint=1 ,
runnc = 1,
refint = 0 ,

    simuldata=  roger1 ,
    survdata =  roger2 ,
    betadata=  roger3 ,
    covmeandata=  roger4   ,

   resultsdata = myresults ,

   hazardratio = 0,
   intcomp = 0 1 ,
   bootstrap_hazard = 0,

    check_cov_models = 1,
    

    run_remote = 1,
    check_remote = 0,
    use_sascmd=1 ,
    localhost = 1,
    nlocalhosts = 10 ,
    DistLog =   distributel.log ,  
    DistLst =   distributel.lst , 
    DIST_DETAIL = 1,      
    DIST_STATUS = 0,      
    DIST_DEBUG  = 0     
    );

 

    

