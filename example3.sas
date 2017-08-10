
options linesize=88 pagesize=54;
options nonotes;
%include "gformula3.sas";
 
%macro create_sample(event = dia) ;
%let condition = ;
%if %upcase(&event) = DIA %then %let condition = dia or censlost or dead ;
%else %if %upcase(&event) = CONT_E %then %let condition = censlost   ;
%else %if %upcase(&event) = BIN_E %then %let condition = dead or  censlost   ;

 
**SAMPLE Data;
    data sample(drop = i j ahbp1-ahbp12 abmi1-abmi12 /* hbp_b bmi_b */);

        call streaminit(5027);
        
        do i=1 to 1000;
            baseage = int( 35 + 25*rand('uniform'));
            
            array ahbp(12);
            array abmi(12);

            do j=1 to 12;
                ahbp(j) = (0.4>rand('uniform'));
                if j > 1 & ahbp(j-1) = 1 then ahbp(j) = 1 ;
                abmi(j) = round((25+5*(rand('normal'))),0.001);
            end;
            
            do j=3 to 12  until ( &condition   ) ;
                id=i;
                time=j-3;
                                    
                hbp     = ahbp(j);
                hbp_l1  = ahbp(j-1);
                hbp_l2  = ahbp(j-2);
                hbp_b   = ahbp(3);

                bmi     = abmi(j);
                bmi_l1  = abmi(j-1);
                bmi_l2  = abmi(j-2);
                bmi_b   = abmi(3);

                dia = ( (j/500) >rand('uniform'));
                

                if time < 9 then censlost  = (0.05>rand('uniform'));
                else censlost = 0 ;
                %if %upcase(&event) = DIA  or %upcase(&event) = BIN_E %then dead  = (0.05>rand('uniform'));;
                if time = 9 then dead = . ;
                if time = 9 then cont_e =  round((bmi+5*(rand('normal'))),0.01)  ;
                else bmi_e = . ; 
                if time = 9 then bin_e = rand('bernoulli',0.6);
                else bin_e = . ;
                       
                output;

            end;
    end;
       
    run;

    data sample ;
    set sample ;
    %if %upcase(&event)=DIA %then %do;
      if censlost = 1 then do ;
          dead = . ;                           
      end;   
      if censlost = 1 or dead = 1 then do ;
           
           dia = . ;
      end;
    %end;
    %if %upcase(&event) = CONT_E %then %do;
       if time < 9 then cont_e = . ;        
       if censlost = 1  then do ;
          
           cont_e = . ;
      end;
    %end;
    %if %upcase(&event)=BIN_E %then %do;
        if censlost = 1 then do ;
            dead = . ;             
        end ;       
        if time < 9 then bin_e = . ;
        if censlost = 1 or dead = 1 then do ;
            
           bin_e = . ;
       end;
    %end;

    run;


    proc means data=sample;
    title 'Means of SAMPLE data';
    run;
%mend ;

 

ods graphics off ;



**INTERV Calls;
%let interv1  =
    intno     = 1, 
    intlabel  = 'BMI Less Than 25 and No HBP',
    nintvar   = 2,
    intvar1   = bmi,
    inttype1  = 2,
    intmax1   = 25,
    inttimes1 = 0 1 2 3 4 5 6 7 8 9,
    intvar2   = hbp,
    inttype2  = 1,
    intvalue2 = 0,
    inttimes2 = 0 1 2 3 4 5 6 7 8 9 ;

%let interv2  =
    intno     = 2, 
    intlabel  = '50% Chance of 10% BMI Reduction on HBP Dx',
    intcond   = ( hbp= 1 and hbp_l1 = 0),
    nintvar   = 1,
    intvar1   = bmi,
    inttype1  = 3,
    intchg1   = -0.1,
    inttimes1 = 0 1 2 3 4 5 6 7 8 9,
    intpr1    = 0.5;
    
**GFORMULA Call;
title 'GFORMULA SAMPLE';

options nomprint nonotes ;
 
%create_sample(event = bin_e ) ;

 

proc datasets library = work nolist ;
save  sample ;
run;
quit;

%gformula(
    data= sample,
    id=id,
    time=time,
    timepoints = 10, 

    outc= bin_e ,
    outctype= bineofu   ,
    outcinteract = 0*1  ,
    comprisk =   ,          
    fixedcov =  hbp bmi  baseage  , /* using hbp and bmi for fixedcov forces the corresponding baseline variables into each model */
    ncov=2,
    timeptype=concat, 
    timeknots= 1 2 3 4 5 6 7 8 9, 
    cov1 = hbp,    cov1otype = 2, cov1ptype = lag1bin    ,   
    cov2 = bmi,    cov2otype = 3, cov2ptype = lag2cub  ,  
    seed= 9458,
    check_cov_models = 1 ,
    print_cov_means = 0,
    save_raw_covmean = 1,
    /* datasets */
    savelib = work,
    simuldata = simul0 ,
    survdata =    mysurv0,
    covmeandata =  mycovmean0 , 
    intervname =  myinterv ,
    observed_surv= myobssurv0,
    betadata = betadata0 ,
    nsimul= 1000 ,
    nsamples = 20, 
    sample_start = 0 ,
    sample_end = -1 ,
    resultsdata = myresults0,
    numint=2 ,
    rungraphs = 0,
    graphfile=bin_e.pdf ,
     printlogstats = 0 
    );

 

 
 proc datasets library= work ;
 save sample mysurv0  mycovmean0  myobssurv0  mycovmean0_raw   myresults0 betadata0;
 quit;

* example run with continuous outcome measured only at end of follow-up using a truncated normal model. Variable of interest
  is the difference of bmi between end and start of follow-up. Here there is no competing risk, only censoring due to lost of follow-up. ;


%gformula(
    data= sample,
    id=id,
    time=time,
    timepoints = 10, 
    outc=bin_e ,
    outctype= bineofu   ,
    outcinteract = 0*1  ,
    comprisk =  ,        
    fixedcov =  hbp bmi  baseage ,
    ncov=2,
    timeptype=concat, timeknots= 1 2 3 4 5 6 7 8 9 ,
    cov1 = hbp,    cov1otype = 2, cov1ptype = lag1bin,
    cov2 = bmi,    cov2otype = 3, cov2ptype = lag2cub  , 
    seed= 9458,
    check_cov_models = 1 ,
    print_cov_means = 0,
    save_raw_covmean = 1,
    /* datasets */
    savelib = work ,
    survdata =    mysurv,
    covmeandata =  mycovmean , 
    intervname =   myinterv ,
    observed_surv= myobssurv,
    betadata = betadata0a ,
    nsimul= 1000 ,
    nsamples = 20, 
    sample_start = 0 ,
    sample_end = 10 ,
    numint=2 ,
    rungraphs = 0 ,
    printlogstats = 0 
    );

  
 proc datasets library= work ;
 save sample mysurv0  mycovmean0  myobssurv0  mycovmean0_raw  myresults0 betadata0 betadata0a   mysurv_0_10   mycovmean_0_10   myinterv_0_10 ;
 quit;



    %gformula(
    data= sample,
    id=id,
    time=time,
    timepoints = 10, 
    outc=bin_e ,
    outctype= bineofu   ,
    outcinteract = 0*1  ,
    comprisk = ,     
    fixedcov =   hbp bmi  baseage ,
    ncov=2,
    timeptype=concat, timeknots= 1 2 3 4 5 6 7 8 9 ,
    cov1 = hbp,    cov1otype = 2, cov1ptype = lag1bin,
    cov2 = bmi,    cov2otype = 3, cov2ptype = lag2cub  , 
    seed= 9458,
    check_cov_models = 1 ,
    print_cov_means = 0,
    save_raw_covmean = 1,

    /* datasets */
    savelib = work ,
    survdata =    mysurv,
    covmeandata =  mycovmean , 
     intervname =   myinterv ,
    observed_surv= myobssurv,
    betadata = betadata0b ,
    nsimul= 1000 ,
    nsamples = 20, 
    sample_start = 11 ,
    sample_end = 20,

    numint=2 ,
    rungraphs = 0 ,
    printlogstats = 0 
    );


     
 proc datasets library= work ;
 save sample mysurv0  mycovmean0  myobssurv0  mycovmean0_raw  myresults0  betadata0 betadata0a betadata0b   mysurv_0_10   mycovmean_0_10   myinterv_0_10  
  mysurv_11_20   mycovmean_11_20    myinterv_11_20 ;
 quit;


 

      %bootstrap_results(
              bootlib =  work ,
              outc = bin_e,
              comprisk = ,
              outctype = bineofu ,
              bootname = myinterv ,
              check_cov_models = 1 ,
              covmeandata = mycovmean   , /* needed for graphs */
              observed_surv =  myobssurv  , /* needed for graphs */
              combine_survdata = 1 , /* for call to construct graphs */
              survdata=mysurv , /* needed for graphs */
              print_cov_means = 0,
              savecovmean = 0,
              time = time ,
              timepoints = 10,
              ncov = 2,
              numparts = 2,
              samplestart = 0 11 ,
              sampleend =   10 20 ,
              numboot = 20,
              numint = 2 ,
              refint = 0 ,
              resultsdata = myresults1 ,
              rungraphs = 1,              
              graphfile=bin_e2.pdf

              );

 

proc compare base = myresults0 compare= myresults1 ;
run;


data mybeta ;
set betadata0a betadata0b ;
run;

proc compare base = betadata0 compare = mybeta ;
run;

 
data mycovmean1_raw ;
set mycovmean_0_10 mycovmean_11_20 ;
run;

proc compare base= mycovmean0_raw compare =mycovmean1_raw ;
run;
 
proc compare base = mysurv0 compare= mysurv ;
run;

 
proc compare base = myobssurv0 compare= myobssurv ;
run;

proc compare base = mycovmean0 compare= mycovmean ;
run;
  
 
