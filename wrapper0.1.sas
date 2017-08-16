%include '/proj/sas_macros/gformula/May-2017/gformula3.sas ' ;



%macro create_sample ;
%let condition = ;
%let condition = dia or censlost or dead ;
 




**SAMPLE Data;
    data sample(drop = i j ahbp1-ahbp8 aact1-aact8 acont1 - acont8);

        call streaminit(5027);

        do i=1 to 1000;
            baseage = int( 35 + 25*rand('uniform'));

            array ahbp(8);
            array aact(8);
            array acont(8);

            do j=1 to 8;
                ahbp(j) = (0.2>rand('uniform'));
                acont(j) = rand('normal');
                if j > 1 & ahbp(j-1) = 1 then ahbp(j) = 1 ;

        aact(j)=(0.7>rand('uniform'));
        if aact(j)=1 then do;
                   aact(j) = int(exp(3.5+0.4*(rand('normal'))));
        end;
                end;

            do j=3 to 8  until ( &condition   ) ;
                id=i;
                time=j-3;
                

                hbp     = ahbp(j);
                hbp_l1  = ahbp(j-1);
                hbp_l2  = ahbp(j-2);
                

                act     = aact(j);
                act_l1  = aact(j-1);
                act_l2  = aact(j-2);


                cont = acont(j);
                cont_l1 = acont(j-1);
                cont_l2 = acont(j-2);
                

              dia = ( (j/500) >rand('uniform'));
                censlost  = (0.05>rand('uniform'));
                dead      = (0.05>rand('uniform'));
               



                output;

                end;
           end;




    run;


data sample;
set sample; 
call streaminit(1234);
if censlost=1 then do;
       dia= .;
       dead= .;
end;
else do;
   if dead=1 then dia= .;
end;
arm = rand('bernouli',0.6);
run;

proc means data=sample;
title 'Means of SAMPLE data';
run;
%mend ;

/***********/


%macro gformula_wrapper(
    data=,           
    id=,      
    time=,       
    timeptype = ,
    timeknots = ,
    timeinc = 1,
    timefuncgen = ,   
    timepoints=,    
    interval=1,      

    outc=,          
    outctype=binsurv,  
    outcinteract=,   
    outcwherem = (1=1) ,  
    outcwherenosim=(1=0),  
    outcnosimelsemacro =,  
    comprisk=,        
    compriskinteract=,   
    compriskwherem = (1=1) , 
    compriskwherenosim=(1=0),   
    comprisknosimelsemacro =, 
  
    arm = , /* variable for separating the data into two arms, need to have values {0,1} for  {a,b},   */
    numinta= 0 ,
    refinta = 0 ,
    runnca = 1,
    fixedcova=,      
    ncova=,          
    usebetadataa = 0 ,
    betadataa = ,
    cova1=,cova1otype=1,cova1ptype=,cova1mtype= all ,cova1cumint= ,cova1skip=-1,cova1inc=0,cova1knots=,cova1interact=,cova1wherem=(1=1),cova1wherenosim=(1=0),
    cova1nosimelsemacro=, cova1class=, cova1classelse=, cova1addvars=,cova1genmacro=,cova1modusermacro =,cova1moddatausermacro =,cova1setblvar =,
    cova1simusermacro = ,cova1barray = ,cova1sarray =,cova1randomvisitp=,cova1visitpmaxgap=9e10,cova1visitpwherem=(1=1),cova1visitpcount = ,

/** for second call **/
 
   numintb= 0 ,
   refintb = 0 ,
   runncb = 1,
   fixedcovb=, 
   ncovb= ,
   usebetadatab = 0 ,
   betadatab = ,
   covb1=,covb1otype=1,covb1ptype=,covb1mtype= all ,covb1cumint= ,covb1skip=-1,covb1inc=0,covb1knots=,covb1interact=,covb1wherem=(1=1),covb1wherenosim=(1=0),
   covb1nosimelsemacro=, covb1class=, covb1classelse=, covb1addvars=,covb1genmacro=,covb1modusermacro =,covb1moddatausermacro =,covb1setblvar =,
   covb1simusermacro = ,covb1barray = ,covb1sarray =,covb1randomvisitp=,covb1visitpmaxgap=9e10,covb1visitpwherem=(1=1),covb1visitpcount = ,

    /*other override options used for each gformula call */

    wherevars=,        /* list of variables referenced in any of the cov#wherem conditions, change JGY*/ 
    keepsimuldata=,    /*list of variables not created by a given ptype or otype that will be needed in simulated data set new change JGY*/
    equalitiessimuldata=,/*user defined macro that equates pre-baseline simulated vars to observed new change JGY*/
    eventaddvars=, /*list of variables to be added to event predictor list new change JGY*/
     
    compriskaddvars=,
    
    censladdvars=,
    

 
    
    simuldata=,     /* data set to store simulated data */
    resultsdata=,   /* data set to store results */
    survdata=,       /*data set to store cumulative survival probabilities at each time point under all interventions*/
    outputs = yes,  /* whether to print regression results */
    print_stats = 1 ,
    check_cov_models = 0, /* create data set for difference of mean of observed covs and mean of simulated covs under natural 
                              course */
    print_cov_means = 0,  /* print out tables of comparison of observed and simulated variables */
    covmeandata =  ,
    save_raw_covmean = 0,
    observed_surv = ,
    intervname = ,
    
   
    seed= 7834,          /* random numbers seed */
    nsamples=50,    /* number of bootstrap samples (default is 50) */
    nsimul=,        /* size (# subjects) of simulated sample, default is sample size */
    nparam=,        /* size (# subjects) of parameter sample, default is sample size */
    hazardratio=0 , /* calculate the hazard ratio for two interventions. This will increase the run time for the macro*/
    intcomp =    ,  /* needs to be a list with two numbers for comparing the hazard ratio when hazardratio = 1 */
    bootstrap_hazard = 0, /* when running bootstrap samples also include calculation of hazard ratios */
    hazardname =  , /* name of data set to hold hazard ratio when runnig bootstraps in parts, data will be saved in savelib. */
   
    sample_start = 0 ,  /* first sample to use in bootstraps (can be 0 for original data ) */
    sample_end = -1,     /* last sample to use in bootstraps (should be at most equal to nsamples) */
    savelib = work,  /* location for saving intermediate results for chunks */


    rungraphs = 0 ,
    title1a=,
    title2a=,
    title3a=,
    titledataa= ,
    graphfilea=gfilea.pdf ,
    tsizea=1 ,
    title1b=,
    title2b=,
    title3b=,
    titledatab= ,
    graphfileb=gfileb.pdf ,
    tsizeb=1 ,
    runnc = 1 ,
    weight = ,
    printlogstats = 1,
    usespline = 1 ,
    checkaddvars = 1,
    minimalistic = no, /* only keep results for outcome variables */
    testing = no /* keep each simulated data set for each intervention. will be named savelib.simulated&intno */
    )  /parmbuff;

 

    %let suffixlist=  otype  ptype  mtype   cumint  skip inc knots interact wherem wherenosim
     nosimelsemacro class classelse addvars genmacro modusermacro  moddatausermacro  setblvar simusermacro barray  
     sarray  randomvisitp visitpmaxgap visitpwherem  visitpcount ;


    %do i = 1 %to &ncova ;
        %let covlist = cova&i ;        
            %do j = 1 %to 25 ; 
                %let word = %scan(&suffixlist,&j) ;                
                %let covlist = &covlist cova&i.&word  ;
               
            %end;             
            %local &covlist ;
            %let cova&i.otype=1  ;
            %let cova&i.ptype=  ;
            %let cova&i.cumint=   ;    
            %let cova&i.mtype= all  ;
            %let cova&i.knots= ;     
            %let cova&i.skip=-1 ;     
            %let cova&i.inc=0 ;    
            %let cova&i.interact= ;  
            %let cova&i.wherem= (1=1) ;                   
            %let cova&i.wherenosim=(1=0) ; 
            %let cova&i.nosimelsemacro = ; 
            %let cova&i.class= ;    
            %let cova&i.classelse= ;     
            %let cova&i.addvars = ;
            %let cova&i.genmacro = ; 
            %let cova&i.modusermacro = ; 
            %let cova&i.moddatausermacro = ;  
            %let cova&i.setblvar = ;  
            %let cova&i.simusermacro = ;  
            %let cova&i.barray = ;   
            %let cova&i.sarray = ; 
            %let cova&i.randomvisitp = ;  
            %let cova&i.visitpmaxgap=9e10 ;     
            %let cova&i.visitpwherem=(1=1) ; 
            %let cova&i.visitpcount=  ;  
      %end;
 
 
     
    %do i = 1 %to &ncovb ;
        %let covlist = covb&i ;        
            %do j = 1 %to 25 ; 
                %let word = %scan(&suffixlist,&j) ;                
                %let covlist = &covlist covb&i.&word  ;               
            %end;             
            %local &covlist ;
 
            %let covb&i.otype=1  ;
            %let covb&i.ptype=  ;
            %let covb&i.cumint=   ;    
            %let covb&i.mtype= all  ;
            %let covb&i.knots= ;     
            %let covb&i.skip=-1 ;     
            %let covb&i.inc=0 ;    
            %let covb&i.interact= ;  
            %let covb&i.wherem= (1=1) ;                   
            %let covb&i.wherenosim=(1=0) ; 
            %let covb&i.nosimelsemacro = ; 
            %let covb&i.class= ;    
            %let covb&i.classelse= ;     
            %let covb&i.addvars = ;
            %let covb&i.genmacro = ; 
            %let covb&i.modusermacro = ; 
            %let covb&i.moddatausermacro = ;  
            %let covb&i.setblvar = ;  
            %let covb&i.simusermacro = ;  
            %let covb&i.barray = ;   
            %let covb&i.sarray = ; 
            %let covb&i.randomvisitp = ;  
            %let covb&i.visitpmaxgap=9e10 ;     
            %let covb&i.visitpwherem=(1=1) ; 
            %let covb&i.visitpcount=  ;  
    %end;
    
     


   %let mycount = %sysfunc(countw(&syspbuff,','));
  
 %do i = 1 %to &mycount ;  
    %let mydef = %scan(&syspbuff,&i,%str((),));   
    %let myvar = %scan(&mydef,1,%str(=));
    %let myval = %scan(&mydef,2,%str(=));
    %let &myvar = &myval ;
 %end;

 %if &ncova > 0 and &ncovb > 0 %then %do;
    %if %bquote(&survdata) = %then %let survdata = survholder ;
 %end;
    
 %if &ncova > 0 %then %modelA ;
 %if &ncovb > 0 %then %modelB;


 %let genresults = 0;
 %if &numinta = 0 and &numintb = 0 %then %let genresults = 1;
 %if &numinta = 0 and &numintb > 0 and &runncb = 1 %then %let genresults = 1 ;
 %if &numinta > 0 and &runnca = 1 and &numintb = 0  %then %let genresults = 1;
 %if &numinta > 0 and &numintb > 0 %then %let genresults = 1;

 %if &genresults = 1 %then %results2;
 %if &genresults = 0 %then %do ;
    %put NO COMPARISON DONE DUE TO NO MATCHING INTERVENTIONS ;
    %put( numinta = &numinta , runnca = &runnca , startint = %eval(1-&runnca) ) ( numintb = &numintb , runncb = &runncb , startint = %eval(1 - &runncb)) ;
%end; 
     
%mend ;
%macro modelA ;

 data dataa ;
 set &data(where = (&arm = 0));
 run;


     
    %let simuldataa = ;
    %let resultsdataa= ;
    %let survdataa = ;
    %let covmeandataa = ;
    %let observed_surva = ;
    %let intervnamea = ;
    

    %if %bquote(&simuldata)^= %then %let simuldataa = &simuldata.a ;          
    %if %bquote(&resultsdata)^=  %then %let resultsdataa = &resultsdata.a;
    %if %bquote(&survdata)^=  %then %let survdataa = &survdata.a  ;
    %if %bquote(&covmeandata)^=  %then %let covmeandataa = &covmeandata.a ;      
    %if %bquote(&observed_surv)^=  %then %let observed_surva = &observed_surv.a; 
    %if %bquote(&intervname)^=  %then %let intervnamea = &intervname.a ;
    
    
    %if &numinta > 0 %then %do;
        %do i = 1 %to &numinta ;
            %let interv&i = &&interva&i ;
            %put interv&i for modelA ::: &&interv&i ;
        %end;
    %end;
 
 %gformula(
    data= dataa,           
    id= &id,      
    time= &time,       
    timeptype = &timeptype ,
    timeknots = &timeknots,
    timeinc = &timeinc,
    timefuncgen = &timefuncgen,   
    timepoints= &timepoints,    
    interval=&interval,      

    outc= &outc,          
    outctype=&outctype,  
    outcinteract=&outcinteract,   
    outcwherem = &outcwherem ,  
    outcwherenosim= &outcwherenosim,  
    outcnosimelsemacro = &outcnosimelsemacro,  
    comprisk= &comprisk,        
    compriskinteract= &compriskinteract,   
    compriskwherem = &compriskwherem , 
    compriskwherenosim= &compriskwherenosim,   
    comprisknosimelsemacro = &comprisknosimelsemacro, 
  
     
    numint= &numinta ,
    fixedcov= &fixedcova,      
    ncov= &ncova,          

    
    %do i = 1 %to &ncova ;
        cov&i = &&cova&i ,
        %do j = 1 %to 25 ;
            %let word = %scan(&suffixlist,&j);
            cov&i.&word = &&cova&i.&word ,
        %end;
    %end;
    /*other override options used for each gformula call */

    wherevars= &wherevars,        
    keepsimuldata= &keepsimuldata,     
    equalitiessimuldata= &equalitiessimuldata, 
    eventaddvars= &eventaddvars,  
     
    compriskaddvars= &compriskaddvars,
    
     
 
    usebetadata =&usebetadataa,
    betadata= &betadataa  ,      
    simuldata=  &simuldataa,     
    resultsdata= &resultsdataa,    
    survdata= &survdataa,       
    outputs = &outputs,   
    print_stats = &print_stats ,
    check_cov_models = &check_cov_models,  
    print_cov_means = &print_cov_means,  
    covmeandata = &covmeandataa ,
    save_raw_covmean = &save_raw_covmean,
    observed_surv = &observed_surva,
    intervname = &intervnamea,
    runnc = &runnca ,
    refint = &refinta,
    seed= &seed,          
    nsamples=&nsamples,     
    nsimul=&nsimul,         
    nparam= &nparam,        
    hazardratio=&hazardratio ,  
    intcomp =   &intcomp ,  /* needs to be a list with two numbers for comparing the hazard ratio when hazardratio = 1 */
    bootstrap_hazard = &bootstrap_hazard,  
    hazardname = &hazardname.a ,  
   
    sample_start = &sample_start ,  
    sample_end = &sample_end,      
    savelib = &savelib,   


    rungraphs = &rungraphs ,
    title1= &title1a,
    title2= &title2a,
    title3= &title3a,
    titledata= &titledataa ,
    graphfile=&graphfilea ,
    tsize=&tsizea ,     
    weight = &weight,
    printlogstats = &printlogstats,
    usespline = &usespline ,
    checkaddvars = &checkaddvars,
    minimalistic = &minimalistic, 
    testing = &testing  
    ) 


%mend;

%macro modelB ;

 data datab ;
 set &data(where = (&arm = 1));
 run;

    %let simuldatab = ;
    %let resultsdatab= ;
    %let survdatab = ;
    %let covmeandatab = ;
    %let observed_survb = ;
    %let intervnameb = ;
    

    %if %bquote(&simuldata)^= %then %let simuldatab = &simuldata.b ;          
    %if %bquote(&resultsdata)^=  %then %let resultsdatab = &resultsdata.b;
    %if %bquote(&survdata)^=  %then %let survdatab = &survdata.b  ;
    %if %bquote(&covmeandata)^=  %then %let covmeandatab = &covmeandata.b ;      
    %if %bquote(&observed_surv)^=  %then %let observed_survb = &observed_surv.b; 
    %if %bquote(&intervname)^=  %then %let intervnameb = &intervname.b ;
    
    
    %if &numintb > 0 %then %do;
        %do i = 1 %to &numintb ;
            %let interv&i = &&intervb&i ;
            %put interv&i for modelB ::: &&interv&i ;
        %end;
    %end;

 %gformula(
    data= datab,           
    id= &id,      
    time= &time,       
    timeptype = &timeptype ,
    timeknots = &timeknots,
    timeinc = &timeinc,
    timefuncgen = &timefuncgen,   
    timepoints= &timepoints,    
    interval=&interval,      

    outc= &outc,          
    outctype=&outctype,  
    outcinteract=&outcinteract,   
    outcwherem = &outcwherem ,  
    outcwherenosim= &outcwherenosim,  
    outcnosimelsemacro = &outcnosimelsemacro,  
    comprisk= &comprisk,        
    compriskinteract= &compriskinteract,   
    compriskwherem = &compriskwherem , 
    compriskwherenosim= &compriskwherenosim,   
    comprisknosimelsemacro = &comprisknosimelsemacro, 
  
     
    numint= &numintb ,
    fixedcov= &fixedcovb,      
    ncov= &ncovb,          

    
    %do i = 1 %to &ncovb ;
        cov&i = &&covb&i ,
        %do j = 1 %to 25 ;
            %let word = %scan(&suffixlist,&j);
            cov&i.&word = &&covb&i.&word ,
        %end;
    %end;
    /*other override options used for each gformula call */

    wherevars= &wherevars,        
    keepsimuldata= &keepsimuldata,     
    equalitiessimuldata= &equalitiessimuldata, 
    eventaddvars= &eventaddvars,  
     
    compriskaddvars= &compriskaddvars,
    
     
 
    usebetadata =&usebetadatab,
    betadata= &betadatab ,      
    simuldata=  &simuldatab,     
    resultsdata= &resultsdatab,    
    survdata= &survdatab,       
    outputs = &outputs,   
    print_stats = &print_stats ,
    check_cov_models = &check_cov_models,  
    print_cov_means = &print_cov_means,  
    covmeandata = &covmeandatab ,
    save_raw_covmean = &save_raw_covmean,
    observed_surv = &observed_survb,
    intervname = &intervnameb,
    runnc = &runncb ,
    refint = &refintb,
    seed= &seed,          
    nsamples=&nsamples,     
    nsimul=&nsimul,         
    nparam= &nparam,        
    hazardratio=&hazardratio ,  
    intcomp =   &intcomp ,  /* needs to be a list with two numbers for comparing the hazard ratio when hazardratio = 1 */
    bootstrap_hazard = &bootstrap_hazard,  
    hazardname = &hazardname.a ,  
   
    sample_start = &sample_start ,  
    sample_end = &sample_end,      
    savelib = &savelib,   


    rungraphs = &rungraphs ,
    title1= &title1b,
    title2= &title2b,
    title3= &title3b,
    titledata= &titledatab ,
    graphfile=&graphfileb ,
    tsize=&tsizeb ,  
    weight = &weight,
    printlogstats = &printlogstats,
    usespline = &usespline ,
    checkaddvars = &checkaddvars,
    minimalistic = &minimalistic, 
    testing = &testing  
    ) 

%mend ;


%macro rescaleround2;
   
%if &outctype=bineofu or &outctype=binsurv %then %do;
   
  /* FOR BINARY OUTCOME, MULTIPLY RISKS BY 100 TO GET % AND ROUND OFF TO TWO DECIMAL PLACES.*/
   

   riskA = round(riskA*10000)/100;
   riakA_mean = round(riskA_mean*10000)/100;
   riskA_std  = round(riskA_std*10000)/100;  
   riskA_ulim95 = round(riskA_ulim95*10000)/100;
   riskA_llim95 = round(riskA_llim95*10000)/100;


   
   riskB = round(riskB*10000)/100;
   riakB_mean = round(riskB_mean*10000)/100;
   riskB_std  = round(riskB_std*10000)/100;  
   riskB_ulim95 = round(riskB_ulim95*10000)/100;
   riskB_llim95 = round(riskB_llim95*10000)/100;
   
   RD = round(RD*10000)/100;

   RD_mean = round(RD_mean*10000)/100;
   RD_std  = round(RD_std*10000)/100;

   
   RD_ulim95 = round(RD_ulim95*10000)/100;
   RD_llim95 = round(RD_llim95*10000)/100;
   

   NNT = round(NNT);

   NNT_mean = round(NNT_mean);
   NNT_std  = round(NNT_std);

   

   
      NNT_ulim95 = round(NNT_ulim95);
      NNT_llim95 = round(NNT_llim95);
   
   
%end;
%else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
   
   obsp= round(&obsp*100)/100;
   call symput('obspm',obsp);
   s&outc = round(s&outc*100)/100;

   s&outc._mean = round(s&outc._mean*100)/100;
   s&outc._std  = round(s&outc._std*100)/100;
   
   s&outc._ulim95 = round(s&outc._ulim95*100)/100;
   s&outc._llim95 = round(s&outc._llim95*100)/100;
   

   RD = round(RD*100)/100;

   RD_mean = round(RD_mean*100)/100;
   RD_std  = round(RD_std*100)/100;
   
   if int^=0 then do;
      RD_ulim95 = round(RD_ulim95*100)/100;
      RD_llim95 = round(RD_llim95*100)/100;
      end;
   
%end;

/* in both cases do the following: */
   RR = round(RR*100)/100;

      RR_ulim95 = round(RR_ulim95*100)/100;
      RR_llim95 = round(RR_llim95*100)/100;
   
   intervened = round(intervened*10000)/100;
   averinterv = round(averinterv*10000)/100;
   
%mend rescaleround;

%macro labels2;
   label int2_A= 'Interventions under A';
   label int2_B= 'Interventions under B';

   label riskA_ulim95 = 'Upper limit 95% CI under A';
   label riskA_llim95 = 'Lower limit 95% CI under A';
   label riskB_ulim95 = 'Upper limit 95% CI under B';
   label riskB_llim95 = 'Lower limit 95% CI under B';
   label RR_ulim95 = 'Upper limit 95% CI';
   label RR_llim95 = 'Lower limit 95% CI';
   label RD_ulim95 = 'Upper limit 95% CI';
   label RD_llim95 = 'Lower limit 95% CI';
   label intervened = '% Intervened On';
   label averinterv = 'Aver % Intervened On';
   label int       = 'Interv.';
   label int2      = 'Description';
   
%if &outctype=binsurv or &outctype=bineofu %then %do;
   label riskA        = 'Risk (%) under A';
   label riskA_std    = 'Bootstrap Risk SE under A';
   label riskA_mean   = 'Bootstrap Risk Mean under A';
   label riskB        = 'Risk (%) under B';
   label riskB_std    = 'Bootstrap Risk SE under B';
   label riskB_mean   = 'Bootstrap Risk Mean under B';
   label rr        = 'Risk ratio';
   label RD        = 'Risk difference';
   label NNT       = '# Needed to Treat';
   label NNT_ulim95 = 'Upper limit 95% CI';
   label NNT_llim95 = 'Lower limit 95% CI';
   %end;

%else %if &outctype=conteofu %then %do;
   label pD        = 'Mean';
   label pD_std    = 'Bootstrap Mean SE';
   label pD_mean   = 'Bootstrap Mean Mean';
   label rr        = 'Ratio of means';
   label RD        = 'Difference of Means';
   %end;
   
%mend labels;

%macro results2 ;
  
    %let numint = %sysfunc(min(&numinta,&numintb));
    %let runnc = %sysfunc(min(&runnca , &runncb)) ;

    %put (numint, runnc)  &numint &runnc ;
    data fin ;  run ;
    data temp ; run ;

     %do i= (1-&runnc) %to &numint;
      
        data intA_&i ;
        set mysurva ;
        where int = &i ;
        keep risk&timepoints int int2 _sample_ ;
        rename risk&timepoints = riskA int2 = int2_A ;
        run;

        data intB_&i ;
        set mysurvb ;
        where int = &i ;
        keep risk&timepoints int int2 _sample_ ;
        rename risk&timepoints = riskB int2 = int2_B ;
        run;
 

        data interv&i; 
        merge intA_&i intB_&i;
        by _sample_; 
        if riskA^=0 then rr=riskB / riskA;
        if riskA^=0 then rd=riskB - riskA ;
        if rd^=. and rd^=0 then nnt = 1/rd;       
        run;


     
          %*Appending intervention datasets;
          data fin; 
          set fin interv&i; 
          if _sample_=0;
          if int=. then int=&i;
          run;
      


          %*Calculating bootstrap mean, variance and confidence intervals;
          proc univariate data=interv&i noprint;
          where _sample_ ne 0;
          var riskA riskB rr rd nnt;
          output out = temp&i
          mean = riskA_mean riskB_mean RR_mean RD_mean NNT_mean
          std =  riskA_std  riskB_std RR_std  RD_std  NNT_std
          pctlpre = riskA_ riskB_ RR_     RD_     NNT_
          pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
          run;

          data temp&i;
          set temp&i;
          int = &i;
          run;

          data temp ;
          set temp temp&i;
     %end;

     data fin ;
     merge fin temp ;
     by int ;
     if int ne . ;
     run;


     %*Cleaning up results to print nicely; 

     data finfin;
     set fin;      
     
     %rescaleround2; /* RESCALE AND ROUND OFF THE OUTPUT */
     %labels2;       /* LABEL THE OUTPUT */
     run;



     %*Outputting results dataset;
     %if %bquote(&resultsdata)^= %then %do;
          %if &printlogstats = 1 %then %put ;
          %if &printlogstats = 1 %then %put  Outputting results to &resultsdata;
          data &resultsdata;
          set finfin;
          run;
     %end;
     %if &printlogstats = 1 %then %put ;

     %*Printing results;

     %if &outctype=binsurv or &outctype=bineofu %then %do;    
          title4 "PREDICTED RISK UNDER SEVERAL INTERVENTIONS UNDER TWO ARMS :  A (&arm = 0)  and B (&arm = 1 )";
     %end;
     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
          title4 "PREDICTED MEAN &outc UNDER SEVERAL INTERVENTIONS UNDER TWO ARMS : A (&arm = 0) and B (&arm = 1)";
     %end;
     proc print data=finfin noobs label double;
     var int  int2_A int2_B ;
     run;

/***
     %let additional_text = ;
     %if &hazardratio = 1 %then %do;
         %let additional_text =, Hazard Ratio between interventions %scan(&intcomp,1) and %scan(&intcomp,2) is &sample_hazard ;
         %if &bootstrap_hazard = 1 %then %let additional_text = &additional_text (&hrlb , &hrub) ;
     %end;
**/
     
     title7 "Data= &data ";
     title8 "Number of bootstrap samples= &nsamples";
     title9 "Reference is using arm A";
/****/

     proc print data=finfin noobs label double; 
     var int riskA  riskA_mean riskA_llim95 riskA_ulim95 riskB riskB_mean  riskB_llim95 riskB_ulim95  /* intervened averinterv */; 
     run;

     proc print data=finfin noobs label double; 
     var int riskA    riskB   rr rr_llim95 rr_ulim95  /* intervened averinterv */; 
     run;

     proc print data=finfin noobs label double; 
     %if &outctype=binsurv or &outctype=bineofu %then %do;
          var int riskA riskB rd rd_llim95 rd_ulim95 nnt nnt_llim95 nnt_ulim95;
     %end;
     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
          var int s&outc s&outc._llim95 s&outc._ulim95 rd rd_llim95 rd_ulim95 ;    
     %end;
     run;


     title;
/***
     %* Deleting no longer needed datasets;

     proc datasets library=work nolist; 
     delete    _paramdata_  _simuldata_ _inputd_ _beta_  _ref_ fin finfin
     %if &outctype = binsurv %then ausc drmst drmst_stat rmst_stat drmst_out  ;
     %if &check_cov_models %then _diff_mean  _cov_std2 
     %do ii = 0 %to &ncov ; _mean_&ii %if &&usevisitp&ii = 1 %then _mean_vp&ii ; %end ;
     ;        
     ;
     quit; 

***/

%mend ;

%create_sample; 

options mprint notes nomlogic nosymbolgen ;




%let interva1 = intno=1 ,  nintvar=1,
    intlabel='All subjects exercise at least 30 minutes per day in all intervals',
    intvar1 = act, inttype1 = 2, intmin1=30, intpr1=1, inttimes1 = 0 1 2 3 4 5 ;





%let intervb1 = intno=1 ,  nintvar=1,
    intlabel='All subjects exercise at least 45 minutes per day in all intervals',
    intvar1 = act, inttype1 = 2, intmin1=45, intpr1=1, inttimes1 = 0 1 2 3 4 5 ;




proc datasets library = work nolist ;
save sample ;
quit;

options notes mprint mprintnest   ;
*options nonotes nomprint ;
%let numint = ; 
%gformula_wrapper(
data= sample,
id=id,
time=time,
timepoints = 6,
outc=dia,
outctype=binsurv,
comprisk =  dead  ,
arm = arm ,

timeptype= concat, 
timeknots = 1 2 3 4 5,

/* for first arm */
ncova=2,
numinta=1 ,
runnca = 1 ,
refinta = 1,
fixedcova = baseage,
cova1  = hbp,    cova1otype  = 2, cova1ptype = tsswitch1,
cova2  = act,    cova2otype  = 4, cova2ptype = lag2cub,

/* for second arm */
ncovb = 2 ,
numintb = 1,
runncb = 1,
refintb = 0,
fixedcovb= baseage ,
covb1  = act,    covb1otype  = 4, covb1ptype = lag2cub,
covb2 = cont, covb2otype = 3 , covb2ptyoe = lag1bin ,
seed= 9458, nsamples = 10 ,

survdata = mysurv,
intervname = myinterv,
resultsdata = myresults
 
);


 