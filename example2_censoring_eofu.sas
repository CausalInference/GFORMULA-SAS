
options linesize=88 pagesize=54;
options nonotes;

 

%include '/proj/sas_macros/gformula/master_branch/new_censoring/GFORMULA-SAS/gformula3.sas';
options linesize=88 pagesize=54;


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
                

                if time < 9 then censlost  = (0.01>rand('uniform'));
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

title 'GFORMULA SAMPLE';

options nomprint nonotes ;
options notes mprint ; 

%let myoutc = cont_e ;

%create_sample(event = &myoutc ) ;

 

proc datasets library = work nolist ;
save  sample ;
run;
quit;




%gformula(
    data= sample,
    id=id,
    time=time,
    timepoints = 10, 

    outc= &myoutc ,
    outctype= conteofu   ,
    outcinteract =   ,
    comprisk =   , 
	censor = censlost ,
	compevent =   ,
	maxipw = p99 , 
    fixedcov =  baseage  , /* using hbp and bmi for fixedcov forces the corresponding baseline variables into each model */
    ncov=2,
    timeptype=concat, 
    timeknots= 1 2 3 4 5 6 7 8 9, 
    cov1 = hbp,    cov1otype = 2, cov1ptype = lag1bin    ,   
    cov2 = bmi,    cov2otype = 3, cov2ptype = lag2cub  ,  
    seed= 9458,
    check_cov_models = 1 ,
    print_cov_means = 0,
    save_raw_covmean = 0,
    /* datasets */
     
    nsamples = 0, 
    
    numint=0 ,
    rungraphs = 0,
    graphfile=bin_e.pdf ,
     printlogstats = 1 
    );



data obsrisk ;
set _observed_surv_all (keep= &myoutc );
obrisk = &myoutc ;
t0 = 9 ; 
keep t0 obrisk ;
run ;

data gform;
set  survdata_all (keep = s&myoutc);
grisk = s&myoutc ;
t0 = 9 ;
keep t0 grisk ;
run;

data both ;
merge gform obsrisk ;
by t0 ;
reldiff = abs(grisk - obrisk)/grisk ;
run;

proc print data = both ;
format obrisk grisk f16.8 ;
var t0 obrisk grisk  reldiff ;
run;



