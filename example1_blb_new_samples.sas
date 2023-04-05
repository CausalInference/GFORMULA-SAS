/*Example: outctype=binsurv*/

%include '/proj/sas_macros/gformula/master_branch/new_censoring/GFORMULA-SAS/gformula4.1.sas';

options linesize=88 pagesize=54;

*options mprint mprintnest;
*options mlogic mlogicnest;
*options fullstimer notes;
options notes;
*options threads;



%let interv1 = intno=1 ,  nintvar=1,
    intlabel='All subjects exercise at least 30 minutes per day in all intervals',
    intvar1 = act, inttype1 = 2, intmin1=30, intpr1=1, inttimes1 = 0 1 2 3 4 5 ;



%macro create_sample ;
%let condition = ;
%let condition = dia or censlost or dead ;
 




**SAMPLE Data;
    data sample(drop = i j ahbp1-ahbp18 aact1-aact18);

        call streaminit(5027);

        do i=1 to 10000;
            baseage = int( 35 + 25*rand('uniform'));

            array ahbp(18);
            array aact(18);

            do j=1 to 18;
                ahbp(j) = (0.2>rand('uniform'));
                if j > 1 & ahbp(j-1) = 1 then ahbp(j) = 1 ;

        aact(j)=(0.7>rand('uniform'));
        if aact(j)=1 then do;
                   aact(j) = int(exp(3.5+0.4*(rand('normal'))));
        end;
                end;

            do j=3 to 18  until ( &condition   ) ;
                id=i;
                time=j-3;
                
                hbp     = ahbp(j);
                hbp_l1  = ahbp(j-1);
                hbp_l2  = ahbp(j-2);
                
                act     = aact(j);
                act_l1  = aact(j-1);
                act_l2  = aact(j-2);
                
              dia = ( (j/500) >rand('uniform'));
                censlost  = (0.05>rand('uniform'));
                dead      = (0.05>rand('uniform'));              
                output;

                end;
           end;
    run;


data sample;
set sample (rename = (act = act_orig)); 
call streaminit(1) ;
if censlost=1 then do;
       dia= .;
       dead= .;
end;
else do;
   if dead=1 then dia= .;
end;

visit_act = rand('bernoulli',0.4);
retain act ts_last_act;
if time = 0 then do ;
   visit_act = 1 ;
end;
if visit_act = 1 then do ;
	act = act_orig ;
	ts_last_act = 0 ;
end;
else do ;
	ts_last_act = ts_last_act + 1 ;
end;
if ts_last_act = 3 then do ;
	visit_act = 1 ;
	act = act_orig ;
	ts_last_act = 0 ;
end;

*drop act_l1 act_l2 ; 
run;

data sample ;
set sample ;

act_holder0 = lag(act);
act_holder1 = lag2(act);

if time ge 1 then act_l1 = act_holder0 ;
if time ge 2 then act_l2 = act_holder1 ;

drop act_holder0 act_holder1 ;
run;

proc means data=sample;
title 'Means of SAMPLE data';
run;
%mend ;

%create_sample;

/***/
proc datasets library = work nolist ;
save sample  /* roger1 simul1 param1 step2a */;
quit;
/***/

*proc greplay igout = GSEG nofs ;
* delete _all_ ;
* run;
* quit;

proc freq data = sample  ;
tables time * visit_act ;
run;


options ls = 120 ;
**GFORMULA Call;
title 'GFORMULA SAMPLE';
options mprint mprintnest notes spool mlogic ;
options nomlogic ;
options nonotes nomprint ;

%let use_samples_orig = 1 ;
%let sample_check = -1 ;
%let BLB_s_start = 1 ; 

%gformula(
data= sample,
id=id,
time=time,
timepoints = 16,
outc=dia,
outctype=binsurv,
compevent=dead,
compevent_cens  = 0   ,
censor = censlost ,

fixedcov = baseage,
timeptype= concat, 
timeknots = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15  ,

ncov=2,
cov1  = hbp,    cov1otype  = 2, cov1ptype = tsswitch1,
cov2  = act,    cov2otype  = 4, cov2ptype = lag2cub, cov2randomvisitp = visit_act , cov2visitpmaxgap=2 , cov2visitpcount =ts_last_act , 

hazardratio = 0 , 
bootstrap_hazard = 0 ,
intcomp = 0 1 ,
seed= 9458, numint=0,
bootstrap_method =2,
nsamples = 2 ,
BLB_b = 1000 ,
BLB_s = 5 /* for outter loop */,
BLB_r = 5, /* for inner loop, takes place of nsamples */

BLB_s_start = 5 ,
BLB_s_max = 10 ,
BLB_s_delta = 3 ,
BLB_s_trend = 3 ,
BLB_s_epsilon = 0.01,
BLB_s_test_method  = 4,


BLB_r_start = 10 ,
BLB_r_max = 30 ,
BLB_r_delta = 10 ,
BLB_r_trend = 15 ,
BLB_r_epsilon = 0.01,
BLB_r_test_method  = 4 ,
BLB_use_seeds = ,
rungraphs = 0,
print_cov_means = 0,
printlogstats = 0,
expand_param_counts = 1 ,
expand_simul_counts = 1  
);


