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
    data sample(drop = i j ahbp1-ahbp8 aact1-aact8);

        call streaminit(5027);

        do i=1 to 1000;
            baseage = int( 35 + 25*rand('uniform'));

            array ahbp(8);
            array aact(8);

            do j=1 to 8;
                ahbp(j) = (0.2>rand('uniform'));
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
                

              dia = ( (j/500) >rand('uniform'));
                censlost  = (0.05>rand('uniform'));
                dead      = (0.05>rand('uniform'));
               



                output;

                end;
           end;




    run;


data sample;
set sample; 
if censlost=1 then do;
       dia= .;
       dead= .;
end;
else do;
   if dead=1 then dia= .;
end;
run;

proc means data=sample;
title 'Means of SAMPLE data';
run;
%mend ;

%create_sample;

libname mydata '/proj/users/rwlogan';

proc copy inlib = mydata outlib=work ;
select surv3_mar_d surv4_mar_d ;
run;


**GFORMULA Call;
title 'GFORMULA SAMPLE';
options mprint notes nomlogic nosymbolgen;
options nonotes nomprint ;





proc datasets library = work nolist ;
save surv4_mar_d surv3_mar_d sample ;
quit;

proc sort data = surv3_mar_d out = fortesting (drop = f1) ;
by id t0 ;
run;

data fortesting ;
set fortesting ;
by id ;
A_l1 = lag(A);
L_l1 = lag(L);

A_l2 = lag2(A);
L_l2 = lag2(L);
if first.id = 1 then do;
	a_l1 = A ;
	a_l2 = A ;
	l_l1 = 0 ;
	l_l2 = 0 ;


end;
if t0 = 1 then do ;
	a_l2 = a_l1 ;
	l_l2 = 0 ;
end;
if c = 1 or d = 1 then censor = 1 ;
run;
 
%let todaysDate = %sysfunc(today(), date9.);

*options notes mprint mprintnest ;

%gformula(
data= fortesting,
id=id,
time=t0,
timepoints = 10,

outc=y,
outctype=binsurv,
eventaddvars = L L_l1 A A_l1 ,

comprisk =  ,
compriskaddvars = L L_l1 A A_l1 ,

censor =    ,
censoraddvars = L L_l1 A A_l1 ,
maxipw = 1000 ,
compevent =  ,

fixedcov = ,
timeptype= conbin, 
timeknots = 1 2 3 4 5 6 7 8 9 ,

ncov=2,
cov1  = L,    cov1otype  = 1, cov1ptype = lag2bin, cov1mtype=nocheck, cov1addvars =  L_l1 A_l1,
cov2  = A,    cov2otype  = 3, cov2ptype = lag2bin, cov2mtype=nocheck, cov2addvars =  L L_l1 A_l1 ,

hazardratio = 0 ,
intcomp = 0 1 ,
seed= 9458, nsamples = 10, numint=0 , 
check_cov_models = 1 ,
save_raw_covmean = 1,
print_cov_means = 1,
rungraphs = 1 ,

 graphfile= /proj/sas_macros/gformula/master_branch/new_censoring/GFORMULA-SAS/scenario0_&todaysDate..pdf,
checkaddvars = 0 
);

/******

data obsrisk ;
set _observed_surv_all (keep = obrisk: );
array aobrisk{10} obrisk0-obrisk9 ;
do i = 1 to 10 ;
	t0 = i - 1 ;
	obrisk = aobrisk[i];
	output ;
end;
keep t0 obrisk ;
run ;

data gform;
set  survdata_all (keep = risk:);
array arisk{10} risk1 - risk10 ;
do i = 1 to 10 ;
     t0 = i - 1 ;
	 grisk = arisk[i];
	 output ;
end;
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


*****/

