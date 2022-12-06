/*Example: outctype=binsurv*/

%include '/proj/sas_macros/gformula/master_branch/new_censoring/GFORMULA-SAS/gformula4.0.sas';
options linesize=88 pagesize=54;

*options mprint mprintnest;
*options mlogic mlogicnest;
*options fullstimer notes;
options notes;
*options threads;



%let interv1 = intno=1 ,  nintvar=1, intvisittype = 1 ,
    intlabel='All subjects exercise at least 10 minutes: with intvisittype= 1',
    intvar1 = act, inttype1 = 2, intmin1=10, intpr1=1, inttimes1 = 0 1 2 3 4 5 ;
%let interv2 = intno=2 ,  nintvar=1, intvisittype = 2 ,
    intlabel='All subjects exercise at least 10 minutes: with intvisittype= 2',
    intvar1 = act, inttype1 = 2, intmin1=10, intpr1=1, inttimes1 = 0 1 2 3 4 5 ;



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



data sample ;
set sample (rename = (act = act_orig));
call streaminit(1234);
*if time = 1 or time = 3 then act = act_l1 ;
visit_act = rand('bernoulli',0.4);
if time = 0 then visit_act = 1 ;
retain ts_last_act ;
if visit_act = 1 then ts_last_act = 0;
else ts_last_act = ts_last_act + 1 ;
if ts_last_act = 3 then do ;
    visit_act = 1 ;
	ts_last_act = 0 ;
end;
retain act ;
if visit_act = 1 then act = act_orig ;
drop act_l1 act_l2 ;
run;

data sample ;
set sample ;
act_l1 = lag(act);
act_l2 = lag2(act );
ts_last_act_l1 = lag(ts_last_act);

if time = 0 then do ;
    act_l1 = act ;
	act_l2 = act ;
	ts_last_act_l1 = . ;
end;
if time = 1 then act_l2 = act_l1 ;
drop act_orig ;
run;

proc freq data = sample ;
tables time * visit_act ts_last_act ;
run;


**GFORMULA Call;
title 'GFORMULA SAMPLE';
options mprint notes ;
options nomprint nonotes ;
%gformula(
data= sample,
id=id,
time=time,
timepoints = 6,
outc=dia,
outctype=binsurv,
compevent=dead,
compevent_cens  = 0   ,
censor = censlost ,

fixedcov = baseage,
timeptype= concat, 
timeknots = 1 2 3 4 5,

ncov=2,
cov1  = hbp,    cov1otype  = 2, cov1ptype = tsswitch1,
cov2  = act,    cov2otype  = 3, cov2ptype =lag1bin, 
   cov2randomvisitp = visit_act, cov2visitpmaxgap = 2 , cov2visitpcount = ts_last_act ,

hazardratio = 0 ,
intcomp = 0 1 ,
seed= 9458, nsamples = 10, numint=2 , usespline=0 
);


