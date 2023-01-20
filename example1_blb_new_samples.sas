/*Example: outctype=binsurv*/

%include '/proj/sas_macros/gformula/master_branch/new_censoring/GFORMULA-SAS/gformula4_0_blb_testing.sas';
%include '/proj/sas_macros/gformula/master_branch/new_censoring/GFORMULA-SAS/bootstrap_blb_version0.sas';
%include '/proj/sas_macros/gformula/master_branch/new_censoring/GFORMULA-SAS/base_sample.sas';
%include '/proj/sas_macros/gformula/master_branch/new_censoring/GFORMULA-SAS/results_blb.sas';
%include '/proj/sas_macros/gformula/master_branch/new_censoring/GFORMULA-SAS/for_covmeans_graphs.sas';

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

/***/
proc datasets library = work nolist ;
save sample  /* roger1 simul1 param1 step2a */;
quit;
/***/

**GFORMULA Call;
title 'GFORMULA SAMPLE';
options mprint mprintnest notes spool ;
options nomlogic ;
*options nonotes nomprint ;

%let use_samples_orig = 1 ;
%let sample_check = -1 ;

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
cov2  = act,    cov2otype  = 4, cov2ptype = lag2cub,

hazardratio = 0 , /* this is broken */
bootstrap_hazard = 1 ,
intcomp = 0 1 ,
seed= 9458, numint=1,
nsimul = 1000 ,
bootstrap_method =1,
nsamples = 10 ,
BLB_b = 500 ,
BLB_r = 2, /* for inner loop, takes place of nsamples */
BLB_s = 2 /* for outter loop */,
rungraphs = 1,
print_cov_means = 0 
);


