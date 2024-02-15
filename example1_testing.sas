/*Example: outctype=binsurv*/
%let macropath = C:\Users\Roger\OneDrive - Harvard University\Gformula work\GitHub files\GFORMULA-SAS ; 
%include "&macropath\gformula4.0.sas";
options linesize=88 pagesize=54;


libname mytmp 'C:\' ;


*options mprint mprintnest;
*options mlogic mlogicnest;
*options fullstimer notes;
options notes;
options nonotes nomprint ;
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

proc freq data = sample ;
table 
time ;
run;


data sample ;
set sample ;
conteof = . ;
call streaminit(1234321);
if time = 5 then do ;
   conteof = -1+2*rand('uniform');
   bineof = rand('bernoulli',0.4);
end;
if time = 3 then do;
   act = act_l1 ;
   act_l1 = act_l2;
end;
if time = 4 then act_l1 = act_l2 ;

censor = (dead = 1 or censlost = 1) ;
run;

proc datasets library = work nolist ;
save sample ;
quit;

**GFORMULA Call;
title 'GFORMULA SAMPLE';
options mprint notes ;
*options nomprint nonotes ;
options mlogic symbolgen ;
options nomlogic nosymbolgen ;
%gformula(
data= sample,
id=id,
time=time,
timepoints = 6,
outc=bineof ,
outctype=bineofu ,
compevent=,
compevent_cens  = 0   ,
censor = censor ,

usehistory_eof = 1 ,


fixedcov = baseage,
timeptype= concat, 
timeknots = 1 2 3 4 5,

ncov=2,
cov1  = hbp,    cov1otype  = 2, cov1ptype = tsswitch1, cov1etype = tsswitch1 none ,
cov2  = act,    cov2otype  = 4, cov2ptype = skpcub, cov2skip = 3 , cov2etype = skpzqdc all , cov2knots  = 10 20 30   , 

hazardratio = 0 ,
intcomp = 0 1 ,
seed= 9458, nsamples = 0, numint=1 ,
rungraphs = 0 ,
usespline = 0 ,
testing = no ,
simuldata =  ,
savelib = work 
);


