/* Adapted from example1.sas (CausalInference/GFORMULA-SAS):
   the %create_sample data-generation macro and its PROC MEANS, run
   standalone. The %include of gformula4.0.sas and the %gformula() call
   that follow it upstream are omitted so this exercises just the sample
   construction. Logic is verbatim from example1.sas. */

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
