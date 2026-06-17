/* Adapted from example1_testing.sas (CausalInference/GFORMULA-SAS):
   the extended %create_sample macro (adds x and z covariate arrays via
   rand('normal') and rand('binomial')) followed by PROC FREQ on time, a
   second DATA step adding end-of-followup outcomes via rand('bernoulli'),
   rand('binomial') and rand('uniform'), and PROC DATASETS SAVE. The
   %include and trailing %gformula() call are omitted. Deterministic via
   seeds 5027 and 1234321. */

%macro create_sample ;
%let condition = ;
%let condition = dia or censlost or dead ;

**SAMPLE Data;
    data sample(drop = i j ahbp1-ahbp8 aact1-aact8 ax1-ax8 az1-az8 );

        call streaminit(5027);

        do i=1 to 1000;
            baseage = int( 35 + 25*rand('uniform'));

            array ahbp(8);
            array aact(8);
			array ax{8} ;
			array az{8};

            do j=1 to 8;
                ahbp(j) = (0.2>rand('uniform'));
                if j > 1 & ahbp(j-1) = 1 then ahbp(j) = 1 ;

        		aact(j)=(0.7>rand('uniform'));
        		if aact(j)=1 then do;
                   aact(j) = int(exp(3.5+0.4*(rand('normal'))));
        		end;
				ax(j ) = 2* rand('normal') ;
				az(j) = rand('binomial',0.4,5) + 1 ;
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

				x = ax(j);
                x_l1 = ax(j-1);
				x_l2 = ax(j-2);

				z = az(j);
				z_l1 = az(j-1);
				z_l2 = az(j - 2);

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
table time ;
run;

data sample ;
set sample ;
conteof = . ;
call streaminit(1234321);
if time = 5 then do ;
   conteof = -1+2*rand('uniform');
   bineof = rand('bernoulli',0.4);
   myconteof4 = rand('bernoulli',0.7);
   if myconteof4 = 1 then myconteof4 = 1 + 5*abs(rand('normal')) ;
   mybinom = rand('binomial',0.5,4) + 1 ;
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
