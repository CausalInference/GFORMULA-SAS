/* Adapted from example_censor.sas (CausalInference/GFORMULA-SAS):
   the longitudinal-data preparation that precedes the %gformula() call.
   The original imports the full surv4_MAR_D.csv (a many-thousand-row
   survival dataset); here we read a small real sample of the same file
   shipped under input/ (header and columns unchanged). The lag/lag2,
   first.id baseline carry-forward, NA->missing recoding, and INPUT()
   conversions are verbatim. A PROC MEANS shows the prepared data. */

proc import datafile="input/surv4_sample.csv"
        out=want
        dbms=csv
        replace;
     getnames=yes;
run;

data want1; set want;
A_l1=lag(A);
L_l1=lag(L);
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

Y_=Y;
D_=D;
if Y_="NA" then Y_=.;
if D_="NA" then D_=.;
drop Y D ;
run;

data bytimes; set want1;
Y=input(Y_,8.);
D=input(D_,8.);
drop  Y_ D_ var1;
run;

proc means data=bytimes;
title 'Prepared longitudinal data (bytimes)';
run;
