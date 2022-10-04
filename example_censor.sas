%include 'gformula4.0.sas';


options linesize=88 pagesize=54;


options notes mprint;

proc import datafile="surv4_MAR_D.csv"
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

%let interv1 = intno=1, nintvar=1,
    intlabel='All subjects maintain A to at least 0.18 in all intervals',
    intvar1 = A, inttype1 = 2, intmin1=0.18, intpr1=1, inttimes1 = 0 1 2 3 4 5 6 7 8 9 10 ;





*Scenario 1: When there are competing events (competing events are censoring event) and censoring due to loss to follow up, 


%gformula(
data= bytimes,
id=id,
time=t0,
timepoints = 10,

outc=y,
outctype=binsurv,
compevent = d ,
maxipw = 1000 ,
censor=c,
compevent_cens  = 1   , 

fixedcov = ,
timeptype= conbin, 
timeknots = 1 2 3 4 5 6 7 8 9 ,

ncov=2,
cov1  = L,    cov1otype  = 1, cov1ptype = lag2bin,  
cov2  = A,    cov2otype  = 3, cov2ptype = lag2bin,  

hazardratio = 0 ,
intcomp = 0 1 ,
seed= 9458, 
numint=1 , 
check_cov_models = 0 ,
save_raw_covmean = 0,
print_cov_means = 0,
rungraphs = 1 ,
graphfile=test.pdf ,
checkaddvars = 0 ,
nsamples = 20
);



*Scenario 2: When there are competing events (competing events are not considered as censoring events) and censoring due to loss to follow up, 


%gformula(
data= bytimes,
id=id,
time=t0,
timepoints = 10,

outc=y,
outctype=binsurv,
compevent = d ,
maxipw = 1000 ,
censor=c,
compevent_cens  = 0  , 

fixedcov = ,
timeptype= conbin, 
timeknots = 1 2 3 4 5 6 7 8 9 ,

ncov=2,
cov1  = L,    cov1otype  = 1, cov1ptype = lag2bin,  
cov2  = A,    cov2otype  = 3, cov2ptype = lag2bin,  

hazardratio = 0 ,
intcomp = 0 1 ,
seed= 9458, 
numint=1 , 
check_cov_models = 0 ,
save_raw_covmean = 0,
print_cov_means = 0,
rungraphs = 1 ,
graphfile=test.pdf ,
checkaddvars = 0 ,
nsamples = 20
);



*Scenario 3;
*when there is no competing events, parameter compevent can be left empty.
*when there is complete follow up,  parameter censor can be left empty.
    
    

    
    

    
    
    
    


 





