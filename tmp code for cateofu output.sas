/****
 proc print data=finfin noobs label double; 
	     var int &pd &pd._llim95 &pd._ulim95 rr rr_llim95 rr_ulim95 &pd._mean &pd._std intervened averinterv; 
	     run;

	     proc print data=finfin noobs label double; 
	     %if &outctype=binsurv or &outctype=bineofu %then %do;
	          var int &pd &pd._llim95 &pd._ulim95 rd rd_llim95 rd_ulim95 nnt nnt_llim95 nnt_ulim95;
	     %end;
	     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 or &outctype = cateofu %then %do;
	          var int s&outc s&outc._llim95 s&outc._ulim95 rd rd_llim95 rd_ulim95 ;    
	     %end;
	     run;

***/

/*****

It may be helpful to update the results printout for the "cateofu" outcome type. 

What do you think of the following proposal?
For each J=j (where j represents a level that the categorical / ordinal outcome can take from J = 1 to J = J), print the following:

Observed proportion for level j =  __%
   A table that contains the following information for each simulated intervention

     Proportion for level j       Upper limit 95% CI       Lower limit 95% CI Ratio for level j  Upper limit 95% CI       Lower limit 95% CI

           e.g., the ratio of the proportion with J= j comparing intervention 1 vs. intervention 0

    Difference for level j       Upper limit 95% CI       Lower limit 95% CI

          e.g., the difference of the proportion with J= j comparing intervention 1 vs. intervention 0

	Bootstrap proportion SE for level j

	Bootstrap proportion mean for level j

	% intervened on
		Identical for all levels j
	Aver % intervened 
		Identical for all levels j

In addition, print the usual information that applies to all levels j:
PREDICTED PROPORTIONS UNDER SEVERAL INTERVENTIONS
Data =  __, Sample size = __, Monte Carlo sample size = __
Number of bootstrap samples = __
Reference intervention is __

****/

%macro silly ;

%let pd = smybinom ;
%let outc = mybinom ;
data test1(keep = int &outc &pd &pd._llim95 &pd._ulim95 rd rd_llim95 rd_ulim95 nnt nnt_llim95 nnt_ulim95 
                     rr rr_llim95 rr_ulim95 &pd._mean &pd._std intervened averinterv)  ;
set finfin;
label &pd ="Proportion (%)"  &outc = "&outc level"
         &pd._llim95= "Lower limit 95% CI"
		 &pd._ulim95= "Upper limit 95% CI"
		 &pd._mean="Proportion mean"
		 &pd._std= "Proportion SE"

;
%do ii = 1 %to 5 ;
     &outc = &ii ;
     &pd = &pd._&ii ;
	 &pd._mean = &pd._&ii._mean ;
	 &pd._std = &pd._&ii._std ;
	 &pd._llim95 = &pd._&ii._llim95 ;
	 &pd._ulim95 = &pd._&ii._ulim95 ;
	 rd = rd_&ii ;
	 rd_llim95 = rd_&ii._llim95 ;
	 rd_ulim95 = rd_&ii._ulim95 ;
	 rr = rr_&ii ;
	 rr_llim95 = rr_&ii._llim95 ;
	 rr_ulim95 = rr_&ii._ulim95 ;
	 nnt = nnt_&ii ;
	 nnt_llim95 = nnt_&ii._llim95 ;
	 nnt_ulim95 = nnt_&ii._ulim95 ;
     output ;
%end;
run;


ods pdf file="c:\users\roger\cateofu_output.pdf";

proc print data=test1 noobs label double; 
	     var int &outc &pd &pd._llim95 &pd._ulim95 rr rr_llim95 rr_ulim95 &pd._mean &pd._std intervened averinterv; 
	     run;


proc print data = test1 noobs label double  ;
var int &outc &pd &pd._llim95 &pd._ulim95 rd rd_llim95 rd_ulim95 nnt nnt_llim95 nnt_ulim95;
run ;

ods pdf close ;

%mend ;

%silly ;


    
