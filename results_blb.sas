
%macro results_blb (blb_samples = 1 );

     data fin_s;
	 run;

	 %do sample = 1 %to &blb_samples ;
     %*Generating reference cumulative incidence;        
     data _ref_;
     set interv&refint._all (where = (_sample_s = &sample   )); /*change jgy*/
     %if &outctype=binsurv %then   pDref=pD;
     %else pDref=s&outc ;;
     keep _sample_r pDref;
     run;



     %do i=&intstart %to &numint;

          %*Outputting summary of intervention;
          %*proc means data = interv&i._all;
          %*run;

          %*Comparing intervention to reference;
                 
          %let pd = pd ;
          %if &outctype ^= binsurv %then %let pd = s&outc ;

          data interv&i; 
          merge interv&i._all (where = (_sample_s = &sample )) _ref_;
          by _sample_r; 
          if pDref^=0 then rr=&pd /pDref;
          if pDref^=0 then rd=&pd - pDref;
          if rd^=. and rd>0 then nnt = 1/rd;
          *logrr=log(rr); /* commented out since this was the only occurrance of this variable */
          run;

		  /***************

          %*Appending intervention datasets;
		  * this data set will be used for printing the output. has rows for each intervention ;

          data fin; 
          set fin interv&i; 
          if _sample =0  ;
          if int=. then int=&i;
          run;
            
		  *********/


          %*Calculating bootstrap mean, variance and confidence intervals;
          proc univariate data=interv&i noprint;
          where _sample_r ne 0;
          var &pd rr rd nnt;
          output out = temp&i
          mean = &pd._mean RR_mean RD_mean NNT_mean
          std =  &pd._std  RR_std  RD_std  NNT_std
          pctlpre = &pd._  RR_     RD_     NNT_
          pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
          run;

          data temp&i;
          set temp&i;
          int = &i;
		  _sample_s = &sample ;
          run;


		  /*** used for printing output 
          data fin;
          merge fin temp&i;
          by int;
          run;
          *************/
 
		  data fin_s ;
		  set fin_s temp&i ;
		  if _sample_s ne . ;
		  run;


          %*Deleting no longer needed datasets;
        *  proc datasets library=work nolist; 
        *  delete  interv&i interv&i._all temp&i;
        *  quit;

     %end;


     %if &hazardratio = 1 AND &bootstrap_hazard = 1 %then %do;
             proc univariate data = &hazardname (where = (_sample_ > 0)) ;
             var hazardratio ;
             output out = _inthrstat_ 
             mean = hr_mean  
             std =  hr_std   
             pctlpre = hr_   
             pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
             run;

             proc sql ;
             select round(hr_llim95,0.01) into :hrlb from _inthrstat_ ;
             select round(hr_ulim95,0.01) into :hrub from _inthrstat_;
             quit;

             %let hrub = %sysfunc(compress(&hrub));
             %let hrlb = %sysfunc(compress(&hrlb));

     %end;
%end;


proc sort data = fin_s ; by int _sample_s ;run;



%let mylist = pd_mean rr_mean rd_mean nnt_mean 
              pd_std  rr_std  rd_std  nnt_std 
              pd_llim95 rr_llim95 rd_llim95 nnt_llim95
			  pd_ulim95 rr_ulim95 rd_ulim95 nnt_ulim95 ;

proc means data = fin_s noprint ;
var &mylist ;
by int ;
output out=fin_results mean(&mylist ) = ;
run;

data fin ;
run;

data _ref_ ;
set interv&refint._all (where = (_sample_r = 0 and _sample_s = 0));
%if &outctype=binsurv %then   pDref=pD;
%else pDref=s&outc ;;
keep _sample_s _sample_r pDref;
run ;


  %do i=&intstart %to &numint;

          %*Outputting summary of intervention;
          %*proc means data = interv&i._all;
          %*run;

          %*Comparing intervention to reference;
                 
          %let pd = pd ;
          %if &outctype ^= binsurv %then %let pd = s&outc ;

          data interv_base&i; 
          merge interv&i._all (where = (_sample_s = 0 and _sample_r = 0)) _ref_;
          by _sample_s _sample_r; 
          if pDref^=0 then rr=&pd /pDref;
          if pDref^=0 then rd=&pd - pDref;
          if rd^=. and rd>0 then nnt = 1/rd;
          *logrr=log(rr); /* commented out since this was the only occurrance of this variable */
          run;


          %*Appending intervention datasets;
          data fin; 
          set fin interv_base&i; 
          if _sample_r=0 and _sample_s = 0;
          if int=. then int=&i;
          run;

		  
         
     %end;


 	data fin;
    merge fin fin_results;
    by int;
    run;

     %*Cleaning up results to print nicely; 

     data finfin;
     set fin;      
     %rescaleround; /* RESCALE AND ROUND OFF THE OUTPUT */
     %labels;       /* LABEL THE OUTPUT */
     run;


	 
     %*Outputting results dataset;
     %if %bquote(&resultsdata)^= %then %do;
          %if &printlogstats = 1 %then %put ;
          %if &printlogstats = 1 %then %put  Outputting results to &resultsdata;
          data &resultsdata;
          set finfin;
          run;
     %end;
     %if &printlogstats = 1 %then %put ;

     %*Printing results;

     %if &outctype=binsurv or &outctype=bineofu %then %do;    
          title4 'PREDICTED RISK UNDER SEVERAL INTERVENTIONS';
     %end;
     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
          title4 "PREDICTED MEAN &outc UNDER SEVERAL INTERVENTIONS";
     %end;
     proc print data=finfin noobs label double;
     var int int2;
     run;


     %let additional_text = ;
     %if &hazardratio = 1 %then %do;
         %let additional_text =, Hazard Ratio between interventions %scan(&intcomp,1) and %scan(&intcomp,2) is &sample_hazard ;
         %if &bootstrap_hazard = 1 %then %let additional_text = &additional_text (&hrlb , &hrub) ;
     %end;

     %if &outctype=binsurv or &outctype=bineofu %then %do;
	     %if %bquote(&censor) ^= %then %do ;
		     title6 "IP-weighted natural course risk= %sysevalf(&obspm) % &additional_text "; 	
		 %end;
		 %else %do;
             title6 "Observed risk= %sysevalf(&obspm) % &additional_text ";
		 %end;
     %end;      
     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
          title6 "Observed mean= %sysevalf(&obspm) ";
     %end;
     title7 "Data= &data, Sample size= &ssize, Monte Carlo sample size= &BLB_b";
     title8 "Number of bootstrap samples using Bag of Little Bootstraps method with &BLB_s samples of size &BLB_b and &BLB_r samples of size &ssize";
     title9 "Reference intervention is &refint";

     proc print data=finfin noobs label double; 
     var int &pd &pd._llim95 &pd._ulim95 rr rr_llim95 rr_ulim95 &pd._mean &pd._std intervened averinterv; 
     run;

     proc print data=finfin noobs label double; 
     %if &outctype=binsurv or &outctype=bineofu %then %do;
          var int &pd &pd._llim95 &pd._ulim95 rd rd_llim95 rd_ulim95 nnt nnt_llim95 nnt_ulim95;
     %end;
     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
          var int s&outc s&outc._llim95 s&outc._ulim95 rd rd_llim95 rd_ulim95 ;    
     %end;
     run;


%mend ;
