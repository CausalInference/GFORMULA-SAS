
%macro results_blb (blb_samples = 1 );
/* blb_samples is BLB_s, the number of outer samples of size BLB_b used in
   bootstrap samples */
/*****/

%if (&check_cov_models = 1 OR &rungraphs = 1  ) AND &runnc = 1 %then %do;
    /* following macro includes the first parts of the results macro */
    
	%blb_graphs_part1(blb_samples = &blb_samples  ) ;
%end;

/******/
     data fin_s;
	 run;

	 data _inthrstat_s ;
	 run;

/* for each sample of size BLB_b of BLB_s samples (outer limit ) use one 
	 loop of the results macro and then average over the BLB_s subsets */

	 %do sample = 1 %to &blb_samples ;

	     %*Generating reference cumulative incidence;        
	     data _ref_s ;
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
	          merge interv&i._all (where = (_sample_s = &sample )) _ref_s ;
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
	             proc univariate data = &hazardname (where = (_sample_s = &sample   ))  ;
	             var hazardratio ;
	             output out = temp_inthrstat_ 
	             mean = hr_mean  
	             std =  hr_std   
	             pctlpre = hr_   
	             pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
	             run;

				data temp_inthrstat_ ;
				set temp_inthrstat_ ;
				_sample_s = &sample ;
				run;

				data _inthrstat_s ;
				set _inthrstat_s temp_inthrstat_ ;
				if _sample_s ne . ;
				run;



	     %end;
%end;


proc sort data = fin_s ; 
by int _sample_s ;
run;



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



	   %if &hazardratio = 1 AND &bootstrap_hazard = 1 %then %do;

	            %let mylist2 = hr_mean hr_std hr_llim95 hr_ulim95 ;
	   
				proc means data = _inthrstat_s noprint ;
				var &mylist2 ;
			*	by int ;
				output out=_inthrstat_ mean(&mylist2 ) = ;
				run;

	             proc sql ;
	             select round(hr_llim95,0.01) into :hrlb from _inthrstat_ ;
	             select round(hr_ulim95,0.01) into :hrub from _inthrstat_ ;
	             quit;

	             %let hrub = %sysfunc(compress(&hrub));
	             %let hrlb = %sysfunc(compress(&hrlb));
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
     title7 "Data= &data, Sample size= &ssize, Monte Carlo sample size= &nsimul ";
	 %if &use_disjoint_blb_samples = 0 %then %do;
     	title8 "Number of bootstrap samples using Bag of Little Bootstraps method with &BLB_s samples of size &BLB_b, each with &BLB_r samples of size &nsimul";
	 %end;
	 %else %if &use_disjoint_blb_samples = 1 %then %do;
     	title8 "Number of bootstrap samples using Bag of Little Bootstraps method with &BLB_s disjoint samples of size &BLB_b, each with &BLB_r samples of size &nsimul";
	 %end;
	 %if &bootstrap_method = 2 %then %do;
			title8 "Number of bootstrap samples using Bag of Little Bootstraps method using &BLB_s samples of size &BLB_b with adaptive selection of r using &bsample_counter bootstrap samples"; 
	 %end;
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


/******/


%macro blb_pct_helper(datain = , varlistin =  , dataout= );

proc sort data = &datain out = _tmp_ ;
by _sample_s _sample_r ;
run ;



proc univariate data = _tmp_ noprint ;
var &varlistin ;
output out = _tmp2_ pctlpre = &varlistin 
          pctlname = _pct025 _pct975 
          pctlpts = 2.5 97.5 ;
by _sample_s ;
run;


%let varlisttmp = ;

%let ntmp = %numargs(&varlistin ) ;


%do itmp = 1 %to &ntmp ;
	%let word = %scan(&varlistin , &itmp ) ;
	%let varlisttmp = &varlisttmp &word._pct025 &word._pct975 ;
%end;

proc means data = _tmp2_ noprint ;
var &varlisttmp ;
output out = &dataout ( keep = &varlisttmp  ) mean( &varlisttmp  ) = ;
run;

%mend ;

%macro check_rconvergence (sample = 1) ;
/*
sample = the sth collection of bootstrap samples being tested for convergence. Data holds all previously run samples 

*/

data test_s ; run ;

     %*Generating reference cumulative incidence;        
	     data _ref_s ;
	     set interv&refint._all (where = (_sample_s = &sample   )); /*change jgy*/
	     %if &outctype=binsurv %then   pDref=pD;
	     %else pDref=s&outc ;;
	     keep _sample_r pDref;
	     run;

		%local i rindex ;

		
	                 
	    %let pd = pd ;
	    %if &outctype ^= binsurv %then %let pd = s&outc ;


	     %do i = 1 %to &BLB_r_trend ;

	       
	          data rsubset&i; 
	          set interv&refint._all (where = (_sample_s = &sample and _sample_r <= %eval(&rend - (&i - 1) )) ) ;
			  keep &pd  _sample_r _sample_s ;
	          run;

	          %*Calculating bootstrap mean, variance and confidence intervals;
	          proc univariate data=rsubset&i noprint;
	          where _sample_r ne 0;
	          var &pd ;
	          output out = temp&i	 
			  mean = &pd._mean 
	          std =  &pd._std   
	          pctlpre = &pd._ 
	          pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
	          run;
				
			  data temp&i;
			  set temp&i ;
			  test_sample = %eval(&rend - (&i - 1)) ;
			  run;

			  %if &i = 1 %then %do;
			  	data test_s ;
				set temp&i;
				run;
			  %end;
			  %else %do;
			  	data test_s ;
			  	set test_s temp&i ;			  
			  	run;
			  %end;
     	%end;

  


	proc sort data = test_s ;
	by descending test_sample ;
	run;

	data test_s ;
	set test_s ;
	retain ref_llimit ref_ulimit ref_mean ref_std ;
	if test_sample = &rend then do ;
		ref_llimit = pd_llim95 ;
		ref_ulimit = pd_ulim95 ;
		ref_mean   = pd_mean;
		ref_std    = pd_std ;
	end;
	test_llimit = abs((pd_llim95 - ref_llimit)/ref_llimit) ;
	test_ulimit = abs((pd_ulim95 - ref_ulimit)/ref_ulimit) ;
	test_mean = abs((pd_mean - ref_mean)/ref_mean );
	test_std = abs((pd_std - ref_std )/ref_std );
	run;

	%if &BLB_r_test_method = 1 %then %do;
		proc means data = test_s noprint;
		var test_llimit test_ulimit ;
		output out = test_s2 sum(test_llimit test_ulimit)= ;
		run;

		data test_s2 ;
		set test_s2 ;
		conv_check = test_llimit + test_ulimit ;
		if conv_check > &BLB_r_epsilon then converged  = 0 ;
		else converged = 1 ;
		call symput('rconverged',compress(converged));
		rename _freq_ = trend ;	
		put _all_ ;
		run;
	%end;
	%else %if &BLB_r_test_method = 2 %then %do;
		data test_s ;
		set test_s end = _end_ ;
		retain maxcheck  ;
        if _n_ = 1 then do;
            
			 maxcheck = -1 ;
		end;
		check = 0.5 * (test_llimit + test_ulimit ) ;
		maxcheck = max(maxcheck,check);
        
		if _end_ then do ;
		    if maxcheck > &BLB_r_epsilon then converged = 0 ;
			else converged = 1 ;
		    put maxcheck= converged= ;
			call symput('rconverged',compress(converged));
		end;
		run;
	%end;
	%else %if &BLB_r_test_method = 3 %then %do;
		data test_s ;
		set test_s end = _end_ ;
		retain maxcheck  ;
        if _n_ = 1 then do;
            
			 maxcheck = -1 ;
		end;
		check = 0.5 * (test_mean + test_std ) ;
		maxcheck = max(maxcheck,check);
        
		if _end_ then do ;
		    if maxcheck > &BLB_r_epsilon then converged = 0 ;
			else converged = 1 ;
		    put maxcheck= converged= ;
			call symput('rconverged',compress(converged));
		end;
		run;
	%end;

/*	%let rconverged = 1 ; */



%mend ;
