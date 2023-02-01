%macro blb_graphs_part1(blb_samples = 1 );
  /*   %if (&check_cov_models = 1 OR &rungraphs  ) AND &runnc = 1 %then %do; */


          %if &save_raw_covmean = 1  %then %do;

               data &covmeanname._raw ;
               set &covmeanname ;
               run;

          %end;


		  data res_s ;
		  run;

		 

           /* include base sample for _sample_s = 0 */
		  %do sample = 0 %to &blb_samples ;


	          %if &outctype = binsurv %then %do;
					/* results from observed data, used for observed_risk data : cuminc and sumsurv  */
			        /* for use in natural course, int = 0 */
 
	               data rr ;              
	               set &covmeanname (keep = _sample_r _sample_s  _time_   meanoutc   meancompevent where = (_sample_s = &sample) ) ;
	               by _sample_r _time_ ;
	               retain cumsurv  cuminc newn ;
	               if first._sample_r then do ;
	                    cumsurv = 1;
	                    cuminc = 0 ;
	                    newn = 0 ;
	               end;
	               newn = newn + 1 ;                            
	               inc = cumsurv * meanoutc * (1 - meancompevent) ;
	               cuminc = cuminc + inc;                
	               surv = (1 - meanoutc) * ( 1 - meancompevent ) ;
	               cumsurv = cumsurv * surv;
	               if newn <= &timepoints then output;
	               keep _sample_r  _time_  cuminc cumsurv ;
	               run;
		 

	               proc transpose data = rr  out=_stmp_ prefix= obsurv ;
	               var cumsurv ;
	               id _time_ ;
	               by _sample_r; 
	               run;

	               proc transpose data = rr out = _itmp_ prefix= obrisk ;
	               var cuminc ;
	               id _time_ ;
	               by _sample_r ;
	               run;

	               data _observed_surv_s  ;
	               merge _stmp_ _itmp_ ;
	               by _sample_r ;
	               _int_ = 0 ;
				   _sample_s = &sample ;
	               drop _NAME_ ;
	               run;

				   %if &sample = 0 %then %do;
				   		data &observed_surv ;
						set _observed_surv_s ;
						run;
					%end;
					%else %do;
						data &observed_surv ;
						set &observed_surv _observed_surv_s ;
						run;
				   %end;
	               proc datasets library = work nolist ;
	               delete rr _stmp_ _itmp_ ;
	               quit;

	          %end; /* binsurv */
	          %else %do; /* for all other outcome types  */
	               data _obsoutcmean_ &observed_surv ;
	               set &covmeanname (keep = _sample_r _sample_s &time &outc where = (_sample_s = &sample ) );
	               if &time = &timepoints -1  ;  
	               run;
	 
	          %end;

	          %do i = 1 %to &ncov ;

	               proc transpose data = &covmeanname  (where = (_sample_s = &sample))  out = _mean_&i(drop = _NAME_ ) prefix = &&cov&i ;
	               var &&cov&i ;
	               by _sample_r ;
	               run;

	               %if &&usevisitp&i = 1 %then %do;
	                    proc transpose data = &covmeanname (where = (_sample_s = &sample )) out = _mean_vp&i (drop = _NAME_) prefix=&&cov&i.randomvisitp ;
	                    var &&cov&i.randomvisitp;
	                    by _sample_r;
	                    run;

	               %end;
	          %end ;

	          data &covmeanname._s ;
	          merge %do i = 1 %to &ncov ;
	                    _mean_&i 
	                    %if &&usevisitp&i = 1 %then _mean_vp&i ;
	               %end;
	               ;
	          by _sample_r ;
	          run;

	          data _diff_mean ;
	          merge &covmeanname._s  (keep = _sample_r 
	               %do i = 1 %to &ncov ;
	                    %do k = 1 %to &timepoints ;
	                         &&cov&i..&k
	                         %if &&usevisitp&i = 1 %then &&cov&i.randomvisitp.&k ;
	                    %end;
	               %end;)                
	               interv0_all(keep = _sample_r _sample_s 
	               %do i = 1 %to &ncov ;
	                    %do k = 1 %to %eval(&timepoints);
	                       /***  s&&cov&i..&k ***/
						      ncs&&cov&i..&k 
	                         %if &&usevisitp&i = 1 %then %do ;
							    /***	s&&cov&i.randomvisitp.&k  ***/
							     ncs&&cov&i.randomvisitp.&k
	                         %end;     
	                    %end;
	               %end; 
				
					  rename = (
						%do i = 1 %to &ncov ;
		                    %do k = 1 %to %eval(&timepoints);	                          
							   ncs&&cov&i..&k = s&&cov&i..&k 
							   %if &&usevisitp&i = 1 %then  ncs&&cov&i.randomvisitp.&k = s&&cov&i.randomvisitp.&k ;  
		                    %end;               				  
				   		%end; 
						 )
						 where = (_sample_s = &sample )
				);

	           %do i = 1 %to &ncov ;
	                    %do k = 1 %to &timepoints ;
	                         d&&cov&i..&k  = &&cov&i..&k  - s&&cov&i..&k ;
	                         %if &&usevisitp&i = 1 %then d&&cov&i.randomvisitp.&k = &&cov&i.randomvisitp.&k - s&&cov&i.randomvisitp.&k;;
	                    %end;
	           %end;
			   by _sample_r ;
			   
	           run;

			  %if &sample = 0 %then %do;
			  		data _diff_mean_0 ;
					set _diff_mean ;				
					run;
			  %end;
			  %else %do;

		          proc univariate data = _diff_mean(where = (_sample_r > 0))  noprint  ;
		          var %do i = 1 %to &ncov ;
		                    %do k = 1 %to &timepoints  ;
		                         d&&cov&i..&k 
		                         %if &&usevisitp&i = 1 %then d&&cov&i.randomvisitp.&k ;
		                    %end;
		               %end;
		               ;

		          output out = _cov_std2 
		               std = %do i = 1 %to &ncov ;
		                         %do k = 1 %to &timepoints  ;
		                              d&&cov&i..&k._stddev 
		                              %if &&usevisitp&i = 1 %then d&&cov&i.randomvisitp.&k._stddev ;
		                         %end;
		                    %end;
		               pctlpre = %do i = 1 %to &ncov ;
		                              %do k = 1 %to &timepoints  ;
		                                   d&&cov&i..&k 
		                                   %if &&usevisitp&i = 1 %then d&&cov&i.randomvisitp.&k ;
		                              %end;
		                         %end;
		               pctlname = _pct025 _pct975 
		               pctlpts = 2.5 97.5 ;
		          run;


			  data _cov_std2 ;
			  set _cov_std2 ;
			  _sample_s = &sample ;
			  run;


			data res_s ;
			set res_s  _cov_std2 ;
			if _sample_s > . ;
			run;
	        %end;

 		%end; /* sample loop */




		%let myvarlist = ;
		%do i = 1 %to &ncov ;
			%do j = 1 %to &timepoints ;
				%let myvarlist = &myvarlist d&&cov&i..&j._pct025 d&&cov&i..&j._pct975 d&&cov&i..&j._stddev ;
				%if &&usevisitp&i = 1 %then %do;
					%let myvarlist = &myvarlist d&&cov&i.randomvisitp.&j._pct025 d&&cov&i.randomvisitp.&j._pct975 d&&cov&i.randomvisitp.&j._stddev ; 
				%end;
			 %end;
		%end;


		
	proc means data = res_s noprint ;
	var &myvarlist ;
	output out=cov_results mean(&myvarlist ) = ;
	run;


		  /*****************/

          data &covmeanname._test ;
          merge _diff_mean_0  cov_results  ; 
          drop  _sample_r ;
          label %do i = 1 %to &ncov ;
                    &&cov&i = "mean of observed &&cov&i" 
                    s&&cov&i = "mean of simulated &&cov&i"
                    &&cov&i.._diff = "difference of mean observed and simulated"
                    &&cov&i.._stddev = "Standard deviation of difference"
                    &&cov&i.._lb = "95% lower bound for difference"
                    &&cov&i.._ub = "95% upper bound for difference"
                    &&cov&i.._lbp="2.5 percentile for difference"
                    &&cov&i.._ubp="97.5 percentile for difference"
                    %if &&usevisitp&i = 1 %then %do;
                         &&cov&i.randomvisitp = "mean of observed &&cov&i.randomvisitp" 
                         s&&cov&i.randomvisitp = "mean of simulated &&cov&i"
                         &&cov&i.randomvisitp._diff = "difference of mean observed and simulated"
                         &&cov&i.randomvisitp._stddev = "Standard deviation of difference"
                         &&cov&i.randomvisitp._lb = "95% lower bound for difference"
                         &&cov&i.randomvisitp._ub = "95% upper bound for difference"
                         &&cov&i.randomvisitp._lbp="2.5 percentile for difference"
                         &&cov&i.randomvisitp._ubp="97.5 percentile for difference"
                    %end;
               %end ;
               &time = "n-th time point" ;
          %do k = 1 %to &timepoints  ;
                %let j = &k ;
                &time = %eval( &k - 1) ;
                %do i = 1 %to &ncov ;
                         &&cov&i = &&cov&i..&k ;
                         s&&cov&i = s&&cov&i..&j;
                         &&cov&i.._diff = d&&cov&i..&k ;
                         &&cov&i.._stddev = d&&cov&i..&k._stddev ;
                         &&cov&i.._lb = &&cov&i.._diff - 1.96 * &&cov&i.._stddev ;
                         &&cov&i.._ub = &&cov&i.._diff + 1.96 * &&cov&i.._stddev ;
                         &&cov&i.._lbp = d&&cov&i..&k._pct025 ;
                         &&cov&i.._ubp = d&&cov&i..&k._pct975 ;
                         %if &&usevisitp&i = 1 %then %do;
                              &&cov&i.randomvisitp = &&cov&i.randomvisitp.&k ;
                              s&&cov&i.randomvisitp = s&&cov&i.randomvisitp.&j;
                              &&cov&i.randomvisitp._diff = d&&cov&i.randomvisitp.&k ;
                              &&cov&i.randomvisitp._stddev = d&&cov&i.randomvisitp.&k._stddev ;
                              &&cov&i.randomvisitp._lb = &&cov&i.randomvisitp._diff - 1.96 * &&cov&i.randomvisitp._stddev ;
                              &&cov&i.randomvisitp._ub = &&cov&i.randomvisitp._diff + 1.96 * &&cov&i.randomvisitp._stddev ;
                              &&cov&i.randomvisitp._lbp = d&&cov&i.randomvisitp.&k._pct025 ;
                              &&cov&i.randomvisitp._ubp = d&&cov&i.randomvisitp.&k._pct975 ;
                         %end;
                         drop &&cov&i..&k s&&cov&i..&j d&&cov&i..&k d&&cov&i..&k._stddev  d&&cov&i..&k._pct025  d&&cov&i..&k._pct975
                              %if &&usevisitp&i = 1 %then &&cov&i.randomvisitp.&k s&&cov&i.randomvisitp.&j d&&cov&i.randomvisitp.&k
                              d&&cov&i.randomvisitp.&k._stddev  d&&cov&i..&k._pct025  d&&cov&i..&k._pct975 ;
                         ;
                 %end;
                    output ;
            %end;
            run;
            /********************/


            %if &print_cov_means = 1 %then %do ;
                    %do i = 1 %to &ncov ;
                         proc print data = &covmeanname._test noobs label;
                         title "Comparison of means of observed &&cov&i and simulated &&cov&i over &time ";
                         var &time   &&cov&i   s&&cov&i  &&cov&i.._diff   &&cov&i.._stddev  &&cov&i.._lb    &&cov&i.._ub  &&cov&i.._lbp &&cov&i.._ubp           ;               
                         run;    
                         title ;

                         %if &&usevisitp&i = 1 %then %do;
                              proc print data = &covmeanname._test noobs label;
                              title "Comparison of means of observed &&cov&i.randomvisitp and simulated &&cov&i.randomvisitp over &time ";
                              var &time   &&cov&i.randomvisitp   s&&cov&i.randomvisitp  &&cov&i.randomvisitp._diff   &&cov&i.randomvisitp._stddev 
                                   &&cov&i.randomvisitp._lb    &&cov&i.randomvisitp._ub  &&cov&i.randomvisitp._lbp &&cov&i.randomvisitp._ubp           ;               
                              run;    
                              title ;

                         %end;
                    %end;
           %end;


           %if &rungraphs = 1 %then %do;
		            %let useboot = 0 ;
					%if &bootstrap_method = 0 and &nsamples > 0 %then %let useboot = 1 ;
					%if &bootstrap_method = 1 and &BLB_s > 0 and &BLB_r > 0 %then %let useboot = 1 ;
 
                    
                    %construct_graphs(
                         time=&time ,
                         outcome=&outc ,
                         compevent = &compevent ,
                         outctype = &outctype,
                         covmean=&covmeandata._test ,
                         obssurv = &observed_surv ,
                         simsurv = &survdata ,
                         sixgraphs = 0 ,
                         gfilename= &graphfile ,
                         title1= &title1 ,
                         title2= &title2,
                         title3= &title3, 
                         titledata=  &titledata , 
                         tsize=&tsize ,
                         frombootstrap = &useboot ) ;
          %end;

		 
   /*  %end; */

%mend ;
