
%macro samples_blb(sample_s = , sample_r = );

/* want to create samples for innter loop of samples (blb_s, blb_r) */
/* not sure where to do this. */
    %local subseed i ;

    %************ CREATING BOOTSTRAP SAMPLES;  
    
       /* %if &printlogstats = 1 %then */ %put  Creating bootstrap sample sample_s = &sample_s , sample_r = &sample_r , seed = &seed , ssize=&ssize , nsimul = &nsimul , nparam = &nparam  ;
       /* %if &printlogstats = 1 %then %put    ; */

        %*Generating random sample of ids to be used in bootstrap sample; 

		data _idsamples ;
		set step2(where = (sample_si = &sample_s and sample_rj = &sample_r )) ;
		keep newid BLB_count0 ;
		run;

        %*Merging with parameter estimation data;
   
      * paramdata is original data for parameter models. paramsample is bootstrap sample with newid , numberhits, and the current bootstrap sample > 0
      this overrites the value in paramdata where newid and time agree. ;
   
        data _paramsample_ ;
        set _idsamples ;
        _sample_r = &sample_r  ;
        _sample_s = &sample_s ;
        run; 

  
    * add in the variable numberhits for the number of times a subject is selected into the bootstrap sample ;
    data param ;
        merge _paramdata_ (in= p) _paramsample_;
        by newid ;
        if BLB_count0 > 0 ; *delete those not selected into sample ;
		%if &use_bootstrap_counts = 0 %then %do;
	        do _copy_ = 1 to BLB_count0 ; * make numberhits copies of each remaining subject ;
			    BLB_counts = 1 ;
				%if &bootstrap_method = 0 %then bootstrap_counts = 1 ;;
	            output ;
	        end;
		%end;
        drop %if &use_bootstrap_counts = 0 %then BLB_count0 ; 
             %if %bquote(&censor) = %then _copy_ ; ;
        run;


	%if %bquote(&censor) ^= %then %do;
		proc sort data = param ;
		by newid _copy_ &time ;
		run;

		data param ;
		set param (drop = newid _copy_ );
		retain newid ;
		if _n_ = 1 then newid = 0 ;
		if &time = 0 then newid = newid + 1 ;
	    run;
	%end;


       * reset the outcome and covariate bounds to that models and simulated 
        values depend on what would be the observed bounds ;

        data _covbounds_ /*_null_ */ ;
        set param end=_end_;
        %do i=0 %to &ncov;
            %if &&cov&i.otype=3 or &&cov&i.otype=4 or &&cov&i.otype=6 or &&cov&i.otype=7   or &&cov&i.otype = -1  %then %do;
                
                retain &&cov&i.._min &&cov&i.._max;

                if _n_ = 1 then do;
                    &&cov&i.._min =  1.0e100;
                    &&cov&i.._max = -1.0e100;
                    end; 

                *if &&cov&i ne . and &&cov&i < &&cov&i.._min then &&cov&i.._min = &&cov&i;
                /* 6-2015 rwl need to treat otype = 4 - 0-lin different. want min of non-zero values, otherwise min will be 0 */
                %if &&cov&i.otype ^= 4 %then if &&cov&i ne . and &&cov&i < &&cov&i.._min then &&cov&i.._min = &&cov&i;;
                %if &&cov&i.otype  = 4 %then if &&cov&i > 0  and &&cov&i < &&cov&i.._min then &&cov&i.._min = &&cov&i;;
                if &&cov&i ne . and &&cov&i > &&cov&i.._max then &&cov&i.._max = &&cov&i;
                
                if _end_ then do; 

                    * for truncated normal model we will extend the bounds based on the oringal sample (rwl 7/2013) ; 
                  
                      

                    call symput("cov&i.min", trim(left(&&cov&i.._min))); 
                    call symput("cov&i.max", trim(left(&&cov&i.._max)));
                 end;
                
                %end;
            %end;
        %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3  %then %do;

             retain &outc._min &outc._max;

                if _n_ = 1 then do;
                    &outc._min =  1.0e100;
                    &outc._max = -1.0e100;
                    end; 

                if &outc ne . and &outc < &outc._min then &outc._min = &outc;
                if &outc ne . and &outc > &outc._max then &outc._max = &outc;
                
                if _end_ then do;
                    call symput("outcmin",trim(left(&outc._min)) );
                    call symput("outcmax",trim(left(&outc._max)) );
                    end;
        %end;
		%if &bootstrap_method = 0 %then %do ;
			_sample_ = &bsample ;
		%end;
		%else %if &bootstrap_method = 1 %then %do;
         	_sample_r = &sample_r ;
		 	_sample_s = &sample_s ;
		%end;
        keep %if &bootstrap_method = 0 %then _sample_ ;
             %else %if &bootstrap_method = 1 %then _sample_r _sample_s ;
              %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3  %then  &outc._min &outc._max ;
                %do i = 0 %to &ncov ; 
                    %if &&cov&i.otype=3 or &&cov&i.otype=4 or &&cov&i.otype=6 or &&cov&i.otype=7   or &&cov&i.otype = -1  %then &&cov&i.._min &&cov&i.._max ;
                %end ;
                 ;
        if _end_ then output ;
    run;

    %do i = 0 %to &ncov;
        %if &&cov&i.otype=3 or &&cov&i.otype=4 or  &&cov&i.otype=6 or &&cov&i.otype=7  or &&cov&i.otype=-1 %then %do;
            %if &printlogstats = 1 %then %put  bootstrap sample &sample_s &sample_r &&cov&i range is &&cov&i.min to &&cov&i.max;                    
            %end;
        %end;
    %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu %then %do;
      %if &printlogstats = 1 %then %put  bootstrap sample &sample_s &sample_r &outc range is &outcmin to &outcmax;
      %end;
   
        

        %*Merging with simulation data;

        * default simul data set is to take all subjects in param data set, number of subjects in nparam ;
  
    * _simuldata_ has one observation per person for time = 0 ;
        /* idsamples containes number of times newid selected into current param data set */
        data simul ;
        merge _simuldata_  _idsamples;
        by newid; 
		if BLB_count0 > 0 ;
		%if &bootstrap_method = 0 %then %do;
			_sample_ = &bsample ;
		%end;
		%else %if &bootstrap_method = 1 %then %do;
        	_sample_r = &sample_r  ; 
			_sample_s = &sample_s ;
		%end;
        run;


        %if &use_bootstrap_counts = 0 %then %do;
        	data simul ;
        	set simul ;
        	do _copy0_ = 1 to BLB_count0;
			    
			    %if &bootstrap_method = 0 %then bootstrap_counts = 1 ;;
				%if &bootstrap_method = 1 %then BLB_counts = 1 ;;
                output ;
        	end;
        	drop BLB_count0 _copy0_ ;
        	run;
		%end;

        /* want to sample when nsimul is different from size of param data set */
        %if &nsimul ne &ssize %then %do;
            * now take a random sample of size nsimul  when nsimul ne nparam ;
            %let subseed = %eval(2*&seed);

			%if &use_bootstrap_counts = 1 %then %do;
        		data simul_expanded ;
        		set simul ;
			 
        		do _copy0_ = 1 to BLB_count0  ;
                	output ;
        		end;
        		drop BLB_count0 _copy0_ ;
        		run;
			%end;

            proc surveyselect data = %if &use_bootstrap_counts = 0 %then simul ;
			                         %else %if &use_bootstrap_counts = 1 %then sumul_expanded ;
                        out = simul_tmp noprint
                        method=urs  
                        sampsize=&nsimul seed=&subseed;			
            run;

			%if &use_bootstrap_counts = 0 %then %do;
            	data simul ;
            	set simul_tmp ;
				%if &bootstrap_method = 1 %then %do;
            		_sample_r = &sample_r ;
					_sample_s = &sample_s ;
				 %end;
				 %else %if &bootstrap_method = 0 %then %do;
					_sample_ = &bsample ;
				 %end;
             	do _copy_ = 1 to numberhits ;
				    BLB_counts = 1 ;
                    output ;
             	end;
				bootstrap_counts = 1 ;
             	drop _copy_ numberhits ;
             	run;
			%end;
			%else %if &use_bootstrap_counts = 1 %then %do;

				proc sql ;
				create table mytest3 as select
				 id, count(newid) as id_count ,
				         sum(numberhits) as BLB_counts_new 
				 from simul_tmp
				 group by newid ;
				 quit;

 				proc sort data = mytest3   ;
 				by newid;
 				run;

				data simul ;
 				merge simul mytest3  ;
 				by newid ;
				BLB_counts = BLB_counts_new ;
				drop BLB_counts_new ;
 				run;				

			%end;
         %end;

 
   
     %if &hazardratio = 1  AND &bootstrap_hazard = 0  %then %do;     
         data _calchazard_ ;
         calchazard = 0 ;
         _sample_r = &sample_r  ;
		 _sample_s = &sample_s ;
         run;
    %end;

    proc datasets library=work nolist ;
        delete _idsamples ;
    quit;
   
    proc contents data = simul ; run ;
     
%mend samples_blb;

%macro bootstrap_blb ;
	/* variables for BLB method:
         BLB_b = size of smaller sample of original data 
		 BLB_s = number of samples of sze BLB_b to create  : used for outer loop 
		 BLB_r = number of samples of size n to use for each sample of size b : used in inner loop

      The outer loop will be for BLB_s for each sample of size b, then an inner loop over 
      BLB_r for each sample of size n

     In this case sample_start and sample_end are for BLB_s
  ****/

       /* as a first test we will generate ids for all samples used in BLB method */ 
		

		%if &sample_start = 0 %then %let sample_start = 1 ;

        %if &BLB_s > 0  AND  &BLB_r > 0 %then %do; 

			proc surveyselect data = tmpids(keep =newid ) out = step1 (rename = (replicate = sample_si))
	   		method = srs 
	   		sampsize = &BLB_b 
	   		reps = &BLB_s 
			seed = &seed 
	   		noprint 
	   			;
			run ;

			proc surveyselect data = step1 
	  		out=step2 (drop = ExpectedHits SamplingWeight rename = (replicate = sample_rj 
	                                                                numberhits = BLB_count0))
	  		method = urs 
	  		sampsize = &ssize 
	  		reps = &BLB_r 
			seed = %eval(&seed * &seed )
	  		noprint 
	 			;
	  		strata sample_si ;
			run;

	 data step2b ;
	 set step2 ;
	 run;

    %*Looping over bootstrap samples;

    %let bsample = 0 ;
    %do sample_s = 1 %to &BLB_s;
        %if (&outputs ^= yes or %eval(&sample_s) ^= 0) %then %do;
               %let ods_logit = ods select none ;
               %let ods_reg = ods select none ;
               %let ods_qlim = ods select none ;
        %end;

		%do sample_r = 1 %to &BLB_r ;
			%let bsample = %eval(&bsample + 1 );
        		%*Generating samples;

				%if &printlogstats = 1 %then %put  before sample =  ( &sample_s , &sample_r )   seed = &seed ;
				%if &printlogstats = 1 %then %put  ;

		        %samples_blb(sample_s = &sample_s , sample_r = &sample_r ); 


			%if &sample_r = &sample_check and &sample_s = 1  %then %do;
				data param2 ; set param ; run ;
				data simul2 ; set simul ; run ;
			 %end; 
		        %*Estimating parameters;
		        %if &usebetadata = 0 %then %do;            
		                %if &uselabelc = 0  and &uselabelo = 0 %then options nolabel ;;
		                %parameters;     
		                options &label_orig ;  
						%if &sample_r = &sample_check and &sample_s = 1  %then %do;
							data param2 ; set param ; run ;
							data simul2 ; set simul ; run ;
						
			 			%end;  
		        %end;

		        %else %do; /*feb2013*/

		            data _beta_;
		            set &betadata;
		            run;

		        %end;   

 
	        ods select all ; 

	        data _betar_ ;
	        set _beta_  (where = ( _sample_r = &sample_r and _sample_s = &sample_s )) ;
	        run;

	        data _seedr_ ;
	        _seedr_ = %eval(&seed);
	        _sample_r = &sample_r ;
			_sample_s = &sample_s ;
	        run;

       		%if &created_dataviews = 0  %then %do;
	             /* initialize data views for interventions */

	             %if &runnc = 1 %then %do;
	                 %interv_init(intno=0, intlabel='Natural course' ); 
	             %end;   
	        

	             %do intnum = 1 %to &numint;
	                %interv_init(&&interv&intnum);
	             %end;             
 
				 %let created_dataviews = 1 ;
        	%end;            
       
	        %if &hazardratio = 1 %then %do;
	             %if &sample_s = 0 %then  %createhazard ;
	             %else %if &bootstrap_hazard = 1 %then %createhazard ;
	        %end;

        	%if &runnc = 1 %then %do;  /*** for natural course results ***/
            	%*No intervention cumulative incidence;

            	%if &sample_s = 0 AND  %bquote(&simuldata)^=  %then %do;

	                %*Outputting/creating  simulated dataset;

	                data &simuldata ;
	                set simulated0 ;
	                run;
              
	                %if &testing = yes %then %do;
	                     %do intcount = 0 %to &numint ;
	                        data &savelib..simuldata&intcount  ;
	                        set simulated&intcount ;
	                        run;
	                      %end;
	                %end;

	                %*Outputting the mean of covariates and probability of event;
	                proc means data= &simuldata  mean min max ;                                                                      
	                var %if &outctype = binsurv %then cuminc ;
	                    intervened averinterv  
	                    %if &outctype = binsurv %then %do;
	                        %do j = 1 %to %eval(&timepoints);
	                            s&outc.&j
	                        %end;
	                    %end;
	                    %else &outc ;
	                    %if &minimalistic = no %then %do;
	                        %do i = 1 %to &ncov;
	                            %do j = 1 %to %eval(&timepoints);
	                                %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j ;
	                                s&&cov&i..&j
									/**  %if &intno = 0  AND %bquote(&censor)^= %then %do ; ***/
										ncs&&cov&i..&j %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
									/** %end;  ****/
	                            %end;
	                        %end; 
	                    %end;
	                    %if &outctype = binsurv %then %do;
	                        %do n = 1 %to %eval(&timepoints) ;
	                            cumincr&n cumsurvr&n cumcompeventr&n   
	                        %end;
	                    %end;
	                ;
            
                    output out=interv0 mean= %if &outctype = binsurv %then pd  ;
                        intervened averinterv 
                        %if &outctype = binsurv %then %do;
                            %do j = 1 %to %eval(&timepoints);
                                s&outc.&j
                            %end;
                        %end;
                        %else s&outc ;
                        %if &minimalistic = no %then %do;
                            %do i = 1 %to &ncov;
                                %do j = 1 %to %eval(&timepoints);
                                    %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j ;
                                    s&&cov&i..&j
									/** %if &intno = 0  AND %bquote(&censor)^=   %then  %do;  ***/
										ncs&&cov&i..&j %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
									/** %end ;  **/
                                %end;
                            %end;  
                        %end;
                        %if &outctype = binsurv %then %do;
                            %do n = 1 %to %eval(&timepoints) ;
                                cumincr&n cumsurvr&n cumcompeventr&n  
                            %end;
                        %end;
                        ;
                    
                    title "mean, min, max under intervention 0";
	                run;


	                data interv0;
	                set interv0;
	                _sample_r = &sample_r;
					_sample_s = &sample_s ;
	                length int2 $70 ;
	                int=0;
	                int2="Natural course";
	                n=_FREQ_;
	                keep int int2 _sample_ n   intervened averinterv 
	                %if &outctype=binsurv %then %do;
	                    pd
	                    %do j = 1 %to %eval(&timepoints);
	                        s&outc.&j
	                    %end;
	                %end;
	                %else s&outc ;
	                 %if &minimalistic = no %then %do;
	                    %do j = 1 %to %eval(&timepoints);
	                        %do i = 1 %to &ncov;
	                            %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j ;
	                            s&&cov&i..&j
								/*** %if &intno = 0  AND %bquote(&censor)^=   %then  %do; ***/
									ncs&&cov&i..&j %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
								/*** %end; ****/
	                        %end;
	                    %end;  
	                %end;
	                dataname ssize obsp 
	                ;
	                dataname = "&data";
	                    ssize = &ssize ;
	                    obsp = &obsp ;
	                run;             

                %if &outctype = binsurv %then %do;

                   proc means data = &simuldata noprint ;
                   var  %do n = 1 %to %eval(&timepoints);
                        cumincr&n cumsurvr&n cumcompeventr&n 
                    %end; ;
                   output out = survprobs0 mean = %do n = 1 %to &timepoints; risk&n surv&n compevent&n  %end; ;
                   run;

                    /* initialize survdata */

                    data surv_tmp0 ;
                        set survprobs0 ;
                        length int2 $70 ;
                        int = 0 ;
                        int2 = "Natural course";
                        _sample_r = &sample_r;
						_sample_s = &sample_s ;
                        surv0 = 1;
                        n = _freq_ ;
                        keep  int int2 _sample_r _sample_s n  surv0
                        %do n = 1 %to &timepoints;
                            risk&n surv&n compevent&n  
                        %end;
                    ;
                    run;



					data interv0 ;
					merge interv0 surv_tmp0 (keep = surv1 - surv&timepoints );
					array surv{&timepoints } ;
					%do myi = 1 %to &ncov ;
						array ncs&&cov&myi {&timepoints } ;
						%if &&usevisitp&myi = 1 %then array ncs&&cov&myi.randomvisitp { &timepoints } ; ;
					%end;

    				do j = 2 to &timepoints ; 
	    				%do myi = 1 %to &ncov ;
							ncs&&cov&myi [j ] = ncs&&cov&myi [ j ] / surv[j - 1 ] ; * should this be j or j-1;
							%if &&usevisitp&myi = 1 %then ncs&&cov&myi.randomvisitp [ j ] = ncs&&cov&myi.randomvisitp [ j ] / surv[j - 1] ;;
						%end;
     				end;
	 				drop j ;
     				run;
                %end;   
                %else %do;

                        /* initialize survdata */

                    data surv_tmp0 ;
                        set interv0 ;
                        length int2 $70 ;
                        int = 0 ;
                        int2 = "Natural course";
                        _sample_r = &sample_r;
						_sample_s = &sample_s ;                                        
                        keep  int int2 _sample_r _sample_s  s&outc  ;
                    run;
                %end; 
            %end;  /* end of saving simuldata , needed outside of interv and interv_init macros */
            %else %do;                   
                %interv (intno=0 , intlabel='Natural course'); /* do not need to save simulated data set */                                 
            %end;   

         
                %if &sample_s = &sample_start AND &sample_r = 1 AND &chunked = 1  %then %do;  
                    data &survdata ;
                    set  surv_tmp0 ;
                    run;
                %end;
                %else %do ;
                    proc append base = &survdata data = surv_tmp0;
                    run;
                %end;
           
        %end ; /* run natural course = 1 **/
        %*Looping over interventions;


        %do intnum = 1 %to &numint;
            %interv(&&interv&intnum);                    
                %if &intnum = 1 AND &runnc = 0 AND &sample_s = &sample_start AND &sample_r = 1 AND &chunked = 1 %then %do;
                      data &survdata ;
                      set surv_tmp&intnum ;
                      run;
                %end;
                %else %do;
                    proc append base = &survdata data = surv_tmp&intnum ;
                    run;
                %end;            
        %end;

        %*Outputing intervention results into temporary data sets ;

        %let intstart = 0 ;
        %if &runnc = 0 %then %let intstart = 1;

        %if %eval(&sample_s) = &sample_start and &sample_r = 1 and &chunked = 1 %then %do;  
           
            %do int = &intstart %to %eval(&numint);
                
                data interv&int._all;
                    set interv&int;
                run;
            %end;
        %end;
        %else %do;
            %do int = &intstart %to %eval(&numint);
                proc append base = interv&int._all data = interv&int;
                run;
            %end;
        %end;

        %if &chunked = 1 %then %do;
             /* save all interventions into a permanent data set for each chunk */
             %if &sample_s = &sample_start and &sample_r = 1 and &chunked = 1  %then %do; 

                  /* save interventions */
                  %if &runnc = 1 %then %do;
                     data &intervdata ; 
                     set interv0 ;
                     run;
                  %end;
                  %else %do;
                      data &intervdata ;
                      set interv1; 
                      run;
                  %end;
                   
                  %do int = %eval(&intstart +1) %to &numint ;
                      *proc append base = &intervdata    data = interv&int ;
                      *run;
		      data &intervdata ;
		      set &intervdata interv&int ;
		      run;
                  %end;  
             %end;
             %else %do;
                  %do int = &intstart %to &numint ;
                      *proc append base = &intervdata    data = interv&int  ;
                      *run;
		      data &intervdata ;
		      set &intervdata interv&int ;
		      run;
                  %end; 
             %end;
        %end;
 
        %if &printlogstats = 1 %then %put  sample_s = &sample_s , sample_r = &sample_r  , seed = &seed  ;

        %let seed = %eval(&seed+3);


       proc datasets library = work nolist ;
        delete _betar_ _seedr_   %if &sample_r > 0 and &sample_s > 0 %then _simulsample_ ; param simul  _covbounds_ 
               %if &runnc = 1 %then surv_tmp0 ;
               %do i = 1 %to &numint ; surv_tmp&i %end;
               ;
        quit;
   
    %end;
 %end;

    
 %end;


     %if &runnc = 1 OR (&runnc = 0 AND &numint > 0 AND &refint > 0 ) %then %do;
        %if &chunked = 0 %then %do;
           %*Summarizing final results;
           %results_blb(blb_samples = &BLB_s );
        %end;
        %else %if &chunked = 1 AND &sample_start = 0 %then %do;
           %local visitlist ;
           %let visitlist = ;
           %do i = 1 %to &ncov;
              %if &&usevisitp&i = 1 %then %let visitlist = &visitlist &i ;
           %end;
           %bootstrap_results(
               bootlib = &savelib,
               bootname = &intervname ,
               outc = &outc ,
               outctype = &outctype,
               compevent = &compevent ,
               censor = &censor ,
               check_cov_models = &check_cov_models,
               covmeandata = &covmeanname0 ,
               ncov = &ncov,
               usevisitp = &visitlist ,
               survdata=&survdata,
               observed_surv = &observed_surv,
               print_cov_means = &print_cov_means,
               savecovmean = 0,
               bootstrap_hazard = 0,
               hazardratio=&hazardratio,
               hazardname = &hazardname0 ,
               intcomp = &intcomp ,
               time = &time ,
               timepoints = &timepoints,
               numparts = 1,
               samplestart = 0 ,
               sampleend = &sample_end,
               numboot = &sample_end,
               numint = &numint ,
               refint = &refint ,
               rungraphs = &rungraphs ,
               graphfile = &graphfile,
               resultsdata = &resultsdata,
               runnc = &runnc ,              
               titledata= &titledata,
               title1= &title1,
               title2= &title2,
               title3= &title3,
               tsize = &tsize,
               printlogstats = &printlogstats
               );
              
           
       %end;
    %end;
%mend ;


