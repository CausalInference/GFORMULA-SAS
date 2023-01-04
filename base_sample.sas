%macro base_sample ;
/*******
   want code for original data , no looping 
   for bsample = 0 
************/

    %if &sample_start = 0 %then %let bsample = 0 ;

    %*Looping over bootstrap samples;
   
        %if (&outputs ^= yes or %eval(&bsample) ^= 0) %then %do;
               %let ods_logit = ods select none ;
               %let ods_reg = ods select none ;
               %let ods_qlim = ods select none ;
        %end;
        %*Generating samples;

		%if &printlogstats = 1 %then %put  before sample = &bsample seed = &seed ;
		%if &printlogstats = 1 %then %put  ;

	 

        %*Base data for estimation and simulation (sample 0);
        
        data param;
        set _paramdata_;
        _sample_ = 0;
		%if &bootstrap_method = 0 %then %do;
			bootstrap_counts = 1 ;
		%end;
		%else %if &bootstrap_method = 1 %then %do;
			BLB_counts = 1 ;
		%end;
        run;
        
       %let subseed = %eval(2*&seed);
       * _simuldata_ contains baseline variables ;
       %if &nsimul = &ssize %then %do;
            data simul;
            set _simuldata_;
            _sample_ = 0;
			%if &bootstrap_method = 0 %then %do;
				bootstrap_counts = 1 ;
			%end;
			%else %if &bootstrap_method = 1 %then %do;
				BLB_counts = 1 ;
			%end;
            run;
        %end;
        %else %do;
       
         

            proc surveyselect data=_simuldata_ out=simul noprint
                method=urs  sampsize=&nsimul seed=&subseed;
            run;

            data simul ;
            set simul ;
             _sample_ = 0 ;
			 %if &use_bootstrap_counts = 0 %then %do;
            	do _copy_ = 1 to numberhits ;
                	output ;
            	end;
           
				%if &bootstrap_method = 0 %then %do;
					bootstrap_counts = 1 ;
				%end;
				%else %if &bootstrap_method = 1 %then %do;
					BLM_counts = 1 ;
				%end;
				drop _copy_ numberhits ;
		     %end;
			 %else %do;
				%if &bootstrap_method = 0 %then %do;
					bootstrap_counts = numberhits ;
				%end;
				%else %if &bootstrap_method = 1 %then %do;
					BLB_counts = numberhits ;
				%end;
				drop numberhits ;
			 %end;            
            run;
        %end;             
                  
   

          
        %*Estimating parameters;
        %if &usebetadata = 0 %then %do;            
                %if &uselabelc = 0  and &uselabelo = 0 %then options nolabel ;;
                %parameters;     
                options &label_orig ;      
        %end;

        %else %do; /*feb2013*/

            data _beta_;
            set &betadata;
            run;

        %end;   

 
        ods select all ; 

        data _betar_ ;
        set _beta_  (where = ( _sample_=&bsample)) ;
        run;

        data _seedr_ ;
        _seedr_ = %eval(&seed);
        _sample_ = &bsample ;
        run;

       %if &bsample = &sample_start  %then %do;
             /* initialize data views for interventions */

             %if &runnc = 1 %then %do;
                 %interv_init(intno=0, intlabel='Natural course' ); 
             %end;   
        

             %do intnum = 1 %to &numint;
                %interv_init(&&interv&intnum);
             %end;             

        %end;            
       
        %if &hazardratio = 1 %then %do;
             %if &bsample = 0 %then  %createhazard ;
             %else %if &bootstrap_hazard = 1 %then %createhazard ;
        %end;

        %if &runnc = 1 %then %do;  /*** for natural course results ***/
            %*No intervention cumulative incidence;

            %if   %bquote(&simuldata)^=  %then %do;

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
                    freq &bootstrap_counts ;
                    title "mean, min, max under intervention 0";
                run;


                data interv0;
                set interv0;
                _sample_ = 0;
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
				   freq &bootstrap_counts ;
                   output out = survprobs0 mean = %do n = 1 %to &timepoints; risk&n surv&n compevent&n  %end; ;
                   run;

                    /* initialize survdata */

                    data surv_tmp0 ;
                        set survprobs0 ;
                        length int2 $70 ;
                        int = 0 ;
                        int2 = "Natural course";
                        _sample_ = 0 ;
                        surv0 = 1;
                        n = _freq_ ;
                        keep  int int2 _sample_ n  surv0
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
                        _sample_ = 0 ;                                        
                        keep  int int2 _sample_ s&outc  ;
                    run;
                %end; 
            %end;  /* end of saving simuldata , needed outside of interv and interv_init macros */
            %else %do;                   
                %interv (intno=0 , intlabel='Natural course'); /* do not need to save simulated data set */                                 
            %end;   

         
                %if &bsample = &sample_start %then %do;  
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
                %if &intnum = 1 AND &runnc = 0 AND &bsample = &sample_start %then %do;
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

        %if %eval(&bsample) = &sample_start %then %do;  
           
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
             %if &bsample = &sample_start %then %do; 

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
 
        %if &printlogstats = 1 %then %put  sample = &bsample , seed = &seed  ;

        %let seed = %eval(&seed+3);


       proc datasets library = work nolist ;
        delete _betar_ _seedr_   %if &bsample > 0 %then _simulsample_ ; param simul  _covbounds_ 
               %if &runnc = 1 %then surv_tmp0 ;
               %do i = 1 %to &numint ; surv_tmp&i %end;
               ;
        quit;
   
   


%mend base_sample  ;
