%let n = 20 ;
%let b = 6;
%let s = 5 ;
%let r = 3 ;

data step0 ;
call streaminit(123);
do id = 1 to &n ;
    x = rand('uniform');
	output ;
end;
run;

proc surveyselect data = step0(keep =id ) out = step1 (rename = (replicate = sample_si))
   method = srs 
   sampsize = &b 
   reps = &s 
   noprint 
   ;
run ;

proc surveyselect data = step1 
  out=step2 (drop = ExpectedHits SamplingWeight rename = (replicate = sample_rj numberhits = BLB_weight))
  method = urs 
  sampsize = &n 
  reps = &r 
  noprint 
 ;
  strata sample_si ;
run;

%macro silly ;

	%do i = 1 %to &s ;
		%do j = 1 %to &r ;
			data blb_sample ;
			merge step0 step2(where = (sample_si = &i and sample_rj = &j ));
			by id ;	
			if blb_weight ne . ; 
			run;

			*title "sample_si = &i and sample_rj = &j";
			proc means data = blb_sample noprint;
			var x ;
			freq blb_weight ;
			output out = mean_si_rj (keep = sample_si sample_rj mean_x) mean = mean_x ;
			run;

			%if &j = 1 %then %do;
				data test_r ;
				set mean_si_rj ;
				run;
			%end;
			%else %do ;
				data test_r ;
				set test_r mean_si_rj ;
				run;
            %end;
		%end;
		   title " for i = &i" ;
			proc means data = test_r ;
			var mean_x ;
			output out = var_meanx (keep = var_meanx ) var = var_meanx ;
			run;

			%if &i = 1 %then %do;
				data test_s ;
				set var_meanx ;
				run;
			%end;
			%else %do;
				data test_s ;
				set test_s var_meanx ;
				run;
			 %end;

	%end;




%mend ;


options mprint notes ;

%silly ;
