/* Adapted from example2.sas (CausalInference/GFORMULA-SAS):
   the parameterized %create_sample(event=) macro and its PROC MEANS,
   invoked with event=dia. The macro uses %if/%upcase/%else to generate
   different DO-UNTIL stop conditions and post-processing per outcome
   type, with 12-period covariate arrays and round()/rand('normal').
   The %include and trailing %gformula() call are omitted. Deterministic
   via seed 5027. */

%macro create_sample(event = dia) ;
%let condition = ;
%if %upcase(&event) = DIA %then %let condition = dia or censlost or dead ;
%else %if %upcase(&event) = CONT_E %then %let condition = censlost   ;
%else %if %upcase(&event) = BIN_E %then %let condition = dead or  censlost   ;

**SAMPLE Data;
    data sample(drop = i j ahbp1-ahbp12 abmi1-abmi12 );

        call streaminit(5027);

        do i=1 to 1000;
            baseage = int( 35 + 25*rand('uniform'));

            array ahbp(12);
            array abmi(12);

            do j=1 to 12;
                ahbp(j) = (0.4>rand('uniform'));
                if j > 1 & ahbp(j-1) = 1 then ahbp(j) = 1 ;
                abmi(j) = round((25+5*(rand('normal'))),0.001);
            end;

            do j=3 to 12  until ( &condition   ) ;
                id=i;
                time=j-3;

                hbp     = ahbp(j);
                hbp_l1  = ahbp(j-1);
                hbp_l2  = ahbp(j-2);
                hbp_b   = ahbp(3);

                bmi     = abmi(j);
                bmi_l1  = abmi(j-1);
                bmi_l2  = abmi(j-2);
                bmi_b   = abmi(3);

                dia = ( (j/500) >rand('uniform'));

                if time < 9 then censlost  = (0.05>rand('uniform'));
                else censlost = 0 ;
                %if %upcase(&event) = DIA  or %upcase(&event) = BIN_E %then dead  = (0.05>rand('uniform'));;
                if time = 9 then dead = . ;
                if time = 9 then cont_e =  round((bmi+5*(rand('normal'))),0.01)  ;
                else bmi_e = . ;
                if time = 9 then bin_e = rand('bernoulli',0.6);
                else bin_e = . ;

                output;

            end;
    end;

    run;

    data sample ;
    set sample ;
    %if %upcase(&event)=DIA %then %do;
      if censlost = 1 then do ;
          dead = . ;
      end;
      if censlost = 1 or dead = 1 then do ;
           dia = . ;
      end;
    %end;
    %if %upcase(&event) = CONT_E %then %do;
       if time < 9 then cont_e = . ;
       if censlost = 1  then do ;
           cont_e = . ;
      end;
    %end;
    %if %upcase(&event)=BIN_E %then %do;
        if censlost = 1 then do ;
            dead = . ;
        end ;
        if time < 9 then bin_e = . ;
        if censlost = 1 or dead = 1 then do ;
           bin_e = . ;
       end;
    %end;

    run;

    proc means data=sample;
    title 'Means of SAMPLE data';
    run;
%mend ;

%create_sample(event = dia) ;
