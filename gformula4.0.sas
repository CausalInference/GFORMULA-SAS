/*
  
GFORMULA SAS MACRO

Authors: Roger W. Logan, Jessica  G. Young, Sarah L. Taubman, Yu-Han Chiu,  Emma McGee , Sally Picciotto, Goodarz Danaei, Miguel A. Hernán

*** EM edited ***
Version September 2024. This version includes additions and fixes for (i) the parametric outcome models available for end of follow-up 
outcomes, (ii) the functional forms of covariate histories allowed in these models, (iii) restricted cubic spline variables, and 
(iv) options for carrying forward covariates during skipped times (i.e., when no measurement of a covariate occurs).

Version April 2022. This version includes additions and fixes for the inclusion of censoring in the calculation of the natural 
course risks and means of covariates under the simulation of the natural course.

Version September 2021. This viersion includes fixes to rcspline macro for using negative knots and listpred for 
including variabls of ptype cumavg twice in the model lists.

Version January 2019. This version includes options and improvements that are not compatible with previous versions 
of the software. For questions and comments, email rwlogan@hsph.harvard.edu or jyoung@hsph.harvard.edu
 
Copyright (c) 2007, 2021, The President and Fellows of Harvard College

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated 
documentation files (the "Software"), to deal in the Software without restriction, including without limitation 

the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.

This software is provided under the standard MIT License:

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED 
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
DEALINGS IN THE SOFTWARE.

*/ 

options mautosource minoperator ;


%macro gformula(
    data =,         /* data set for analysis */
    id =,           /* unique identifier for subjects */
    time =,         /* time point variable */
    timeptype =,
    timeknots =,
    timeinc = 1,
    timefuncgen =,  /* user-defined macro for creating general functions of time beyond the general ptype used in all models. To use these
                        functions of time the user will need to also set covXaddvars for including the variables in the various covariate models
                        and for covXsaddvars for providing predictors in the simulated array for each variable the new funtion is used in. 
                        (not currently supported, additional description to be added). */
    timepoints =,    /* number of time points */
    interval = 1,    /* length of time between time points (assumed constant, default=1 time period) */

    outc =,                   /* outcome variable: if dichotomous, then 1=event, 0=no event */
    outctype = binsurv,        /* outcome type: default=binsurv (binary survival time analysis), 
                               other options: bineofu (binary end of follow-up),
                                              cateofu (categorical or ordinal end of follow-up), *** EM edited ***
                                              conteofu (continuous end of follow-up modeled using linear regression), 
                                              conteofu2 (continuous end of follow-up modeled using truncated normal regression), *** EM edited ***
                                              conteofu3 (continuous end of follow-up modeled using Tobit regression), *** EM edited *** 
                                              conteofu4 (continuous end of follow-up modeled using logistic to log-linear approach), *** EM edited *** */
    outcinteract =,           /* interaction terms (between dichotomous covariates) to include in model for outc */
    outcwherem = (1=1),     /* outcome: (optional) condition under which to model outcome */
    outcwherenosim =(1=0),  /* outcome: (optional) condition under which not to simulate the outcome under intervention but assign some fixed value at a 
                               given time*/
    outcnosimelsemacro =,   /* user defined macro to evaluate when  outcwherenosim holds*/
    usehistory_eof = 0 ,   /* for eof type outcomes use complete history of all modeled covariates */

    compevent =,             /* competing risk event */   
    compeventinteract =,     /* interaction terms (between dichotomous covariates) to include in model for censoring by death from other causes */
    compeventwherem = (1=1), /* compevent: (optional) condition under which to model outcome */
    compeventwherenosim =(1=0), /* compevent: (optional) condition under which not to simulate the compevent under intervention but assign some fixed value at a 
                                given time*/
    compeventnosimelsemacro =, /* user defined macro to evaluate simulated compevent when  compeventeadwherenosim holds*/

    censor = , 
	maxipw = p99 , /* pX or number for maximum value, calculated using proc means (needs to be valid value ) */
	censorinteract= ,
    censorwherem = (1=1), /* compevent: (optional) condition under which to model outcome */
   
	compevent_cens = 0 , 
    
    fixedcov =,                /* predictors that are not predicted themselves */
    ncov =,                    /* number of (time-varying) covariates to be parametrically estimated */

    cov1 =,             /* covariate 1: name */
    cov1otype = 1,      /* covariate 1: outcome type (default = binary) */
    cov1ptype =,        /* covariate 1: predictor type */
    cov1etype= ,       /* covariate 1 predictor type for eof type outcome models. Not used for covariate models. Two parts: first is the general
	                      predictor type and second is how many lags to include in model (all none or k) when second is k this will include the lags from 1
	                      to k. */
    cov1cumint = ,      /* covariate to include as factor in calculation of average. should be {0,1}-variable */
    cov1mtype = all,    /* defines how variable is included in models, all for all models and some for including variable in some of the models */ 
    cov1knots =,        /* covariate 1 (-cat, -spl types only): list of knots for categories or splines */
	cov1eknots=,
    cov1skip = -1,      /* covariate 1: time points in which variable is not to be predicted (unmeasured) */
    cov1inc = 0,        /* covariate 1 (type 0 only): fixed increment for linearly increasing variables */ 
    cov1interact =,     /* interaction terms (between dichotomous covariates) to be included in the model for cov1 */
    cov1wherem = (1=1), /* covariate 1: (optional) condition under which to model the covariate 
                           (in terms of measured covariates) */
   cov1wherenosim=(1=0),  /* covariate 1: (optional) condition under which not to simulate a covariate under intervention but assign some fixed value at a 
                            given time*/
    cov1nosimelsemacro =, /* user defined macro to evaluate when  cov1wherenosim holds*/
    cov1class =,          /* covariate 1 (type 4 only): seperate models by lagged cov1 */
    cov1classelse =,      /* covariate 1 (type 7 only): default value to assign to cov1class in simulation if not simulated - change JGY*/ 

    /* user-defined macros that may be used for overriding automated features of main macro for cov1 */
    cov1addvars = ,         /* list of additional variables to use as predictors of cov1 new change JGY*/  
    cov1genmacro = ,        /* user defined macro to generate a variable in simulated dataset based on simulated cov1 or below new change JGY*/  
    cov1modusermacro =,     /* user defined macro implementing model fitting for otype=-1 */
    cov1moddatausermacro =, /* user defined macro that creates data sets with estimated coefficients from code in cov1modusermacro, otype=-1 only*/
    cov1setblvar =,         /* user defined macro to set value of a baseline variable not previously set in macro, otype=-1 only */
    cov1simusermacro = ,    /* user defined macro for simulation from estimated densiites models fit in cov1modusermacro, otype=-1 only */
    cov1barray = ,          /* user defined macro with any extra necessary array statements containing coefficients, otype = -1 only */ 
    cov1sarray = ,          /* user defined macro with any extra necessary array statements containing simulated variable lists, otype = -1 only*/
    cov1randomvisitp = ,    /* random visit process to model for given variable : different from skip which is a fixed visit process */
    cov1visitpmaxgap =9e10, /* maximum gap allowed between visits in visit process  */     
    cov1visitpwherem =(1=1),/* where conditions to use in the visit process */
    cov1visitpcount = ,     /* variable to use for initiating the counter for time since last visit in visit process, taken to be lag1 time since last visit at baseline */  
     
    cov2=,cov2otype=1,cov2ptype=,cov2etype=,cov2mtype= all ,cov2cumint= ,cov2skip=-1,cov2inc=0,cov2knots=,cov2eknots=,cov2interact=,cov2wherem=(1=1),cov2wherenosim=(1=0),
    cov2nosimelsemacro=, cov2class=, cov2classelse=, cov2addvars=,cov2genmacro=,cov2modusermacro =,cov2moddatausermacro =,cov2setblvar =,
    cov2simusermacro = ,cov2barray = ,cov2sarray =,cov2randomvisitp=,cov2visitpmaxgap=9e10,cov2visitpwherem=(1=1),cov2visitpcount = ,

    cov3=,cov3otype=1,cov3ptype=,cov3etype=,cov3mtype= all,cov3cumint= ,cov3skip=-1,cov3inc=0,cov3knots=,cov3eknots=,cov3interact=,cov3wherem=(1=1),cov3wherenosim=(1=0),
    cov3nosimelsemacro=,  cov3class=, cov3classelse=,cov3addvars=,cov3genmacro=,cov3modusermacro =,cov3moddatausermacro =,cov3setblvar =,
    cov3simusermacro = ,cov3barray = ,cov3sarray =,cov3randomvisitp=,cov3visitpmaxgap=9e10,cov3visitpwherem=(1=1),cov3visitpcount = ,
    cov3classstart =,

    cov4=,cov4otype=1,cov4ptype=,cov4etype=,cov4mtype= all ,cov4cumint= ,cov4skip=-1,cov4inc=0,cov4knots=,cov4eknots=,cov4interact=,cov4wherem=(1=1),cov4wherenosim=(1=0),
    cov4nosimelsemacro=,  cov4class=, cov4classelse=,cov4addvars=,cov4genmacro=,cov4modusermacro =,cov4moddatausermacro =,cov4setblvar =,
    cov4simusermacro = ,cov4barray = ,cov4sarray = ,cov4randomvisitp=,cov4visitpmaxgap=9e10,cov4visitpwherem=(1=1),cov4visitpcount = ,

    cov5=,cov5otype=1,cov5ptype=,cov5etype=,cov5mtype= all ,cov5cumint= ,cov5skip=-1,cov5inc=0,cov5knots=,cov5eknots=,cov5interact=,cov5wherem=(1=1),cov5wherenosim=(1=0),
    cov5nosimelsemacro=,  cov5class=, cov5classelse=,cov5addvars=,cov5genmacro=,cov5modusermacro =,cov5moddatausermacro =,cov5setblvar =,
    cov5simusermacro = ,cov5barray = ,cov5sarray =,cov5randomvisitp=,cov5visitpmaxgap=9e10,cov5visitpwherem=(1=1),cov5visitpcount = ,

    cov6=,cov6otype=1,cov6ptype=,cov6etype=,cov6mtype= all ,cov6cumint= ,cov6skip=-1,cov6inc=0,cov6knots=,cov6eknots=,cov6interact=,cov6wherem=(1=1),cov6wherenosim=(1=0),
    cov6nosimelsemacro=,  cov6class=, cov6classelse=,cov6addvars=,cov6genmacro=,cov6modusermacro =,cov6moddatausermacro =,cov6setblvar =,
    cov6simusermacro = ,cov6barray = ,cov6sarray =,cov6randomvisitp=,cov6visitpmaxgap=9e10,cov6visitpwherem=(1=1),cov6visitpcount = ,

    cov7=,cov7otype=1,cov7ptype=,cov7etype=,cov7mtype= all ,cov7cumint= ,cov7skip=-1,cov7inc=0,cov7knots=,cov7eknots=,cov7interact=,cov7wherem=(1=1),cov7wherenosim=(1=0),
    cov7nosimelsemacro=,  cov7class=, cov7classelse=,cov7addvars=,cov7genmacro=,cov7modusermacro =,cov7moddatausermacro =,cov7setblvar =,
    cov7simusermacro = ,cov7barray = ,cov7sarray =,cov7randomvisitp=,cov7visitpmaxgap=9e10,cov7visitpwherem=(1=1),cov7visitpcount = ,

    cov8=,cov8otype=1,cov8ptype=,cov8etype=,cov8mtype= all ,cov8cumint= ,cov8skip=-1,cov8inc=0,cov8knots=,cov8eknots=,cov8interact=,cov8wherem=(1=1),cov8wherenosim=(1=0),
    cov8nosimelsemacro=,  cov8class=, cov8classelse=,cov8addvars=,cov8genmacro=,cov8modusermacro =,cov8moddatausermacro =,cov8setblvar =,
    cov8simusermacro = ,cov8barray = ,cov8sarray =,cov8randomvisitp=,cov8visitpmaxgap=9e10,cov8visitpwherem=(1=1),cov8visitpcount = ,

    cov9=,cov9otype=1,cov9ptype=,cov9etype=,cov9mtype= all ,cov9cumint= ,cov9skip=-1,cov9inc=0,cov9knots=,cov9eknots=,cov9interact=,cov9wherem=(1=1),cov9wherenosim=(1=0),
    cov9nosimelsemacro=,  cov9class=, cov9classelse=,cov9addvars=,cov9genmacro=,cov9modusermacro =,cov9moddatausermacro =,cov9setblvar =,
    cov9simusermacro = ,cov9barray = ,cov9sarray =,cov9randomvisitp=,cov9visitpmaxgap=9e10,cov9visitpwherem=(1=1),cov9visitpcount = ,

    cov10=,cov10otype=1,cov10ptype=,cov10etype=,cov10mtype= all ,cov10cumint= ,cov10skip=-1,cov10inc=0,cov10knots=,cov10eknots=,cov10interact=,cov10wherem=(1=1), 
    cov10wherenosim=(1=0),cov10nosimelsemacro=,  cov10class=, cov10classelse=,cov10addvars=,cov10genmacro=,cov10modusermacro =,
    cov10moddatausermacro =,cov10setblvar =,cov10simusermacro = ,cov10barray = ,cov10sarray =,cov10randomvisitp=,cov10visitpmaxgap=9e10,cov10visitpwherem=(1=1),
    cov10visitpcount= ,

    cov11=,cov11otype=1,cov11ptype=,cov11etype=,cov11mtype= all ,cov11cumint= ,cov11skip=-1,cov11inc=0,cov11knots=,cov11eknots=,cov11interact=,cov11wherem=(1=1), 
    cov11wherenosim=(1=0),cov11nosimelsemacro=,  cov11class=, cov11classelse=,cov11addvars=,cov11genmacro=,cov11modusermacro =,
    cov11moddatausermacro =,cov11setblvar =,cov11simusermacro = ,cov11barray = ,cov11sarray =,cov11randomvisitp=,cov11visitpmaxgap=9e10,cov11visitpwherem=(1=1),
    cov11visitpcount= ,

    cov12=,cov12otype=1,cov12ptype=,cov12etype=,cov12mtype= all ,cov12cumint= ,cov12skip=-1,cov12inc=0,cov12knots=,cov12eknots=,cov12interact=,cov12wherem=(1=1),
    cov12wherenosim=(1=0),cov12nosimelsemacro=,  cov12class=, cov12classelse=,cov12addvars=,cov12genmacro=,cov12modusermacro =,
    cov12moddatausermacro =,cov12setblvar =,cov12simusermacro = ,cov12barray = ,cov12sarray =,cov12randomvisitp=,cov12visitpmaxgap=9e10,cov12visitpwherem=(1=1),
    cov12visitpcount= ,

    cov13=,cov13otype=1,cov13ptype=,cov13etype=,cov13mtype= all ,cov13cumint= ,cov13skip=-1,cov13inc=0,cov13knots=,cov13eknots=,cov13interact=,cov13wherem=(1=1), 
    cov13wherenosim=(1=0),cov13nosimelsemacro=,  cov13class=, cov13classelse=,cov13addvars= ,cov13genmacro=,cov13modusermacro =,
    cov13moddatausermacro =,cov13setblvar =,cov13simusermacro = ,cov13barray = ,cov13sarray =,cov13randomvisitp=,cov13visitpmaxgap=9e10,cov13visitpwherem=(1=1),
    cov13visitpcount= ,

    cov14=,cov14otype=1,cov14ptype=,cov14etype=,cov14mtype= all ,cov14cumint= ,cov14skip=-1,cov14inc=0,cov14knots=,cov14eknots=,cov14interact=,cov14wherem=(1=1), 
    cov14wherenosim=(1=0),cov14nosimelsemacro=,  cov14class=, cov14classelse=,cov14addvars= ,cov14genmacro=,cov14modusermacro =,
    cov14moddatausermacro =,cov14setblvar =,cov14simusermacro = ,cov14barray = ,cov14sarray =, cov14randomvisitp=,cov14visitpmaxgap=9e10,cov14visitpwherem=(1=1),
    cov14visitpcount= ,

    cov15=,cov15otype=1,cov15ptype=,cov15etype=,cov15mtype= all ,cov15cumint= ,cov15skip=-1,cov15inc=0,cov15knots=,cov15eknots=,cov15interact=,cov15wherem=(1=1), 
    cov15wherenosim=(1=0),cov15nosimelsemacro=,  cov15class=, cov15classelse=,cov15addvars= ,cov15genmacro=,cov15modusermacro =,
    cov15moddatausermacro =,cov15setblvar =,cov15simusermacro = ,cov15barray = ,cov15sarray =, cov15randomvisitp=,cov15visitpmaxgap=9e10,cov15visitpwherem=(1=1),
    cov15visitpcount= ,

    cov16=,cov16otype=1,cov16ptype=,cov16etype=,cov16mtype= all ,cov16cumint= ,cov16skip=-1,cov16inc=0,cov16knots=,cov16eknots=,cov16interact=,cov16wherem=(1=1), 
    cov16wherenosim=(1=0),cov16nosimelsemacro=,  cov16class=, cov16classelse=,cov16addvars= ,cov16genmacro=,cov16modusermacro =,
    cov16moddatausermacro =,cov16setblvar =,cov16simusermacro = ,cov16barray = ,cov16sarray =, cov16randomvisitp=,cov16visitpmaxgap=9e10,cov16visitpwherem=(1=1),
    cov16visitpcount= ,

    cov17=,cov17otype=1,cov17ptype=,cov17etype=,cov17mtype= all ,cov17cumint= ,cov17skip=-1,cov17inc=0,cov17knots=,cov17eknots=,cov17interact=,cov17wherem=(1=1), 
    cov17wherenosim=(1=0),cov17nosimelsemacro=,  cov17class=, cov17classelse=,cov17addvars= ,cov17genmacro=,cov17modusermacro =,
    cov17moddatausermacro =,cov17setblvar =,cov17simusermacro = ,cov17barray = ,cov17sarray =, cov17randomvisitp=,cov17visitpmaxgap=9e10,cov17visitpwherem=(1=1),
    cov17visitpcount= , 

    cov18=,cov18otype=1,cov18ptype=,cov18etype=,cov18mtype= all ,cov18cumint= ,cov18skip=-1,cov18inc=0,cov18knots=,cov18eknots=,cov18interact=,cov18wherem=(1=1), 
    cov18wherenosim=(1=0),cov18nosimelsemacro=,  cov18class=, cov18classelse=,cov18addvars= ,cov18genmacro=,cov18modusermacro =,
    cov18moddatausermacro =,cov18setblvar =,cov18simusermacro = ,cov18barray = ,cov18sarray =, cov18randomvisitp=,cov18visitpmaxgap=9e10,cov18visitpwherem=(1=1),
    cov18visitpcount= , 

    cov19=,cov19otype=1,cov19ptype=,cov19etype=,cov19mtype= all ,cov19cumint= ,cov19skip=-1,cov19inc=0,cov19knots=,cov19eknots=,cov19interact=,cov19wherem=(1=1), 
    cov19wherenosim=(1=0),cov19nosimelsemacro=  ,  cov19class=, cov19classelse=,cov19addvars= ,cov19genmacro=,cov19modusermacro =,
    cov19moddatausermacro =,cov19setblvar =,cov19simusermacro = ,cov19barray = ,cov19sarray =, cov19randomvisitp=,cov19visitpmaxgap=9e10,cov19visitpwherem=(1=1),
    cov19visitpcount= , 

    cov20=,cov20otype=1,cov20ptype=,cov20etype=,cov20mtype= all ,cov20cumint= ,cov20skip=-1,cov20inc=0,cov20knots=,cov20eknots=,cov20interact=,cov20wherem=(1=1), 
    cov20wherenosim=(1=0),cov20nosimelsemacro=,  cov20class=, cov20classelse=,cov20addvars= ,cov20genmacro=,cov20modusermacro =,
    cov20moddatausermacro =,cov20setblvar =,cov20simusermacro = ,cov20barray = ,cov20sarray =, cov20randomvisitp=,cov20visitpmaxgap=9e10,cov20visitpwherem=(1=1),
    cov20visitpcount= , 

    cov21=,cov21otype=1,cov21ptype=,cov21etype=,cov21mtype= all ,cov21cumint= ,cov21skip=-1,cov21inc=0,cov21knots=,cov21eknots=,cov21interact=,cov21wherem=(1=1), 
    cov21wherenosim=(1=0),cov21nosimelsemacro=,  cov21class=, cov21classelse=,cov21addvars=,cov21genmacro=,cov21modusermacro =,
    cov21moddatausermacro =,cov21setblvar =,cov21simusermacro = ,cov21barray = ,cov21sarray =, cov21randomvisitp=,cov21visitpmaxgap=9e10,cov21visitpwherem=(1=1),
    cov21visitpcount= , 

    cov22=,cov22otype=1,cov22ptype=,cov22etype=,cov22mtype= all ,cov22cumint= ,cov22skip=-1,cov22inc=0,cov22knots=,cov22eknots=,cov22interact=,cov22wherem=(1=1), 
    cov22wherenosim=(1=0),cov22nosimelsemacro=,  cov22class=, cov22classelse=,cov22addvars=,cov22genmacro=,cov22modusermacro =,
    cov22moddatausermacro =,cov22setblvar =,cov22simusermacro = ,cov22barray = ,cov22sarray =, cov22randomvisitp=,cov22visitpmaxgap=9e10,cov22visitpwherem=(1=1),
    cov22visitpcount= , 

       cov23=,cov23otype=1,cov23ptype=,cov23etype=,cov23mtype= all ,cov23cumint= ,cov23skip=-1,cov23inc=0,cov23knots=,cov23eknots=,cov23interact=,cov23wherem=(1=1), 
    cov23wherenosim=(1=0),cov23nosimelsemacro=,  cov23class=, cov23classelse=,cov23addvars=,cov23genmacro=,cov23modusermacro =,
    cov23moddatausermacro =,cov23setblvar =,cov23simusermacro = ,cov23barray = ,cov23sarray =, cov23randomvisitp=,cov23visitpmaxgap=9e10,cov23visitpwherem=(1=1),
    cov23visitpcount= , 

    cov24=,cov24otype=1,cov24ptype=,cov24etype=,cov24mtype= all ,cov24cumint= ,cov24skip=-1,cov24inc=0,cov24knots=,cov24eknots=,cov24interact=,cov24wherem=(1=1), 
    cov24wherenosim=(1=0),cov24nosimelsemacro=,  cov24class=, cov24classelse=,cov24addvars=,cov24genmacro=,cov24modusermacro =,
     cov24moddatausermacro =,cov24setblvar =,cov24simusermacro = ,cov24barray = ,cov24sarray =, cov24randomvisitp=,cov24visitpmaxgap=9e10,cov24visitpwherem=(1=1),
    cov24visitpcount= , 

    cov25=,cov25otype=1,cov25ptype=,cov25etype=,cov25mtype= all ,cov25cumint= ,cov25skip=-1,cov25inc=0,cov25knots=,cov25eknots=,cov25interact=,cov25wherem=(1=1), 
    cov25wherenosim=(1=0),cov25nosimelsemacro=,  cov25class=, cov25classelse=,cov25addvars=,cov25genmacro=,cov25modusermacro =,
    cov25moddatausermacro =,cov25setblvar =,cov25simusermacro = ,cov25barray = ,cov25sarray =, cov25randomvisitp=,cov25visitpmaxgap=9e10,cov25visitpwherem=(1=1),
    cov25visitpcount= , 

    cov26=,cov26otype=1,cov26ptype=,cov26etype=,cov26mtype= all ,cov26cumint= ,cov26skip=-1,cov26inc=0,cov26knots=,cov26eknots=,cov26interact=,cov26wherem=(1=1), 
    cov26wherenosim=(1=0),cov26nosimelsemacro=,  cov26class=, cov26classelse=,cov26addvars=,cov26genmacro=,cov26modusermacro =,
    cov26moddatausermacro =,cov26setblvar =,cov26simusermacro = ,cov26barray = ,cov26sarray =, cov26randomvisitp=,cov26visitpmaxgap=9e10,cov26visitpwherem=(1=1),
    cov26visitpcount= , 

    cov27=,cov27otype=1,cov27ptype=,cov27etype=,cov27mtype= all ,cov27cumint= ,cov27skip=-1,cov27inc=0,cov27knots=,cov27eknots=,cov27interact=,cov27wherem=(1=1), 
    cov27wherenosim=(1=0),cov27nosimelsemacro=,  cov27class=, cov27classelse=,cov27addvars=,cov27genmacro=,cov27modusermacro =,
    cov27moddatausermacro =,cov27setblvar =,cov27simusermacro = ,cov27barray = ,cov27sarray =, cov27randomvisitp=,cov27visitpmaxgap=9e10,cov27visitpwherem=(1=1),
    cov27visitpcount= , 

    cov28=,cov28otype=1,cov28ptype=,cov28etype=,cov28mtype= all ,cov28cumint= ,cov28skip=-1,cov28inc=0,cov28knots=,cov28eknots=,cov28interact=,cov28wherem=(1=1), 
    cov28wherenosim=(1=0),cov28nosimelsemacro=,  cov28class=, cov28classelse=,cov28addvars=,cov28genmacro=,cov28modusermacro =,
    cov28moddatausermacro =,cov28setblvar =,cov28simusermacro = ,cov28barray = ,cov28sarray =, cov28randomvisitp=,cov28visitpmaxgap=9e10,cov28visitpwherem=(1=1),
    cov28visitpcount= , 

    cov29=,cov29otype=1,cov29ptype=,cov29etype=,cov29mtype= all ,cov29cumint= ,cov29skip=-1,cov29inc=0,cov29knots=,cov29eknots=,cov29interact=,cov29wherem=(1=1), 
    cov29wherenosim=(1=0),cov29nosimelsemacro=,  cov29class=, cov29classelse=,cov29addvars=,cov29genmacro=,cov29modusermacro =,
    cov29moddatausermacro =,cov29setblvar =,cov29simusermacro = ,cov29barray = ,cov29sarray =, cov29randomvisitp=,cov29visitpmaxgap=9e10,cov29visitpwherem=(1=1),
    cov29visitpcount= , 

    cov30=,cov30otype=1,cov30ptype=,cov30etype=,cov30mtype= all ,cov30cumint= ,cov30skip=-1,cov30inc=0,cov30knots=,cov30eknots=,cov30interact=,cov30wherem=(1=1), 
    cov30wherenosim=(1=0),cov30nosimelsemacro=,  cov30class=, cov30classelse=,    cov30addvars=,cov30genmacro=,cov30modusermacro =,
     cov30moddatausermacro =,cov30setblvar =,cov30simusermacro = ,cov30barray = ,cov30sarray =, cov30randomvisitp=,cov30visitpmaxgap=9e10,cov30visitpwherem=(1=1),
    cov30visitpcount= ,
 
    /*other override options*/

    wherevars =,            /* list of variables referenced in any of the cov#wherem conditions, change JGY*/ 
    keepsimuldata =,        /*list of variables not created by a given ptype or otype that will be needed in simulated data set new change JGY*/
    equalitiessimuldata =,    /*user defined macro that equates pre-baseline simulated vars to observed new change JGY*/
    eventaddvars =,         /*list of variables to be added to event predictor list new change JGY*/
    compeventaddvars =,
    censoraddvars= ,
    usebetadata = 0,
    betadata =,       /* data set to store parameter estimates */
    simuldata =,      /* data set to store simulated data */
    resultsdata =,    /* data set to store results */
    survdata =,       /*data set to store cumulative survival probabilities at each time point under all interventions*/
    outputs = yes,    /* whether to print regression results */
    print_stats = 1,
    check_cov_models = 0, /* create data set for difference of mean of observed covs and mean of simulated covs under natural 
                             course */
    print_cov_means = 0,  /* print out tables of comparison of observed and simulated variables */
    covmeandata = ,
    save_raw_covmean = 0,
    observed_surv =,
    intervname =,
    
    refint = 0,
    seed = 7834,          /* random numbers seed */
    nsamples = 50,    /* number of bootstrap samples (default is 50) */
    nsimul =,         /* size (# subjects) of simulated sample, default is sample size */
    nparam =,         /* size (# subjects) of parameter sample, default is sample size */
    hazardratio =0,   /* calculate the hazard ratio for two interventions. This will increase the run time for the macro*/
    intcomp = ,       /* needs to be a list with two numbers for comparing the hazard ratio when hazardratio = 1 */
    bootstrap_hazard = 0, /* when running bootstrap samples also include calculation of hazard ratios */
    hazardname =,     /* name of data set to hold hazard ratio when runnig bootstraps in parts, data will be saved in savelib. */
    numint = 0,       /* number of interventions */
	sim_trunc = 1 ,   /* when simulating variables of otype = 3 (continuous regression) truncate values to be within observed bound */
    sample_start = 0, /* first sample to use in bootstraps (can be 0 for original data ) */
    sample_end = -1,  /* last sample to use in bootstraps (should be at most equal to nsamples) */
    savelib = work,   /* location for saving intermediate results for chunks */

    rungraphs = 0,
    title1 =,
    title2 =,
    title3 =,
    titledata = ,
    graphfile = gfile.pdf,
    tsize = 1,
    runnc = 1,
    weight =,
    printlogstats = 1,   
    checkaddvars = 1,
    minimalistic = no, /* only keep results for outcome variables */
    testing = no       /* keep each simulated data set for each intervention. will be named savelib.simulated&intno */
    );


    %* Setting some local parameters ;
    %local i j n intnum intstart  int bsample  uselabelo uselabelc ods_logit ods_reg ods_qlib    chunked seed_list seed_holder   ;
    %local  ssize obsp outcpred compeventpred censorpred censlpred dimoutc dimcompevent dimcensl    
       outcmin outcmax outcninterx compeventninterx censorninterx censlninterx covmeanname0  uselabelo uselabelc sample_hazard hazardname0 simulkeeplist;
    %local   cov0 cov0otype cov0ptype cov0mtype cov0skip cov0inc cov0knots cov0eknots cov0lev cov0interact cov0wherem cov0cumint
          cov0wherenosim cov0nosimelsemacro cov0class cov0classelse cov0addvars  cov0genmacro  cov0modusermacro 
          cov0moddatausermacro cov0setblvar cov0simusermacro cov0barray cov0sarray cov0randomvisitp cov0visitpmaxgap cov0visitpwherem
          cov0visitpcount cov0visitpelse cov0etype newsimulkeeplist;
    %local covlist ;
    %local anytsswitch1 anylagcumavg anyecumavg anyecumcat anyfindpknots anyfindeknots ; /* any ptype equal to tsswitch1 for creating spline variables */
    %let anytsswitch1 = 0 ;
    %let anylagcumavg = 0 ;
	%let anyecumavg = 0  ;
	%let anyfindpknots = 0 ;
	%let anyfindeknots = 0 ;
    %local created_global  ; /* keep track of global variables created in interaction macros, 
                               will remove at end of current run */
    
	%if %bquote(&compevent ) = %then %let compevent_cens = 0 ;

    %let covlist = ;
    %if %symexist(enforcegap) = 0 %then %let enforcegap = 0 ; 

	%if &minimalistic = 1 %then %let minimalistic = yes ;
	%else %if &minimalistic = 0 %then %let minimalistic = no ;
    %let minimalistic = %lowcase(&minimalistic);
    %let testing = %lowcase(&testing);
    %if &testing = yes %then %do;
         %let nsample = 0 ;
         %let minimalistic = no ;
    %end;
    %if nsamples = 0 %then %let bootstrap_hazard = 0 ;
    %if %upcase(&outctype) ^= BINSURV %then %do;
          %let hazardratio = ;
          %let bootstrap_hazard = 0;
          %let intcomp = ;
          %let hazardname = ;
		  %if %bquote(&censor) ^= AND %bquote(&compevent)^= %then %let compevent_cens = 1 ;
    %end; 
	%if %upcase(&outctype) = BINSURV %then %let usehistory_eof = 0 ;

    %if &minimalistic = yes %then %do;
        %let rungraphs = 0 ;
        %let print_cov_means = 0 ;
        %let check_cov_models = 0 ;
    %end;
    %*Producing error message if no dataset is specified;
    %if %bquote(&data)= %then %do;
        %put ERROR: The DATA= argument must be specified.;
        %goto exit;
    %end;
    %if &runnc = 0 and &refint = 0 %then %do;
         %put ERROR: When not running the natural course analysis, refint must be different from 0;
         %goto exit;
    %end;

    %let uselabelo = 0 ;
    %let uselabelc = 0 ;

    /* create a new covariate to hold time variables  will be cov0 . all loops for ncov should now start at 0 instead of 1 */
    %createcov0vars ;   
   
    %do i=0 %to %eval(&ncov);
       %local  cov&i.simstart cov&i.else cov&i.infixed cov&i.fixedcov cov&i.visitpelse cov&i.findknots cov&i.findeknots ;          
       
       %let cov&i.ptype = %sysfunc(left(&&cov&i.ptype)) ;
       %let cov&i.ptype = %lowcase(&&cov&i.ptype ) ;  
	   %let cov&i.etype = %lowcase(&&cov&i.etype) ;	  
       
       %let cov&i.simstart = 1 ;        
       %let tmp = &&cov&i ;        
       %let cov&i.else =  &tmp._b;  
       %let cov&i.infixed = 0 ;

 		/* create indicator of if the knots need to be found for the ptypes that are splines. Indicator variable is covXfindknots.
	      The knots will be found in the dataprep macro. For any categorical variables we will require the user to supply the knots. */
	   %let cov&i.findknots = 0 ;
	   %let cov&i.findeknots = 0 ;
	   %if ( %index(%upcase(&&cov&i.ptype),SPL) > 0 OR &&cov&i.ptype = tsswitch1  )  %then %do;
	   		%if %quote(%upcase(&&cov&i.knots)) = NONE  or  %numargs(&&cov&i.knots)   =  0  %then %let cov&i.knots = 0 ;

			
			%if %numargs(&&cov&i.knots) = 1 %then %do;
				%if &&cov&i.knots = 1 %then %let cov&i.knots = 3 ;
				%else %if &&cov&i.knots = 2 %then %let cov&i.knots = 4 ;
				%if &&cov&i.knots in (3 4 5) %then %do;
					%let cov&i.findknots = 1 ; 
					%let fineanypknots = 1 ;
				%end;
				%if %quote(%upcase(&&cov&i.knots)) =  0 AND  %index(%upcase(&&cov&i.ptype),SPL) > 0 %then %do  ;
				%put ERROR: WHEN USING COVARIATES OF TYPE LAGSPL OR SKPSPL, USER MUST DEFINE COV&I.KNOTS TO BE NON-EMPTY OR NON-ZERO ;
				%GOTO exit;
	        %end;
			%end;


	   %end;

      %if  %index(%upcase(&&cov&i.ptype),CAT) > 0   %then %do;
	   		%if %quote(%upcase(&&cov&i.knots)) = NONE  or %numargs(&&cov&i.knots)  = 0  %then %do  ;
				%put ERROR: WHEN USING COVARIATES OF TYPE -CAT, USER MUST DEFINE COV&I.KNOTS TO BE NON-EMPTY ;
				%GOTO exit;
	        %end;	 
	  %end;


       %if &&cov&i.otype=0 and &&cov&i.inc=  %then %do; 
           %put ERROR: Otype 0 increment must be specified;
           %goto exit;
       %end;

       %local usevisitp&i ;
       %let usevisitp&i = 0 ;
       %if %bquote(&&cov&i.randomvisitp)^= %then  %do;
             %let usevisitp&i = 1  ;   
             %let tmp = &&cov&i.randomvisitp ;
             %let cov&i.visitpelse = &tmp._b ;
      %end;



       %if not(%upcase(&&cov&i.mtype)  in (  ALL , SOME , NOCHECK ) )   %then %let cov&i.mtype = all ; 

      
	   %if &usehistory_eof = 1 AND &i > 0 %then %do;
			%local cov&i.etype_part1 cov&i.etype_part2 ;
			%let cov&i.etype = %sysfunc(compbl(&&cov&i.etype));		
			%let cov&i.etype_part1 = %scan(&&cov&i.etype,1, %str(' '));
			%let cov&i.etype_part2 = %scan(&&cov&i.etype,2,%str(' '));

			%if %upcase(&&cov&i.etype_part1) = NONE  OR %upcase(&&cov&i.etype_part2) = NONE OR &&cov&i.etype_part2 = 0 %then %do;
			       %let cov&i.etype_part1 = NONE ;
                   %let cov&i.etype_part2 = 0 ;
            %end;

			%let cov&i.etype = &&cov&i.etype_part1 ;

			%if &&cov&i.ptype = tsswitch1 AND (&&cov&i.etype in ( tsswitch1 cumsum cumsumnew )) %then %do;
				%if &&cov&i.etype_part1 = tsswitch1 %then %do;
				   %let cov&i.etype_part1 = cumsum ; 
                   %let cov&i.etype = cumsum ;
                   %put CHANGING cov&i.etype FROM TSSWITCH1 TO CUMSUM. ;
                %end;
			
			%end;

			%if %lowcase(&&cov&i.etype) in (cumsum cumsumnew cumavg cumavgnew ) /*and &usespline = 1 */ %then %do;
			
                %if (%quote(%upcase(&&cov&i.eknots)) = NONE ) OR ( %numargs(&&cov&i.eknots) = 0 )   %then %let cov&i.eknots = 0 ;			
                %if %numargs(&&cov&i.eknots) = 1 %then %do;
					%if &&cov&i.eknots = 1 %then %let cov&i.eknots = 3 ;
			    	%else %if &&cov&i.eknots = 2 %then %let cov&i.eknots = 4 ;
				    %if &&cov&i.eknots in (3 4 5) %then %let cov&i.findknots = 1 ; 
			    %end;

			%end;
               
			%if %index(&&cov&i.ptype,cat ) > 0  %then %do;
					%if %quote(%upcase(&&cov&i.knots)) = NONE  or %numargs(&&cov&i.knots) = 0 %then %do  ;
						%put ERROR: WHEN USING COVARIATES OF ETYPE -CAT, USER MUST DEFINE COV&I.EKNOTS TO BE NON-EMPTY ;
						%GOTO exit;
	        		%end;
			%end;

			%if &&cov&i.etype_part2 = %then %let cov&i.etype_part2 = %eval(&timepoints) ;
            %if %upcase(&&cov&i.etype_part2) = ALL %then %let cov&i.etype_part2 = %eval(&timepoints) ;
            
			%*put &i , &&cov&i , &&cov&i.etype ,&&cov&i.etype_part1 , &&cov&i.etype_part2 ;
	   %end;

       %if &usehistory_eof = 0 %then %do;
             %let cov&i.etype = ;
			 %let cov&i.eknots= 0;
	  %end;
 

       %if %upcase(&&cov&i.ptype) eq TSSWITCH1   %then %let anytsswitch1 = 1;
       %if %upcase(&&cov&i.ptype) = CUMAVG or %upcase(&&cov&i.ptype) eq LAG1CUMAVG or %upcase(&&cov&i.ptype) eq LAG2CUMAVG 
            %then %do;
                %let anylagcumavg = 1;
                %local cov&i._cumavg_l1_knots   ;
			
				%if %quote(%upcase(&&cov&i.knots)) = NONE or %numargs(&&cov&i.knots) = 0 %then %let cov&i.knots = 0 ;
				

                %if %numargs(&&cov&i.knots) > 1  %then %do;
					%let cov&i.lev=%numargs(&&cov&i.knots);
					%let cov&i._cumavg_l1_knots  = &&cov&i.knots ;
				%end;
			    %else %do ;
					%if &&cov&i.knots = 1 %then %let cov&i.knots = 3 ;
					%else %if &&cov&i.knots = 2 %then %let cov&i.knots = 4 ;
					%let cov&i.lev = &&cov&i.knots ;
					%let anyfindpknots = 1 ;
					%if &&cov&i.knots > 0 %then %let cov&i.findknots = 1 ;
			   %end;
			   %if &&cov&i.lev = 0 %then %let cov&i.lev = 2 ;
               %if &printlogstats = 1 %then %put  Number of levels of &&cov&i : cov&i.lev=&&cov&i.lev;                     
       %end;

	   %if &i > 0 AND &usehistory_eof = 1  %then %do;
	   		%if (%upcase(&&cov&i.etype) in CUMSUM CUMSUMNEW CUMAVG CUMAVGNEW )  %then %do; 
			    %if  %numargs(&&cov&i.eknots) = 1 %then %do;
                     %if &&cov&i.eknots > 0 %then %let cov&i.findknots = 1 ;
			    %end ;
	   			%let anyecumavg = 1 ;
			/*	%if &usespline = 0 %then %let covk&i.knots = 0 ; */
			%end;
			%if (%upcase(&&cov&i.etype) in ( CUMSUMCAT CUMSUMCATNEW CUMAVGCAT CUMAVGCATNEW ) )  %then %do; 
			    %if %quote(%upcase(&&cov&i.eknots)) = NONE or %numargs(&&cov&i.eknots)  = 0  %then %do;
					%put ERROR: WHEN USING AN ETYPE IN (CUMSUMCAT CUMSUMCATNEW CUMAVGCAT CUMAVGCATNEW ) REQUIRES THE SETTING OF COV&i.EKNOTS ) ;
					%goto exit ;
				%end;
	   			%let anyecumcat = 1 ;
				
			%end;
	   %end;


	/* %put &&cov&i : &&cov&i.ptype , &&cov&i.knots , &&cov&i.etype , &&cov&i.eknots , &&cov&i.findknots ; */
    %end;      
  
    %createfixedcov ;
    
    %if &printlogstats = 1 %then %put ;

    %let ods_logit =  ods select ModelInfo NObs ResponseProfile ConvergenceStatus 
                     ParameterEstimates  ;
    %let ods_reg = ods select NObs FitStatistics ParameterEstimates ;

    %let ods_qlim = ods select SummaryContResponse FitSummary ConvergenceStatus 
                     ParameterEstimates   ;

    %local label_orig ;
    %let label_orig = %sysfunc(getoption(nolabel));
    /*************************
	%if %bquote(&censor)^= %then %do;
		%let check_cov_models = 1 ; 
		%if %bquote(&covmeandata)= %then %do;
	            %if &printlogstats = 1 %then %put  SETTING COVMEAN = _COVMEAN_ ;
	            %let covmeandata = _covmean_ ;
	    %end;
	%end ;
 *****************************/

    %if &rungraphs = 1 %then %do;
         %if %bquote(&covmeandata)= %then %do;
            %if &printlogstats = 1 %then %put  SETTING COVMEAN = _COVMEAN_ ;
            %let covmeandata = _covmean_ ;
         %end;
         %if %bquote(&survdata)= %then %do;
             %if &printlogstats = 1 %then %put  SETTING SURVDATA=_SURVDATA_ ;
              %let survdata = _survdata_ ;
         %end;
         %if %bquote(&observed_surv)= %then %do;
            %if &printlogstats = 1 %then %put  SETTING OBSERVED_SURV = _OBSSURV_ ;
            %let observed_surv = _obssurv_ ; 
         %end ;

         %let check_cov_models = 1 ;
    %end;

    %if &ncov = 0 %then %do;
        %let check_cov_models = 0;
    %end;
    %if &sample_end =  -1 %then %let sample_end = &nsamples ;
    %if ( &sample_start = 0 AND  &sample_end = &nsamples ) %then %let chunked = 0 ;
    %else %let chunked = 1 ;
 
    /* rename saved data sets to contain first and last sample in the name */

    %if &chunked = 1 %then %do;
         /* rename the survival data if desired. need to consider the case where there
            is a libname. will change this to be the savelib variable */                           

         %if %bquote(&survdata)^= %then %do;
            %let _dot_ = %index(&survdata,%str(.)) ;
            %if &_dot_ > 0 %then   %let survname = %cmpres(%qscan(&survdata,2,%str(.))) ;
            %else %let survname = &survdata ; 
            %let survname = &survname._&sample_start._&sample_end ;
            %let survdata = &savelib..&survname ;
        %end;
        %else %do;
            %if &outctype=binsurv %then %do;
                %put TO RUN THE GFORMULA WITH  ONLY PART OF THE BOOTSTRAP ;
                %PUT SAMPLES YOU NEED TO SUPPLY SURVDATA TO BE NON-MISSING TO BE ABLE TO SAVE THE ;
                %PUT INTERMEDIATE RESULTS IN A CONSISTENT WAY. ;
                 %GOTO exit;
            %end;

            %else %let survdata = survdata_all ;
        %end;

       %if &check_cov_models = 1 %then %do;
             %if &save_raw_covmean = 1 %then %do;
                  %let save_raw_covmean = 0;
                  %if &printlogstats = 1 %then %put  WHEN RUNNING THE GFORMULA MACRO IN PARTS ALL THE RAW DATA FOR COVMEANDATA IS SAVE IN PARTS THAT WILL BE USED;
                  %if &printlogstats = 1 %then %put  WHEN CALLING THE BOOTSTRAP_RESULTS MACRO FOR CREATING THE FINAL RESULTS TABLE. THE COMPLETE RAW DATA CAN BE ;
                  %if &printlogstats = 1 %then %put  OBTAINED BY STACKING THE RESULTING DATA SETS FROM EACH PART. ;
            %end;
                    
            %if %bquote(&covmeandata)^= %then %do;
                %let _dot_ = %index(&covmeandata,%str(.)) ;
                %if &_dot_ > 0 %then   %let covmeanname0 = %cmpres(%qscan(&covmeandata,2,%str(.))) ;
                %else %let covmeanname0 = &covmeandata ;
                %let covmeanname = &covmeanname0._&sample_start._&sample_end ;
            %end;
            %else %do ;
                %PUT ERROR ;
                %put TO RUN THE GFORMULA WITH CHECK_COV_MODELS = 1 AND ONLY PART OF THE BOOTSTRAP ;
                %PUT SAMPLES YOU NEED TO SUPPLY COVMEANDATA TO BE NON-MISSING TO BE ABLE TO SAVE THE ;
                %PUT INTERMEDIATE RESULTS IN A CONSISTENT WAY.;
                %GOTO exit;
            %end;
       %end; 

       %if %bquote(&intervname)= %then %do;
             %PUT ERROR ;
             %put TO RUN THE GFORMULA WITH  ONLY PART OF THE BOOTSTRAP ;
             %PUT SAMPLES YOU NEED TO SUPPLY INTERVNAME TO BE NON-MISSING TO BE ABLE TO SAVE THE ;
             %PUT INTERMEDIATE RESULTS IN A CONSISTENT WAY. ;
             %GOTO exit;
       %end; 
       %else %do;
            %let intervdata = &savelib..&intervname._&sample_start._&sample_end ;
       %end;   
      
       %if &hazardratio = 1  %then %do;
                  
          %if %bquote(&hazardname) =  %then %do;
              %put TO RUN THE GFORMULA WITH  ONLY PART OF THE BOOTSTRAP ;
              %PUT SAMPLES YOU NEED TO SUPPLY HAZARDNAME TO BE NON-MISSING TO BE ABLE TO SAVE THE ;
              %PUT INTERMEDIATE RESULTS IN A CONSISTENT WAY. ;
              %GOTO exit;
           %end;
           %else %do;
               %let hazardname0 = &hazardname ;
               %let hazardname = &savelib..&hazardname._&sample_start._&sample_end;
           %end;
        %end;
    %end ;
    %else %do;
        %if %bquote(&survdata)= %then %let survdata = survdata_all ;
        %let covmeanname = &covmeandata ;
        %if &check_cov_models = 1 AND %bquote(&covmeandata) =   %then %let covmeanname   = _cov_mean_all ;
        %if &check_cov_models = 1 AND %bquote(&observed_surv) = %then %let observed_surv = _observed_surv_all ;
        %if &save_raw_covmean = 1 %then %let covmean_raw = &covmeanname._raw ;;
        %let hazardname = _inthr_ ;
                
    %end;
    
    
    /* set starting seed for run */

    %let seed_list = &seed ;
    %let seed_holder  =  &seed ;
    %do i = 1 %to &nsamples ;
       %let seed =%eval(&seed + 3) ;
       %let seed_list = &seed_list &seed ;
    %end;
   %let seed = %scan(&seed_list,%eval(&sample_start + 1 ),%str( ));

   %if &printlogstats = 1 %then %put  seed = &seed , seedlist = &seed_list ;   
   
   %do i = 1 %to %numargs(&outcinteract) ;
        %let uselabelo = 1 ;
        %local outc_I&i ;
   %end;
   %do i = 1 %to %numargs(&compeventinteract) ;
        %local compevent_I&i ;
   %end;
   %if &outctype = cateofu %then %do;
      proc sql noprint ;
                select max(&outc) as maxlev into :outclev from &data ;
            quit ;
            %let outcknots=%unquote(%makeknots5(&outclev));
            %if &printlogstats = 1 %then %put  Knots defining the categories for &outc : outcknots=&outcknots;

    %end;
    %do i = 0 %to  %eval(&ncov) ;     
       %local dimpred&i dimvar&i dimvar&i.z cov&i.lev cov&i.elev cov&i.min cov&i.max cov&i.array cov&i.zarray 
        cov&i.ninterx dimvar&i.class  ;
      

       %do j = 1 %to %eval(%numargs(&&cov&i.interact)) ;
           %let uselabelc = 1 ;
           %local cov&i._I&j ;
       %end;

	 /*  %put i = &i ,var =  &&cov&i , pytpe= &&cov&i.ptype , knots =  &&cov&i.knots , etype = &&cov&i.etype , eknots =  &&cov&i.eknots ;
       %put &&cov&i -- &&cov&i.etype -- %numargs(&&cov&i.eknots ) ;  */

       /* create cov&i.lev = number of levels (or knots) for ptype -cat  variables.*/
	   
       %if &&cov&i.otype ne 5 and ( %index(&&cov&i.ptype ,cat) > 0  or  %index(&&cov&i.etype, cat) > 0 )           
         %then %do;
		 
             /** set knots and lev for lag1cumavgcat and lag2cumavgcat terms ***/
		     /** THIS CODE IS DEPENDENT ON COV.KNOTS NOT BEING A SINGLE INTEGER FOR THE NUMBER OF PERCENTILES TO USE **/

             %local cov&i._cumavg_l1_knots cov&i._cumavg_l1_lev  ;
             %if &&cov&i.ptype = lag1cumavgcat or &&cov&i.ptype = lag2cumavgcat %then %do;
                %if %index(&&cov&i.knots,%str(:)) = 0 %then %do;
                    %let cov&i._cumavg_l1_knots = &&cov&i.knots ;
                    %let cov&i._cumavg_l1_lev = %eval(%numargs(&&cov&i._cumavg_l1_knots) + 1) ;
                %end ;
                %else %do;
                    %let cov&i._cumavg_l1_knots = %scan(&&cov&i.knots,2,%str(:));
                    %let cov&i._cumavg_l1_lev = %eval(%numargs(&&cov&i._cumavg_l1_knots) + 1);
                    %let cov&i.knots = %scan(&&cov&i.knots,1,%str(:));
                %end;
                %put (knots, lev) for &&cov&i.ptype  &&cov&i : ( &&cov&i.knots , %eval(%numargs(&&cov&i.knots) + 1) )  _cumavg   (&&cov&i._cumavg_l1_knots , &&cov&i._cumavg_l1_lev );
            %end;	
 
			%if %index(&&cov&i.ptype, cat) > 0 %then %do;
                %let cov&i.lev=%eval(%numargs(&&cov&i.knots)+1); /* in this case covXknots is a true list of knots */             	
            %end;
			/* User must supply knots for categorical etype variables. */
			%if %index(&&cov&i.etype, cat) > 0 %then %do;
                %let cov&i.elev=%eval(%numargs(&&cov&i.eknots)+1);              	 
            %end;		
             
             %if &printlogstats = 1 AND %index(&&cov&i.ptype, cat) > 0  %then %put  Number of categories of &&cov&i : cov&i.lev=&&cov&i.lev ;
			 %if &usehistory_eof = 1 AND &i > 0   %then %do;
			    %if ( %index(&&cov&i.etype, spl) > 0 or ( &&cov&i.etype in (cumsum cumsumnew cumavg cumavgnew )) ) %then %do;
				 	%if %numargs(&&cov&i.eknots) > 1 %then %let cov&i.elev = %eval(%numargs(&&cov&i.eknots)+1);
					%else %if %numargs(&&cov&i.eknots) = 1 and &&cov&i.eknots > 0 %then %do;
						%let cov&i.elev = %eval(&&cov&i.eknots + 1) ;
						%let anyfindeknots = 1 ;
						%let cov&i.findeknots = 1 ;
					%end;
				 	%if &printlogstats = 1  %then %put  Number of (e)levels of &&cov&i , &&cov&i.etype : cov&i.elev=&&cov&i.elev ;
				 %end ;
			 %end;
        %end;
	   %if &&cov&i.ptype = tsswitch1  %then %do;			
	        %if %numargs(&&cov&i.knots) > 1 %then %let cov&i.lev = %eval(%numargs(&&cov&i.knots));
			%else %do;
			    %if &&cov&i.knots = 0 %then %do;
					%let cov&i.lev = 1; /* for splines the loop is over 1 to lev - 2, and default is to incriment by +1 */
				%end;
				%else %do ;
  					%let cov&i.lev = &&cov&i.knots ;
					%let anyfindpknots = 1 ;
					%let cov&i.findknots = 1 ;
				%end;
			%end;
			%if &printlogstats = 1 %then %put  Number of categories of &&cov&i(tsswitch1) : cov&i.lev=&&cov&i.lev ;
	   %end;
       %if %index(&&cov&i.ptype, spl) > 0   %then %do;
              %if %numargs(&&cov&i.knots) > 1  %then %let cov&i.lev=%numargs(&&cov&i.knots);
			  %else %do ;
					%let cov&i.lev = &&cov&i.knots ;
					%let anyfindpknots = 1 ;
					%let cov&i.findknots = 1 ;
			  %end;
              %if &printlogstats = 1 %then %put  Number of levels of &&cov&i : cov&i.lev=&&cov&i.lev;         
       %end;
	   %if  %index(&&cov&i.etype, spl) > 0  AND &usehistory_eof = 1  %then %do;
              %if %numargs(&&cov&i.eknots) > 1 %then %let cov&i.elev=%numargs(&&cov&i.eknots);
			  %else %do;
				%let cov&i.elev = &&cov&i.eknots ;
				%let anyfindeknots = 1 ;
				%let cov&i.findeknots = 1 ;
              %end;
              %if &printlogstats = 1 %then %put  Number of (e)level for &&cov&i using &&cov&i.etype : cov&i.elev=&&cov&i.elev;         
       %end;       
	   %if &i > 0  AND &usehistory_eof = 1 %then %do;
	   		%if (%upcase(&&cov&i.etype) in CUMSUM CUMSUMNEW CUMAVG CUMAVGNEW )  %then %do; 
			  	%if %numargs(&&cov&i.eknots) > 1 %then %let cov&i.elev=%numargs(&&cov&i.eknots);
			  	%else %do;
					%let cov&i.elev = &&cov&i.eknots ;
					%let anyfindeknots = 1 ;
					%let cov&i.findeknots = 1 ;
				%end;
              	%if &printlogstats = 1 %then %put  Number of (e)level for &&cov&i using &&cov&i.etype : cov&i.elev=&&cov&i.elev; 
			%end;
			%if (%upcase(&&cov&i.etype) in CUMSUMCAT CUMSUMCATNEW CUMAVGCAT CUMAVGCATNEW )  %then %do; 
	   			%if %numargs(&&cov&i.eknots) > 1 %then %let       cov&i.elev = %eval(%numargs(&&cov&i.eknots)+1);
				%else %if %numargs(&&cov&i.eknots) = 1 %then %do;
					%let cov&i.elev = %eval(&&cov&i.eknots + 1) ;
					%let anyfindeknots = 1 ;
					%let cov&i.findeknots = 1 ;
				%end;
			 	%if &printlogstats = 1  %then %put  Number of (e)levels of &&cov&i , &&cov&i.etype : cov&i.elev=&&cov&i.elev ;
				
			%end;
	   %end;
        
       %if &&cov&i.otype=5 %then %do;
            %if &printlogstats = 1 %then %put  Otype 5?: &&cov&i.otype; 
             
            proc sql noprint ;
                select max(&&cov&i) as maxlev into :cov&i.lev from &data ;
            quit ;
            %let cov&i.knots=%unquote(%makeknots5(&&cov&i.lev));
            %if &printlogstats = 1 %then %put  Knots defining the categories for &&cov&i : cov&i.knots=&&cov&i.knots;

			%if &usehistory_eof = 1 %then %do;
				   %let cov&i.eknots = &&cov&i.knots ;
				   %if &&cov&i.etype ne none %then %do;
                             %let cov&i.etype = cat ;
				   		     %let cov&i.elev = %eval(%numargs(&&cov&i.eknots) + 1) ;

				   		%if &printlogstats = 1 %then %put SETTING &&cov&i otype 5 to have etype = cat and eknots = &&cov&i.eknots  and elev = &&cov&i.elev ;
				   %end ;

			%end;

			
      %end;
              
       %end;
   
    %*Preparing data;    
    %dataprep;
 
 
   %*return ;

    %if %bquote(&nparam)= %then %let nparam=&ssize;
    %if %bquote(&nsimul)= %then %let nsimul=&ssize;

    %if %eval(&nparam - &ssize ) > 0 %then %do;
          %put ERROR: NPARAM SHOULD BE LESS THAN THE SAMPLE SIZE ;
          %goto exit ;
   %end;

    %*Looping over bootstrap samples;
    %do bsample = &sample_start %to &sample_end;
        %if (&outputs ^= yes or %eval(&bsample) ^= 0) %then %do;
               %let ods_logit = ods select none ;
               %let ods_reg = ods select none ;
               %let ods_qlim = ods select none ;
        %end;
        %*Generating samples;

%if &printlogstats = 1 %then %put  before sample = &bsample seed = &seed ;
%if &printlogstats = 1 %then %put  ;

        %samples;    
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

            %if &bsample = 0 AND  %bquote(&simuldata)^=  %then %do;

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

				%if &outctype = cateofu %then %do;
				    title "Frequency of &outc under intervention 0";
					proc freq data = &simuldata ;
					table &outc / out = interv0 (keep = &outc percent ) nocum  ;
					run;

					
					proc transpose data = interv0 out = interv0_&outc  (drop = _name_ _label_)  prefix = smybinom_ ;
					var percent ;
					id mybinom ;
					run;

					

                %end;

                proc means data= &simuldata  mean min max ;                                                                      
                var %if &outctype = binsurv %then cuminc ;
                    intervened averinterv  
                    %if &outctype = binsurv %then %do;
                        %do j = 1 %to %eval(&timepoints);
                            s&outc.&j
                        %end;
                    %end;
                    %else %if &outctype ne cateofu %then %do ;
                          &outc 
					%end;
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
                %if &outctype ne cateofu %then %do;
                	set interv0;
				%end;
				%else %do;
					merge interv0 interv0_&outc ;
                %end; 
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
				%else %if &outctype = cateofu %then %do;
					%do j = 1 %to &outclev ;
					   s&outc._&j
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
                        keep  int int2 _sample_ %if &outctype = cateofu %then %do;
					                                %do j = 1 %to &outclev ;
					                                   s&outc._&j
					                                %end;
				                              %end;
                                              %else  s&outc  ;;
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
   
    %end;


     %if &runnc = 1 OR (&runnc = 0 AND &numint > 0 AND &refint > 0 ) %then %do;
        %if &chunked = 0 %then %do;
           %*Summarizing final results;
           %results;
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
    proc datasets library = work memtype = (data view)  nolist ;
    delete tmpids _calchazard_ %if &runnc = 1 %then simulated0 ;
            %do i = 1 %to &numint ; simulated&i %end;
            %if  &outctype = bineofu or &outctype=conteofu or &outctype=conteofu2 or 
                   &outctype = conteofu3 or &outctype = conteofu4 or &outctype = cateofu %then &data ;
         ;
    quit;

   /* clean up any global variables created in interaction macros */
    %remove_created_global ;

    %*Exiting in case of error;
    %exit:;
   
     
%mend gformula;
/***/
%macro createfixedcov ;
   %local i j k fixedcov1 fixedcov2  fixedcov3  covlist word testword nfixed incovlist ;
   %do i = 1 %to &ncov ;
       %let covlist = &covlist &&cov&i ;
   %end;
   %if &printlogstats = 1 %then %put  covlist = &covlist ;

   %let nfixed = %numargs(&fixedcov);

   %let fixedcov1 = ;
   %let fixedcov2 =  ;
   


   %do i = 1 %to &nfixed ;
       %let incovlist = 0 ;
       %let word = %scan(&fixedcov,&i);
       %do j = 1 %to &ncov ;
            %let testword = %scan(&covlist,&j);
            
            %if %upcase(&word) = %upcase(&testword) %then  %do; 
            /* when a baseline variable is also one of the modeled covariates we need to create the interaction with the 
               desired indicator function to delay the inclusion into the model. This is only for the model of cov&j */

                  %let incovlist = 1;
                  %let cov&j.infixed = 1 ;
                  %let ptype = %substr(&&cov&j.ptype,1,4);
                  
                  %if %upcase(&ptype)=LAG1 or %upcase(&ptype) = CUMA or %upcase(&ptype) = RCUM %then %let cov&j.interact = &&cov&j.interact B&i.*L2 ;
                  %else %if %upcase(&ptype)=LAG2 %then %let cov&j.interact = &&cov&j.interact B&i.*L3 ;
                  %else %if %upcase(&ptype)=LAG3 %then %let cov&j.interact = &&cov&j.interact B&i.*L4 ;                                    
 
            %end;
       %end;
       %if &incovlist = 1 %then %do;           
         %let fixedcov1 = &fixedcov1 &word._b ;                
       %end; 
       %else  %do;
           %let fixedcov1 = &fixedcov1 &word ;                   
       %end;
  %end;
  /* for models, fixedcov may change as follows: for outcome types, fixedcov will be */
  %let fixedcov =  &fixedcov1 ;/* for outcome type models */
  
  %do i = 1 %to &ncov ;
       %if &&cov&i.infixed = 0 %then %let cov&i.fixedcov = &fixedcov ; 
       %else %do ;
             /* baseline variable will be included in the interaction terms */
             %do j = 1 %to &nfixed ;
                   %let word = %scan(&fixedcov,&j);
                   %if %upcase(&word) ^= %upcase(&&cov&i)_B %then 
                          %let cov&i.fixedcov = &&cov&i.fixedcov &word ;
             %end;
      %end;
     
     
  %end; 
           
%mend;
/***/
%macro createcov0vars   ;
  %let cov0= &time ; 
  %let cov0otype=0; 
  %let cov0ptype= &timeptype; 
  %let cov0etype = ;
  %let cov0cumint = ;
  %let cov0mtype= all ; 
  %let cov0skip=-1; 
  %let cov0inc=&timeinc; 
  %let cov0knots=&timeknots; 
  %let cov0eknots = ;
  %let cov0lev= ; 
  %let cov0interact= ; 
  %let cov0wherem=(1=1); 
  %let cov0wherenosim=(1=0); 
  %let cov0nosimelsemacro=; 
  %let cov0class=; 
  %let cov0classelse=; 
  %let cov0addvars=; 
  
  %let cov0genmacro= &timefuncgen; 
  %let cov0modusermacro =; 
  %let cov0moddatausermacro =; 
  %let cov0setblvar =; 
  %let cov0simusermacro = ; 
  %let cov0barray = ; 
  %let cov0sarray =; 
  %let cov0randomvisitp=; 
  %let cov0visitpmaxgap=9e10; 
  %let cov0visitpwherem=(1=1); 
  %let cov0visitpcount= ; 
  %let cov0visitpelse = ;

%mend ;
%macro dataprep;
   /*Creating product ("interaction") terms*/
               
  %local i iii dsid rc droplist  ;
  %let dsid = %sysfunc(open(&data ));
  %do i = 1 %to &ncov ;
      %if %sysfunc(varnum(&dsid,&&cov&i.._b)) %then %let droplist = &droplist &&cov&i.._b ;
  %end;
   
  %let rc = %sysfunc(close(&dsid)) ;

  %if  &outctype = bineofu or &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype = conteofu4 or &outctype = cateofu %then %do;
       data &data._eof ;
	   set &data ;
	   if &time >= &timepoints then delete ; * will use upto time = timepoints - 1 ;
	   run;

	   %let data = &data._eof ;
  %end ;	   

   data _null_(drop=);
   set &data;

   %interactions(outc,,,createlist );
   %if %bquote(&compevent)^= %then %do;
          %interactions(compevent, ,,createlist);
   %end;
   %if %bquote(&censor)^= %then %do;
          %interactions(censor, ,,createlist);
   %end;
   
   
   %do iii=0 %to &ncov;
        %if &&cov&iii.otype^=0 %then %do;
            
            %interactions(cov&iii, ,,createlist);
        %end;
        %else %do;
            %let cov&iii.ninterx=0;
        %end;
        %let covlist = &covlist &&cov&iii ;
   %end;
   
   run;
   
    %*Arraying regression predictors;
     
    %if &checkaddvars = 1 %then %addvarcheck(vartype = 0, outvar = event );

	 /* for eof type models starting listpred at 0 includes time. For eof type models, time is not needed. We could start the loop at 1
	   instead. Would need to check the use of listpred in the interv_init macro for setting up the arrays ***/
	%if &usehistory_eof = 0 %then %do;	     
    		%let outcpred = &fixedcov %listpred(main,,0,&ncov,_l2,_l1) &eventaddvars %qcmpres(%interxarrays(main,outc))   ;
	%end ;
	%else %if &usehistory_eof = 1 %then %do; 
		%let outcpred = &fixedcov %listpred(eof,,0,&ncov,_l2,_l1) &eventaddvars %qcmpres(%interxarrays(main,outc))   ;
	%end ;
    %if &printlogstats = 1 %then %put  Outcome predictors are: &outcpred;

    %if %bquote(&compevent)^=  %then %do;
       %if &checkaddvars = 1 %then %addvarcheck(vartype = 0, outvar = compevent );
       %let compeventpred = &fixedcov %listpred(main,,0,&ncov,_l2,_l1)  &compeventaddvars %qcmpres(%interxarrays(main,compevent))   ;
       %if &printlogstats = 1 %then %put  Censoring by death predictors are: &compeventpred;
    %end;  

	%if %bquote(&censor)^=  %then %do;
       %if &checkaddvars = 1 %then %addvarcheck(vartype = 0, outvar = censor );
       %let censorpred = &fixedcov %listpred(main,,0,&ncov,_l2,_l1)  &censoraddvars %qcmpres(%interxarrays(main,censor))   ;
       %if &printlogstats = 1 %then %put  Censoring by lost-to-followup predictors are: &censorpred;
    %end;   
    


   %do i = 0 %to &ncov;

      %if &checkaddvars = 1 %then %addvarcheck;
      %let changemtype = 0 ;
      %if %upcase(&&cov&i.mtype)  = SOME %then %do;
             %let changemtype = 1 ;
             %let cov&i.mtype = all ; /* need this to include the correct form of the lagged values of cov&i in model for cov&i */
      %end; 
       %if &printlogstats = 1 %then %put  cov&i.ninterx=&&cov&i.ninterx;
       %let cov&i.array =   &&cov&i.fixedcov %listpred(main,&i,0,&ncov,_l3,_l2,_l1) %listpred(contemp,,0,&i-1) &&cov&i.addvars 
        %qcmpres(%interxarrays(main,cov&i) ) ;
        
    
       %if &&cov&i.otype^=0 %then %do;
          %if &printlogstats = 1 %then %put  cov&i (&&cov&i) predictors are: &&cov&i.array;
          %if &printlogstats = 1 %then %put ;       
       %end;
       %if &changemtype = 1 %then %let cov&i.mtype =  some ;
    %end;
   
    
 
    %************ PREPARING DATA;  

    %*Calculating obs pr of outcome (by time points) if survival time analysis;
    %if &outctype=binsurv %then %do;
       %obscuminc(data=&data,time=&time,timepoints=&timepoints,event=&outc, compevent=&compevent);
    %end; 

   
    %*Calculating observed mean of outcome (at end of follow-up) if continuous or dichotomous outcome analysis;
    %else %do;
       %obstotcont(data=&data,time=&time,timepoints=&timepoints,outc=&outc, compevent=&compevent);
    %end; 
    

    %*Generating _orig_ dataset with ID number, a new sequential ID and sample (=0);
    %*Also generating macros for obsp, ssize, min, max, and dim arrays;
    proc sort data=&data nodupkey out=indiv;
        by &id;
    run;


    data _orig_ (keep= _sample_ newid &id);
	   %if &outctype ne cateofu %then %do;
        	merge indiv cuminc  end=_end_;
	   %end;
	   %else %do;
	   	  set indiv end= _end_ ;
	   %end;

        newid=_n_;
        _sample_=0;
                
        if _n_=1 then do;

            array aoutc  &outcpred;
            dimoutc=dim(aoutc);
            call symput('dimoutc',trim(left(put(dimoutc,8.))) );
            
           %if %bquote(&compevent)^= %then %do;  
             array acompevent  &compeventpred;
             dimcompevent=dim(acompevent);
             call symput('dimcompevent',trim(left(put(dimcompevent,8.))) );
           %end;           
        

           %do i = 0 %to &ncov;
                
                array avar&i &&cov&i.array;
                dimvar&i = dim(avar&i);
                call symput('dimvar'||(left(put(&i,8.))),trim(left(put(dimvar&i,8.))) );
                    
             %end;
                        
            %if &outctype ne cateofu %then call symput('obsp',trim(left(cuminc)));;
            
        end;
      
        if _end_ then do;
            call symput('ssize',trim(left(newid)));
        end;
    run;
    

    %*Adds sequential ID to the data and sorts on that;
    proc sort data=_orig_;
        by &id;
    run;

    proc sort data=&data out=_inputd_;
    by &id descending &time;
    run;

    data _inputd_;
    merge _orig_ _inputd_;
    by &id;
    %if %bquote(&droplist)^= %then drop &droplist ;;
    run;

    proc sort data = _inputd_ ;
    by newid &time;
    run;

    %* Generates macros for min and max of covariates and outcome if continuous;
    data _covbounds_;
        set _inputd_ end=_end_;
        %do i=0 %to &ncov;
            %if &&cov&i.otype=3 or &&cov&i.otype=4 or &&cov&i.otype=6 or &&cov&i.otype=7   or &&cov&i.otype = -1 %then %do;
                
                retain &&cov&i.._min &&cov&i.._max;

                if _n_ = 1 then do;
                    &&cov&i.._min =  1.0e100;
                    &&cov&i.._max = -1.0e100;
                    end; 
                
                %if &&cov&i.otype ^= 4 %then if &&cov&i ne . and &&cov&i < &&cov&i.._min then &&cov&i.._min = &&cov&i;;
                %if &&cov&i.otype  = 4 %then if &&cov&i > 0  and &&cov&i < &&cov&i.._min then &&cov&i.._min = &&cov&i;;
                
                if &&cov&i ne . and &&cov&i > &&cov&i.._max then &&cov&i.._max = &&cov&i;
                
                if _end_ then do; 
               
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
		%else %if &outctype=conteofu4    %then %do;

             retain &outc._min &outc._max;

                if _n_ = 1 then do;
                    &outc._min =  1.0e100;
                    &outc._max = 0 ;
                    end; 

                if &outc > 0 and &outc < &outc._min then &outc._min = &outc;
                if &outc > 0 and &outc > &outc._max then &outc._max = &outc;
                
                if _end_ then do;
                    call symput("outcmin",trim(left(&outc._min)) );
                    call symput("outcmax",trim(left(&outc._max)) );
                    end;
        %end;
        _sample_ = 0 ;
        keep  _sample_ %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype = conteofu4 %then  &outc._min &outc._max ;
                %do i = 0 %to &ncov ; 
                    %if &&cov&i.otype=3 or &&cov&i.otype=4 or &&cov&i.otype=6 or &&cov&i.otype=7   or &&cov&i.otype = -1  %then &&cov&i.._min &&cov&i.._max ;
                %end ;
                 ;
        if _end_ then output ;
    run;

    %do i = 0 %to &ncov;
        %if &&cov&i.otype=3 or &&cov&i.otype=4 or  &&cov&i.otype=6 or &&cov&i.otype=7  or &&cov&i.otype=-1 %then %do;
            %if &printlogstats = 1 %then %put  &&cov&i range is &&cov&i.min to &&cov&i.max;                    
            %end;
        %end;
    %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype = conteofu4 %then %do;
      %if &printlogstats = 1 %then %put  &outc range is &outcmin to &outcmax;
      %end;
   

    %*Adding necessary variables for models;
    data _inputd_;
        set _inputd_ ;
        by newid &time;  
        
        &time._l1 = &time-1;
        &time._l2 = &time-2;
        &time._l3 = &time-3;
    
        retain _ind2 _ind3 _ind4  ;
        if first.newid then do ;
            _ind2 = 0 ;
            _ind3 = 0 ;
            _ind4 = 0 ;
        end;
        if &time = 2 then _ind2 = 1 ;
        if &time = 3 then _ind3 = 1 ;
        if &time = 4 then _ind4 = 1 ;
 
        %do i=0 %to &ncov;
              retain &&cov&i.else ;
              %if &&cov&i.otype = 0 %then %do;
                   &&cov&i.._l1 = &&cov&i - 1 ;
                   &&cov&i.._l2 = &&cov&i - 2 ;
                   &&cov&i.._l3 = &&cov&i - 3 ;
              %end ;
              if &time = 0 then &&cov&i.else = &&cov&i ;
              %if &&usevisitp&i = 1 %then %do;
                   retain &&cov&i.visitpelse ;
                   if &time = 0 then &&cov&i.visitpelse = &&cov&i.randomvisitp ;
              %end; 
              %genpred(main);
       %end;


%do iii = 0 %to &ncov ;
 
    %if &&usevisitp&iii = 1 %then %do;
        retain ts_last_&&cov&iii ;
        if first.newid then do ;
            ts_last_&&cov&iii.._l1 = &&cov&iii.visitpcount   ; /* do not need to subtract 1 since visitpcount is the lagged value already */ 
            ts_last_&&cov&iii = &&cov&iii.visitpcount   ;           
        end;
        if &&cov&iii.randomvisitp = 1 then ts_last_&&cov&iii = 0; 
        else ts_last_&&cov&iii = ts_last_&&cov&iii + 1 ;
        lastvisitholder = lag(ts_last_&&cov&iii) ;
        if first.newid = 0 then ts_last_&&cov&iii.._l1 = lastvisitholder ;
        drop lastvisitholder ;   
       
    %end;

%end;

   %if %bquote(&weight) ^= %then %do ;
       _weight_ = &weight ;
   %end;
   %else %do;
      _weight_ = 1.0 ;
   %end;

   %if &usehistory_eof = 1 %then %do ;
		%do i = 1 %to &ncov ;
		    retain %do timeindex = 0 %to %eval(&timepoints - 1) ;
                   		&&cov&i.._&timeindex._eof 
                   %end ;
				    ;
            array a&&cov&i{0:%eval(&timepoints - 1) } %do timeindex = 0 %to %eval(&timepoints - 1) ;
     														&&cov&i.._&timeindex._eof
												      %end;
													   ;
            a&&cov&i [&time ] = &&cov&i ;

			%genpred(type=main,useeof = 1 ) ;
			     
        %end;
          

   %end;

    run;
/***
%put ANYFINDPKNOTS = &anyfindpknots ;
%put ANYFINDEKNOTS = &anyfindeknots ;
****/

/* anyfindpknots is an indicator if we need to find the knots based on the percentiles of at least one variable. In this case, it should be that 
   covXpknots is in {3,4,5} for the number of percentiles to use.

   In the loop below, numberofknots is the number of knots being generated. Will need this for the case where anyfindpknots = 0 so that all knots are being
   specified in the covXknots variable.
     
*/

 %if /* ((&anytsswitch1 = 1 OR &anylagcumavg = 1  ) ) OR */ &anyfindpknots = 1  %then %do;

      %do i = 1 %to &ncov ;
           %local numberofknots&i ;
           %let numberofknots&i = 0 ; 
/* %put &&cov&i , &&cov&i.ptype , &&cov&i.knots , &&cov&i.findknots  ; */
           %if %upcase(&&cov&i.ptype) = TSSWITCH1 AND /* &usespline = 1 */  &&cov&i.findknots = 1 /*%quote(&&cov&i.knots) ^= 0 */  %then %do;

		            %let numberofknots&i = &&cov&i.knots ;
				/*	%let numberofknots&i = 4 ; */ /* temporarily hard code to use 4 percentiles as originally coded */
                /*** EM edited ***/
                    proc univariate data = _inputd_ (where = (&&cov&i = 1))  noprint  ;
                    var ts&&cov&i.._inter ;
                    output out = tsscov_pct   
                         pctlpre =  tsscov  
                         %if &&numberofknots&i = 3 %then 
                            pctlname =  _pct10 _pct50 _pct90 
                             pctlpts =   10 50 90 ;
                          %else %if &&numberofknots&i = 4 %then
                                 pctlname = _pct5 _pct25 _pct75 _pct95 
                                 pctlpts =  5 25 75 95 ;
                           %else %if  &&numberofknots&i = 5 %then
                                 pctlname = _pct5 _pct25 _pct50 _pct75 _pct95 
                                 pctlpts =  5 25 50 75 95 ;
                          ;
                    run;

                    proc transpose data = tsscov_pct out = ttsscov_pct  (keep = col1);
                    run;
                    proc sql  noprint ;
                    select col1 into : cov&i.knots separated by ' ' from ttsscov_pct ;
                    quit ;
                    %if &printlogstats = 1 %then %put knots for &&cov&i.._inter  = &&cov&i.knots ;
                    proc datasets library = work nolist ;
                    delete tsscov_pct ttsscov_pct ;
                    quit;

            %end;
			%if (%upcase(&&cov&i.ptype) in LAG1SPL LAG2SPL LAG3SPL SKPSPL ) AND /* &usespline = 1 */  &&cov&i.findknots = 1 /*%quote(&&cov&i.knots) ^= 0 */  %then %do;

		            %let numberofknots&i = &&cov&i.knots ;
				/*	%let numberofknots&i = 4 ; */ /* temporarily hard code to use 4 percentiles as originally coded */
                    /*** EM edited ***/
                    proc univariate data = _inputd_   noprint  ;
					%if &&cov&i.ptype = skpspl %then where &time not in (&&cov&i.skip) ;;
                    var &&cov&i ;
                    output out = &&cov&i.._pct   
                         pctlpre =  &&cov&i  
                         %if &&numberofknots&i = 3 %then 
                            pctlname =  _pct10 _pct50 _pct90 
                             pctlpts =   10 50 90 ;
                          %else %if &&numberofknots&i = 4 %then
                                 pctlname = _pct5 _pct25 _pct75 _pct95 
                                 pctlpts =  5 25 75 95 ;
                           %else %if  &&numberofknots&i = 5 %then
                                 pctlname = _pct5 _pct25 _pct50 _pct75 _pct95 
                                 pctlpts =  5 25 50 75 95 ;
                          ;
                    run;

                    proc transpose data = &&cov&i.._pct out = t&&cov&i.._pct  (keep = col1);
                    run;
                    proc sql  noprint ;
                    select col1 into : cov&i.knots separated by ' ' from t&&cov&i.._pct ;
                    quit ;
                    %if &printlogstats = 1 %then %put knots for &&cov&i  = &&cov&i.knots ;
                    proc datasets library = work nolist ;
                    delete &&cov&i.._pct t&&cov&i.._pct ;
                    quit;

            %end;


            %if ( %upcase(&&cov&i.ptype) in  CUMAVG LAG1CUMAVG LAG2CUMAVG ) and /*&usespline*/ &&cov&i.findknots = 1 %then %do;
                      %local cumavgwhere cumavgvar numberofknots ;
                      %let numberofknots = &&cov&i.knots ;
                      %let cumavgwhere = ;
                      %let cumavgvar = &&cov&i.._cumavg_l1 ;
                      %if %bquote(&&cov&i.cumint)^= %then %do;
                          %let cumavgwhere = (where = (&&cov&i.cumint = 1)) ;
                          
                      %end;
                                              

                    %if &numberofknots > 0  %then %do;  
                    /*** EM edited ***/
                    proc univariate data = _inputd_   %if &&cov&i.otype = 4  %then (where = (&&cov&i > 0)) ; noprint  ;
                    var %if %upcase(&&cov&i.ptype) ne CUMAVG  %then &&cov&i   ;
                        %else  &&cov&i.._cumavg ; ;
                    output out = &&cov&i.._pct   
                         pctlpre =  &&cov&i  
                         %if &numberofknots = 3 %then 
                            pctlname =  _pct10 _pct50 _pct90 
                             pctlpts =   10 50 90 ;
                          %else %if &numberofknots = 4 %then
                                 pctlname = _pct5 _pct25 _pct75 _pct95 
                                 pctlpts =  5 25 75 95 ;
                           %else %if  &numberofknots = 5 %then
                                 pctlname = _pct5 _pct25 _pct50 _pct75 _pct95 
                                 pctlpts =  5 25 50 75 95 ;
                          ;
                      run;

                     proc transpose data = &&cov&i.._pct out = &&cov&i.._pct  (keep = col1);
                     run;

                     proc sql noprint  ;
                     select col1 into : cov&i.knots separated by ' ' from &&cov&i.._pct ;
                     quit ;
                    %if &printlogstats = 1 %then %put knots for &&cov&i  = &&cov&i.knots ;
   
                    /*** EM edited ***/
                    %if %upcase(&&cov&i.ptype) ne CUMAVG  %then %do;
                        proc univariate data = _inputd_  &cumavgwhere   noprint  ;
                        var  &cumavgvar ;
                        output out = cumavg_l1_pct   
                             pctlpre =  cumavg_l1
                              %if &numberofknots = 3 %then 
                                pctlname =  _pct10 _pct50 _pct90
                                pctlpts =   10 50 90 ;
                              %else %if &numberofknots = 4 %then
                                     pctlname = _pct5 _pct25 _pct75 _pct95 
                                     pctlpts =  5 25 75 95 ;
                               %else %if  &numberofknots = 5 %then
                                     pctlname = _pct5 _pct25 _pct50 _pct75 _pct95 
                                     pctlpts =  5 25 50 75 95 ;
                              ;
                         
                         
                        run;

                        proc transpose data = cumavg_l1_pct out = cumavg_l1_pct  (keep = col1);
                        run;

                        proc sql noprint  ;
                        select col1 into : cov&i._cumavg_l1_knots separated by ' ' from cumavg_l1_pct ;
                        quit ;

                        %if &printlogstats = 1 %then %put knots for &cumavgvar  = &&cov&i._cumavg_l1_knots ;
                    %end ;
                    proc datasets library = work nolist ;
                    delete &&cov&i.._pct  %if %upcase(&&cov&i.ptype) ne CUMAVG %then cumavg_l1_pct ; ;
                    quit;
                %end;
            %end;

		
       %end;

     %end;

	 %if  ((( &anyecumavg = 1 /* AND &usespline = 1 */) OR &anyecumcat = 1 ) AND &usehistory_eof = 1) OR &anyfindeknots = 1 %then %do;



      %do i = 1 %to &ncov ; /* etype is only for non-time variables */
        

		 %if &&cov&i.etype in cumavg /* cumavgcat */ cumsum /* cumsumcat */  %then %do;
			     %if &&cov&i.etype_part2 = &timepoints %then %do;
					%let timeindex = 0;		            
				 %end;
				 %else %do;
				 	%let timeindex = %eval(&timepoints - &&cov&i.etype_part2 ) ;		            
				 %end;
			%end;
			%if &&cov&i.etype in cumavgnew cumavgcatnew cumsumnew cumsumcatnew  %then %do;
	     		%if &&cov&i.etype_part2 = &timepoints %then %do;
					%let timeindex = 0;            	 
			 	%end;
		 		%else %do;
		 			%let timeindex = %eval(&timepoints - &&cov&i.etype_part2 - 1 ) ;            	 
		 		%end;
			%end;	
           

			%local numberofeknots&i ;
			%let numberofeknots&i = 0 ;
           

			%if (&&cov&i.etype in spl skpspl) AND  &&cov&i.findeknots = 1  %then %do;

		            %let numberofeknots&i = &&cov&i.eknots ;
				/*	%let numberofknots&i = 4 ; */ /* temporarily hard code to use 4 percentiles as originally coded */
                    /*** EM edited ***/
                    proc univariate data = _inputd_   noprint  ;
					%if &&cov&i.etype = skpspl %then where &time not in (&&cov&i.skip) ;;
                    var &&cov&i ;
                    output out = &&cov&i.._pct   
                         pctlpre =  &&cov&i  
                         %if &&numberofeknots&i = 3 %then 
                            pctlname =  _pct10 _pct50 _pct90
                             pctlpts =   10 50 90 ;
                          %else %if &&numberofeknots&i = 4 %then
                                 pctlname = _pct5 _pct25 _pct75 _pct95 
                                 pctlpts =  5 25 75 95 ;
                           %else %if  &&numberofeknots&i = 5 %then
                                 pctlname = _pct5 _pct25 _pct50 _pct75 _pct95 
                                 pctlpts =  5 25 50 75 95 ;
                          ;
                    run;

                    proc transpose data = &&cov&i.._pct out = t&&cov&i.._pct  (keep = col1);
                    run;
                    proc sql  noprint ;
                    select col1 into : cov&i.eknots separated by ' ' from t&&cov&i.._pct ;
                    quit ;
                    %if &printlogstats = 1 %then %put eknots for &&cov&i  = &&cov&i.eknots ;
                    proc datasets library = work nolist ;
                    delete &&cov&i.._pct t&&cov&i.._pct ;
                    quit;

            %end;



            %if ( &&cov&i.etype in  cumsum cumsumnew  /* cumsumcat cumsumcatnew */ cumavg cumavgnew /* cumavgcat cumavgcatnew */   )  %then %do;
			   %if /*%numargs(&&cov&i.eknots ) = 1  */  &&cov&i.findeknots = 1 %then %do ;
                      
                      %let numberofeknots&i = &&cov&i.eknots ;   
 
                    %if &&numberofeknots&i > 0  %then %do;  
	                    /*** EM edited ***/
                        proc univariate data = _inputd_   (where = (&time = (%eval(&timepoints - 1))))  noprint  ;
	                    var %if &&cov&i.etype in cumsum cumsumnew cumsumcat cumsumcatnew  %then &&cov&i.._cumsum_&timeindex._eof   ;
	                        %else                                                               &&cov&i.._cumavg_&timeindex._eof  ; ;
	                    output out = &&cov&i.._pct   
	                         pctlpre =  &&cov&i  
	                         %if &&numberofeknots&i = 3 %then 
	                            pctlname =  _pct10 _pct50 _pct90 
                                pctlpts =   10 50 90 ;
	                          %else %if &&numberofeknots&i = 4 %then
	                                 pctlname = _pct5 _pct25 _pct75 _pct95 
	                                 pctlpts =  5 25 75 95 ;
	                           %else %if  &&numberofeknots&i = 5 %then
	                                 pctlname = _pct5 _pct25 _pct50 _pct75 _pct95 
	                                 pctlpts =  5 25 50 75 95 ;
	                          ;
	                      run;

	                     proc transpose data = &&cov&i.._pct out = &&cov&i.._pct  (keep = col1);
	                     run;

	                     proc sql noprint  ;
	                     select col1 into : cov&i.eknots separated by ' ' from &&cov&i.._pct ;
	                     quit ;

	                    %if &printlogstats = 1 %then  %put eknots for &&cov&i , &&cov&i.etype  : &&cov&i.eknots ;
	   
	                    

	                    proc datasets library = work nolist ;
	                    delete &&cov&i.._pct   ;
	                    quit;
                   %end; /* number of knots > 0 */
				%end; /* need to calculate the knots using percentiles, covIknots is a single entry */
            %end; /* etype is one of the cumsum type variables */


		

       %end; /* end of loop over i */

    %end;  /* end of condition for usespline, anyecumavg */


	%if &usehistory_eof = 0 %then %do;
	   %do i = 1 %to &ncov ;
           %local numberofeknots&i numberofknots&i   ;
		   %let numberofknots&i = %numargs(&&cov&i.knots) ;
		   %let numberofeknots&i = 0 ;
	   %end;
	%end;
	%else %if &usehistory_eof = 1  /*AND &anyfindpknots = 0 */  %then %do;
		%do i = 1 %to &ncov ;
			%local numberofknots&i numberofeknots&i ;		
			%let numberofknots&i = %numargs(&&cov&i.knots) ;
			%let numberofeknots&i = %numargs(&&cov&i.eknots );		
		%end;
	%end;
  

       data _inputd_;
       set _inputd_;
       %do i = 1 %to &ncov ;
       %put &&cov&i --- &&cov&i.ptype --- &&numberofknots&i -- &&cov&i.etype --- &&numberofeknots&i ;
         %if (%upcase(&&cov&i.ptype) = TSSWITCH1  ) AND  &&numberofknots&i >= 3 %then %do;
              %rcspline(ts&&cov&i.._inter ,&&cov&i.knots);
              %rcspline(ts&&cov&i.._l1_inter ,&&cov&i.knots);
         %end;

			%if ( &&cov&i.ptype in lag1cat lag2cat lag3cat  ) AND &&numberofknots&i > 0 %then %do;
				%makecat(&&cov&i ,&&cov&i.knots , &&cov&i.lev);
				%makecat(&&cov&i.._l1 ,&&cov&i.knots , &&cov&i.lev);
				%if &&cov&i.ptype = lag2cat %then %makecat(&&cov&i.._l2, &&cov&i.knots, &&cov&i.lev ) ;
				%if &&cov&i.ptype = lag3cat %then %makecat(&&cov&i.._l3, &&cov&i.knots, &&cov&i.lev ) ;

			%end;
			%if ( &&cov&i.ptype in lag1spl lag2spl lag3spl  ) AND &&numberofknots&i > 0 %then %do;
				%rcspline(&&cov&i ,&&cov&i.knots );
				%rcspline(&&cov&i.._l1 ,&&cov&i.knots );
				%if &&cov&i.ptype = lag2spl %then %rcspline(&&cov&i.._l2, &&cov&i.knots ) ;
				%if &&cov&i.ptype = lag3spl %then %rcspline(&&cov&i.._l3, &&cov&i.knots ) ;

			%end;

			%if ( &&cov&i.ptype in skpcat  ) AND &&numberofknots&i > 0 %then %do;
				%makecat(&&cov&i ,&&cov&i.knots , &&cov&i.lev);
				%makecat(&&cov&i.._l1 ,&&cov&i.knots , &&cov&i.lev);
				%do lev = 1 %to %eval(&&cov&i.lev - 1);
               		%if &current = 1 %then  %maketi(&&cov&i.._&lev,&time,&time._l1,&&cov&i.skip, &interval);
               		%if &lagged = 1 %then   %maketi(&&cov&i.._l1_&lev,&time._l1,&time._l2,&&cov&i.skip, &interval);
				%end;
            %end;

			%if ( &&cov&i.ptype in skpspl  ) AND &&numberofknots&i > 0 %then %do;
				%rcspline(&&cov&i ,&&cov&i.knots );
				%rcspline(&&cov&i.._l1 ,&&cov&i.knots );
				%do knot = 1 %to %eval(&&cov&i.lev - 2);
               		%if &current = 1 %then  %maketi(&&cov&i.._spl&knot,&time,&time._l1,&&cov&i.skip, &interval);
               		%if &lagged = 1 %then    %maketi(&&cov&i.._l1_spl&knot,&time._l1,&time._l2,&&cov&i.skip, &interval);
            	%end;
			%end;


         %if %upcase(&&cov&i.ptype) = CUMAVG or %upcase(&&cov&i.ptype) = LAG1CUMAVG or %upcase(&&cov&i.ptype) = LAG2CUMAVG  
                       %then %do;
                %let cumavgvar = &&cov&i.._cumavg_l1 ;
                %if %bquote(&&cov&i.cumint)^= %then %let cumavgvar = &&cov&i.._avg_l1 ;                

              %if %numargs(&&cov&i.knots) > 2 %then %do;
                %if %upcase(&&cov&i.ptype) ne CUMAVG %then %do;
                    %rcspline(&&cov&i ,&&cov&i.knots);
                    %rcspline(&&cov&i.._l1 ,&&cov&i.knots);
                    %rcspline(&&cov&i.._l2 ,&&cov&i.knots);
                 
                    %rcspline(&&cov&i.._cumavg_l1 ,&&cov&i._cumavg_l1_knots);
                    %rcspline(&&cov&i.._cumavg_l2 ,&&cov&i._cumavg_l1_knots);
                    %rcspline(&&cov&i.._cumavg_l3 ,&&cov&i._cumavg_l1_knots);
                %end;
                %else %do;
                    %rcspline(&&cov&i.._cumavg ,    &&cov&i.knots );
                    %rcspline(&&cov&i.._cumavg_l1 , &&cov&i.knots);
                 %end;
             %end;

          %end;

          %if %upcase(&&cov&i.ptype) = CUMAVGCAT %then %do;
            %makecat(&&cov&i.._cumavg, &&cov&i.knots, &&cov&i.lev);
            %makecat(&&cov&i.._cumavg_l1 , &&cov&i.knots, &&cov&i.lev);

          %end;

        
       %if %upcase(&&cov&i.ptype) = LAG1CUMAVGCAT %then %do;
            %makecat(&&cov&i , &&cov&i.knots , &&cov&i.lev);
            %makecat(&&cov&i.._l1, &&cov&i.knots, &&cov&i.lev);
            %makecat(&&cov&i.._cumavg_l1,&&cov&i.knots, &&cov&i.lev);
            %makecat(&&cov&i.._cumavg_l2,&&cov&i.knots, &&cov&i.lev);
       %end;

       %if %upcase(&&cov&i.ptype) = LAG2CUMAVGCAT %then %do;
            %makecat(&&cov&i , &&cov&i.knots , &&cov&i.lev);
            %makecat(&&cov&i.._l1, &&cov&i.knots, &&cov&i.lev);
            %makecat(&&cov&i.._l2, &&cov&i.knots, &&cov&i.lev);
            %makecat(&&cov&i.._cumavg_l1,&&cov&i.knots, &&cov&i.lev);
            %makecat(&&cov&i.._cumavg_l2,&&cov&i.knots, &&cov&i.lev);
            %makecat(&&cov&i.._cumavg_l3,&&cov&i.knots, &&cov&i.lev);
       %end;


        %if &usehistory_eof = 1 AND &&numberofeknots&i >= 3  %then %do; /* at this point covXeknots should be 0 or a list of 3 or more knots */

                 
			%if &&cov&i.etype in cumavg cumavgcat cumsum cumsumcat %then %do;
			     %if &&cov&i.etype_part2 = &timepoints %then %do;
					%let timeindex = 0;		            
				 %end;
				 %else %do;
				 	%let timeindex = %eval(&timepoints - &&cov&i.etype_part2 ) ;		            
				 %end;
			%end;
			%if &&cov&i.etype in cumavgnew cumavgcatnew cumsumnew cumsumcatnew  %then %do;
	     		%if &&cov&i.etype_part2 = &timepoints %then %do;
					%let timeindex = 0;            	 
			 	%end;
		 		%else %do;
		 			%let timeindex = %eval(&timepoints - &&cov&i.etype_part2 - 1 ) ;            	 
		 		%end;
			%end;	

		  %if (&&cov&i.etype in cumavg cumavgnew ) AND  &&numberofeknots&i >= 3   %then %do;
	            %rcspline(&&cov&i.._cumavg_&timeindex._eof , &&cov&i.eknots );	                   
	      %end;		  
		  %else %if &&cov&i.etype in cumavgcat cumavgcatnew %then %do;
		  		%makecat(&&cov&i.._cumavg_&timeindex._eof , &&cov&i.eknots , &&cov&i.elev );
		  %end;
		  %else %if (&&cov&i.etype in cumsum cumsumnew ) /*AND &usespline = 1 */   %then %do;
	             %rcspline(&&cov&i.._cumsum_&timeindex._eof , &&cov&i.eknots);	                   
	      %end;
		  %else %if &&cov&i.etype = cumsumcat or &&cov&i.etype = cumsumcatnew %then %do;
	      		%makecat(&&cov&i.._cumsum_&timeindex._eof , &&cov&i.eknots , &&cov&i.elev);        
		  %end;


		  %if &&cov&i.etype in cat skpcat spl skpspl %then %do;
		  		%if &&cov&i.etype_part2 = &timepoints %then %let  startindex = 0 ;
				%else %if &&cov&i.etype_part2 = 0 %then %let startindex = %eval(&timepoints);
				%else %let startindex = %eval(&timepoints - &&cov&i.etype_part2 ) ;

				%let stopindex = %eval(&timepoints - 1);

                  %do timeindex = &startindex %to &stopindex ;
                      %if &&cov&i.etype in spl skpspl %then %do;
							%rcspline(&&cov&i.._&timeindex._eof , &&cov&i.eknots );
					   %end;
					   %if &&cov&i.etype in cat skpcat %then %do;
							%makecat(&&cov&i.._&timeindex._eof , &&cov&i.eknots , &&cov&i.elev );
					   %end;



				  %end;
               

		  %end ;

	%end;
         
%end; 

       
      
      %if &outctype = conteofu4 %then %do;
		  z&outc = (&outc > 0) ;
		  if z&outc = 1 then l&outc = log(&outc);
	  %end ;
	  %if &outctype = cateofu %then %do;
                %do lev = 1 %to %eval(&outclev - 1);
                     outc_&lev = (&outc = &lev);
                %end;
     %end; 
      
             %interactions(outc,,1,createvar);
             %if %bquote(&compevent)^= %then %do;
                  %interactions(compevent,,1,createvar);
             %end;
             
         
             %do iii=0 %to &ncov;
                %if &&cov&iii.otype^=0 %then %do;            
                   %interactions(cov&iii ,  ,2 ,createvar );
                %end;
                %else %do;
                    %let cov&iii.ninterx=0;
                %end;
            %end;
        
       run;


    %*Dataset for basis for parameter estimates;
    data _paramdata_;
        set _inputd_; 
        keep  _sample_ newid _weight_ &outc &outcpred &compevent %if %bquote(&compevent)^=  %then %do; &compeventpred %end;  
													   &censor %if %bquote(&censor)^=  %then %do; &censorpred %end; 
													 /*  &censorcomp */
                                             &time  &wherevars  
		   %if &outctype = conteofu4 %then z&outc l&outc ;
		   %if &outctype = cateofu %then %do;
                %do lev = 1 %to %eval(&outclev - 1);
                     outc_&lev 
                %end;
           %end;
            %do i = 0 %to &ncov ;
                 &&cov&i  &&cov&i.array 
                %if &&cov&i.otype=2 %then %do;  
                       &&cov&i.._l1
                %end;    
                %if &&cov&i.otype=4 %then %do;
                       z&&cov&i l&&cov&i &&cov&i.class
                %end;
                %if &&cov&i.otype=5 %then %do;
                    %do lev = 1 %to %eval(&&cov&i.lev - 1);
                        &&cov&i.._&lev
                    %end;
                %end;
                %if &&usevisitp&i = 1 %then %do;   
                    &&cov&i.randomvisitp  
				/*	ts_last_&&cov&i.._l1 this is included in the output from listpred and is in the &&cov&i.array list */
                %end;
        
           %end; 
             ;
    run;
    
 
    %*Dataset for basis of simulation;   

     %let newsimulkeeplist = newid &time &fixedcov &keepsimuldata ;
                    %do i= 1 %to &ncov;
                        %let newsimulkeeplist = &newsimulkeeplist  &&cov&i ;                      
                        %if %length(&&cov&i.ptype) > 3  %then %do;
                                %if %substr(&&cov&i.ptype ,1,4) = lag3   %then %let newsimulkeeplist = &newsimulkeeplist &&cov&i.._l3  &&cov&i.._l2 &&cov&i.._l1 ;
                                %if %substr(&&cov&i.ptype ,1,4) = lag2   %then %let newsimulkeeplist = &newsimulkeeplist  &&cov&i.._l2 &&cov&i.._l1 ;
                                %if %substr(&&cov&i.ptype ,1,4) = lag1   %then  %let newsimulkeeplist = &newsimulkeeplist &&cov&i.._l1 ;
                                %if %substr(&&cov&i.ptype ,1,3) = skp    %then  %let newsimulkeeplist = &newsimulkeeplist &&cov&i.._l1 ;
                        %end;
                                          
                        %let newsimulkeeplist = &newsimulkeeplist  &&cov&i.else ;
                        %if &&cov&i.ptype = tsswitch1  %then %let newsimulkeeplist = &newsimulkeeplist  &&cov&i.._l1 ; 
                        %if &&usevisitp&i = 1 %then %let newsimulkeeplist = &newsimulkeeplist  ts_last_&&cov&i.._l1   &&cov&i.visitpcount  &&cov&i.randomvisitp &&cov&i.visitpelse ;  
                     %end; 


                  data _simuldata_ ;
                  set _inputd_ (where = (&time = 0) keep= newid &time &fixedcov  &keepsimuldata
                    %do i= 1 %to &ncov;
                        &&cov&i                       
                        %if %length(&&cov&i.ptype ) > 3 %then %do;
                                %if %substr(&&cov&i.ptype ,1,4) = lag3   %then  &&cov&i.._l3  &&cov&i.._l2 &&cov&i.._l1 ;
                                %if %substr(&&cov&i.ptype ,1,4) = lag2   %then  &&cov&i.._l2 &&cov&i.._l1 ;
                                %if %substr(&&cov&i.ptype ,1,4) = lag1   %then  &&cov&i.._l1 ;
                                %if %substr(&&cov&i.ptype ,1,3) = skp    %then  &&cov&i.._l1 ;
                        %end;
                                               
                        %if &&cov&i.ptype = tsswitch1  %then   &&cov&i.._l1 ;                           
                        &&cov&i.else 
                        %if &&usevisitp&i = 1 %then  ts_last_&&cov&i.._l1   &&cov&i.visitpcount  &&cov&i.randomvisitp &&cov&i.visitpelse ;  
                     %end; 
            ) ;
            run;

    
     data tmpids ;
     do newid = 1 to &ssize ;
          output;
     end;
     run;
   
    data _calchazard_ ;
    calchazard = 0 ;    
    run;

    %if &chunked = 1 AND &hazardratio = 1 AND &bootstrap_hazard = 0 %then %do;
         data &hazardname ;
         run;
    %end;
    %*Deleting no longer needed datasets;
    proc datasets library=work nolist; 
        delete cuminc indiv _orig_ ;
    quit;


%mend dataprep;

%macro samples;

    %local subseed i ;

    %************ CREATING BOOTSTRAP SAMPLES;  

    %if %eval(&bsample) = 0 %then %do;

        %*Base data for estimation and simulation (sample 0);
        
        data param;
        set _paramdata_;
        _sample_ = 0;
        run;
        
       %let subseed = %eval(2*&seed);
       * _simuldata_ contains baseline variables ;
       %if &nsimul = &ssize %then %do;
            data simul;
            set _simuldata_;
            _sample_ = 0;
            run;
        %end;
        %else %do;
       
         

            proc surveyselect data=_simuldata_ out=simul noprint
                method=urs  sampsize=&nsimul seed=&subseed;
            run;

            data simul ;
            set simul ;
             _sample_ = 0 ;
            do _copy_ = 1 to numberhits ;
                output ;
            end;
           
            drop _copy_ numberhits ;
            run;
        %end;             
                  
    %end;

    %else %do;

       /* %if &printlogstats = 1 %then */ %put  Creating bootstrap sample &bsample , seed = &seed , ssize=&ssize , nsimul = &nsimul , nparam = &nparam  ;
       /* %if &printlogstats = 1 %then %put    ; */

        %*Generating random sample of ids to be used in bootstrap sample; 

         proc surveyselect data= tmpids 
         method = urs
         n= &nparam
         seed = &seed  
         out = _idsamples (keep =  newid  numberhits  )  /* contains number of times selected into random sample (numberhits) */
         outall    
         noprint               
           ;
     
        run;

 
        %*Merging with parameter estimation data;
   
      * paramdata is original data for parameter models. paramsample is bootstrap sample with newid , numberhits, and the current bootstrap sample > 0
      this overrites the value in paramdata where newid and time agree. ;
   
        data _paramsample_ ;
        set _idsamples ;
        _sample_ = &bsample ;
        run; 

  
    * add in the variable numberhits for the number of times a subject is selected into the bootstrap sample ;
    data param ;
        merge _paramdata_ (in= p) _paramsample_;
        by newid ;
        if numberhits > 0 ; *delete those not selected into sample ;
        do _copy_ = 1 to numberhits ; * make numberhits copies of each remaining subject ;
            output ;
        end;
        drop numberhits 
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

                      

                    call symput("cov&i.min", trim(left(&&cov&i.._min))); 
                    call symput("cov&i.max", trim(left(&&cov&i.._max)));
                 end;
                
                %end;
            %end;
        %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype = conteofu4  %then %do;

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
		%else %if &outctype=conteofu4    %then %do;

             retain &outc._min &outc._max;

                if _n_ = 1 then do;
                    &outc._min =  1.0e100;
                    &outc._max = 0 ;
                    end; 

                if &outc > 0 and &outc < &outc._min then &outc._min = &outc;
                if &outc > 0 and &outc > &outc._max then &outc._max = &outc;
                
                if _end_ then do;
                    call symput("outcmin",trim(left(&outc._min)) );
                    call symput("outcmax",trim(left(&outc._max)) );
                    end;
        %end;
         _sample_ = &bsample ;
        keep  _sample_ %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype = conteofu4 %then  &outc._min &outc._max ;
                %do i = 0 %to &ncov ; 
                    %if &&cov&i.otype=3 or &&cov&i.otype=4 or &&cov&i.otype=6 or &&cov&i.otype=7   or &&cov&i.otype = -1  %then &&cov&i.._min &&cov&i.._max ;
                %end ;
                 ;
        if _end_ then output ;
    run;

    %do i = 0 %to &ncov;
        %if &&cov&i.otype=3 or &&cov&i.otype=4 or  &&cov&i.otype=6 or &&cov&i.otype=7  or &&cov&i.otype=-1 %then %do;
            %if &printlogstats = 1 %then %put  bootstrap sample &bsample &&cov&i range is &&cov&i.min to &&cov&i.max;                    
            %end;
        %end;
    %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
      %if &printlogstats = 1 %then %put  bootstrap sample &bsample &outc range is &outcmin to &outcmax;
      %end;
   
        

        %*Merging with simulation data;

        * default simul data set is to take all subjects in param data set, number of subjects in nparam ;
  
    * _simuldata_ has one observation per person for time = 0 ;
        /* idsamples containes number of times newid selected into current param data set */
        data simul ;
        merge _simuldata_  _idsamples;
        by newid; 
        _sample_ = &bsample ; 
        run;

        data simul ;
        set simul ;
        do _copy0_ = 1 to numberhits ;
                output ;
        end;
        drop numberhits _copy0_ ;
        run;

        /* want to sample when nsimul is different from size of param data set */
        %if &nsimul ne &nparam %then %do;
            * now take a random sample of size nsimul  when nsimul ne nparam ;
            %let subseed = %eval(2*&seed);
            proc surveyselect data = simul out = simul noprint
                    method=urs  sampsize=&nsimul seed=&subseed;
            run;

            data simul ;
            set simul ;
            _sample_ = &bsample ;
             do _copy_ = 1 to numberhits ;
                    output ;
             end;
             drop _copy_ numberhits ;
             run;
         %end;

 
   
     %if &hazardratio = 1  AND &bootstrap_hazard = 0  %then %do;     
         data _calchazard_ ;
         calchazard = 0 ;
         _sample_ = &bsample ;
         run;
    %end;

    proc datasets library=work nolist ;
        delete _idsamples ;
    quit;
    %end;
    
     
%mend samples;
        
%macro parameters;

    %************ ESTIMATING PARAMETERS FOR EACH SAMPLE;
    
    %local i l ;
 
    %if &printlogstats = 1 %then %put ;
    %if &printlogstats = 1 %then %put  Estimating parameters from bootstrap sample &bsample;


   sasfile param open;


   %if &outctype=binsurv   %then %do;        
    %* Binary Outcome;
        proc logistic descending data=param outest=outc ;
        &ods_logit ;
        %if %bquote(&outcwherem )^= %then where  &outcwherem    ;;
        model &outc = &outcpred  %if &uselabelo = 1 %then / parmlabel ; ;
        weight _weight_ ;
        
        data outc;
        set outc;
        if _type_='PARMS';
        _sample_ = &bsample;
        array avar intercept &outcpred;
        array abeta boutc00-boutc&dimoutc;
        do i=1 to dim(avar);
            abeta(i)=avar(i);
            end;
        keep _sample_ boutc00-boutc&dimoutc;
        run;
   %end;

  %if &outctype = bineofu %then %do;
    data param2 ;
    set param ;
    by newid ;
   
    if &time = &timepoints -1 ;  
     
    run;


    proc logistic descending data=param2 outest=outc ;
    &ods_logit ;
    %if %bquote(&outcwherem )^= %then where  &outcwherem    ;;
    model &outc = &outcpred  %if &uselabelo = 1 %then / parmlabel ; ;
    weight _weight_ ;
    run;
    
    data outc;
    set outc;
    if _type_='PARMS';
     _sample_ = &bsample;
        array avar intercept &outcpred;
        array abeta boutc00-boutc&dimoutc;
        do i=1 to dim(avar);
            abeta(i)=avar(i);
            end;
        keep _sample_ boutc00-boutc&dimoutc;
        run;



 %end;


    %if &outctype=conteofu %then %do;    
    %* Continuous outcome modeled using linear regression *** EM edited ***;

    data param2 ;
    set param ;
    by newid ;   
    if &time = &timepoints -1 ;  
    run;

    proc reg data=param2 outest=outc ;
        &ods_reg ;
         %if %bquote(&outcwherem )^= %then where  &outcwherem    ;;
        model &outc=&outcpred;
        weight _weight_ ;

    data outc;
        set outc;
        if _type_="PARMS";
        _sample_=&bsample;
        array avar intercept &outcpred;
        array abeta boutc00-boutc&dimoutc;
        do i=1 to dim(avar);
           abeta(i)=avar(i);
           end;
        se&outc=_rmse_;
        keep _sample_ boutc00-boutc&dimoutc se&outc;
    run;


    proc datasets library = work nolist ;
    delete param2 ;
    quit;

   %end;

   %else %if &outctype=conteofu2 %then %do ;
               /* Continuous outcome modeled using truncated normal using new proc qlim and bounds shifted by small ammount, similar to otype = 6   *** EM edited *** */
 
            data _null_ ;
            cmin  = &outcmin ;
            cmax =  &outcmax ;
            minoffset = 10.0 ** (log10(abs(cmin)) - 3.0) ;
            maxoffset = 10.0 ** (log10(abs(cmax)) - 3.0) ;
            call symput('outcminoffset',trim(left(minoffset))) ;
            call symput('outcmaxoffset',trim(left(maxoffset))) ;
            run;

            %if &printlogstats = 1 %then %put  ;
            %if &printlogstats = 1 %then %put  using truncated normal for  &outc : starting bounds are = [&outcmin , &outcmax] ;
           
            %let outcmin = %sysevalf(&outcmin - &outcminoffset ) ;
            %let outcmax = %sysevalf(&outcmax + &outcmaxoffset ) ;
            %if &printlogstats = 1 %then %put  enlarge bounds by .1 percent of previous bounds ;
            %if &printlogstats = 1 %then %put  new bounds will be [&outcmin , &outcmax] using offset of &outcminoffset , &outcmaxoffset  ;
            %if &printlogstats = 1 %then %put  ;


            data _covbounds_ ;
            set _covbounds_ ;
            &outc._min = &outcmin ;
            &outc._max = &outcmax ;
            run;


            data param2 ;
            set param ;
            by newid ;            
            if &time = &timepoints -1 ;               
            run;

            proc qlim   data=param2(keep =  &outc &compevent  &time &outcpred _weight_ )
                outest=outc  ;
                &ods_qlim ;   
                %if %bquote(&outcwherem )^= %then where  &outcwherem    ;; 
                model  &outc=&outcpred /  truncated (lb=  %sysfunc(compress(&outcmin))  ub= %sysfunc(compress(&outcmax)) );
                weight _weight_ ;
            run; ;
            
            data outc;
            set outc; 
            if _type_="PARM";
            _sample_=&bsample;
            array avar intercept &outcpred;
            array abeta boutc00-boutc&dimoutc;
            do i=1 to dim(avar);
               abeta(i)=avar(i);
            end;
            se&outc=_Sigma;
            keep _sample_ boutc00-boutc&dimoutc se&outc;
            run;
          
            proc datasets library = work nolist ;
            delete param2 ;
            quit;
         %end ;
      
          %else %if  &outctype=conteofu3 %then %do ;
               /*  Continuous outcome modeled using tobit density, similar to what is being done using otype = 7 *** EM edited ***  */
        

            data param2 ;
            set param ;
            by newid ;           
            if &time = &timepoints -1 ;             
            run;

              proc qlim   data=param2(keep =  &outc &compevent   &time &outc &outcpred _weight_ )
                outest=outc  ;
                &ods_qlim ;    
                %if %bquote(&outcwherem )^= %then where  &outcwherem    ;; 
                model  &outc=&outcpred  /  censored (lb=  %sysfunc(compress(&outcmin))  ub= %sysfunc(compress(&outcmax)) );
                weight _weight_ ;
            run; 

           data outc;
           set outc;
           if _type_="PARM";
        _sample_=&bsample;
        array avar intercept &outcpred;
        array abeta boutc00-boutc&dimoutc;
        do i=1 to dim(avar);
           abeta(i)=avar(i);
           end;
        se&outc=_Sigma;
        keep _sample_ boutc00-boutc&dimoutc se&outc;
    run;

    proc datasets library = work nolist ;
    delete param2 ;
    quit;
        
         %end ;

		%else %if &outctype=conteofu4  %then %do;
        /*  Continuous outcome modeled using logistic to log-linear approach, similar to what is being done using otype = 4 *** EM edited ***  */
             
			    data param2 ;
			    set param ;
			    by newid ;   
			    if &time = &timepoints -1 ;  
			    run;

                proc logistic descending 
                    data=param2(keep =  _weight_  &compevent   &time z&outc   &outcpred )
                    outest=z&outc ;
                    &ods_logit ;
                    %if %bquote(&outcwherem )^= %then where  &outcwherem    ;; 
                    model z&outc = &outcpred %if &uselabelc = 1 %then / parmlabel ;;
                    weight _weight_ ;
                run;
            
        
            proc reg
                data=param2(keep =  _weight_  &compevent   &time z&outc l&outc  &outcpred )
                outest=&outc ;
                &ods_reg ;
                %if %bquote(&outcwherem )^= %then where  &outcwherem    ;;  
                model l&outc = &outcpred;
                weight _weight_ ;
            run;


			data z&outc;
            set z&outc;
            if _type_='PARMS';
            _sample_ = &bsample;
            array avar intercept &outcpred;
            array abeta boutcz_00-boutcz_&dimoutc;
            do j=1 to dim(avar); 
				abeta(j)=avar(j); 
            end;
            keep _sample_ boutcz_00-boutcz_&dimoutc ;
            run;
                
            data &outc;
            set &outc;
            if _type_='PARMS';
            _sample_ = &bsample;
            array avar intercept &outcpred;
            array abeta boutc00-boutc&dimoutc;
            do j=1 to dim(avar); 
				abeta(j)=avar(j); 
			end;
            se&outc=_rmse_;
            keep _sample_ boutc00-boutc&dimoutc;
            run;

            data outc;                 
            merge &outc z&outc;                    
            by _sample_;
            run;

			proc datasets library = work nolist ;
    		delete param2 z&outc &outc ;
    		quit;
        %end;

		%else %if &outctype = cateofu %then %do;
        /*  Categorical or ordinal outcomes modeled using nested logistic regression models *** EM edited *** */
        

		 data param2 ;
            set param ;
            by newid ;   
            if &time = &timepoints -1 ;  
         run;

			%do lev = 1 %to %eval(&outclev-1);
                proc logistic descending data=param2(keep =  _weight_   &compevent   &time   &outcpred
                    outc_&lev 
                    %do lowerlev = 1 %to %eval(&lev-1);
                        outc_&lowerlev
                     %end;
                        )
                    outest=outc_&lev ;
                        &ods_logit ;
                        where  &time > 0  and %unquote(&outcwherem)  
                            %do lowerlev = 1 %to %eval(&lev-1);
                                and outc_&lowerlev ne 1
                            %end;
                               
                                ;
                        model outc_&lev = &outcpred %if &uselabelc = 1 %then / parmlabel ;;
                        weight _weight_ ;
                run;
            %end;

 			%do lev = 1 %to %eval(&outclev - 1);
                data outc_&lev;
                set outc_&lev ;
                if _type_='PARMS';
                _sample_ = &bsample;
                    array avar_&lev intercept &outcpred;
                    array about_&lev about_&lev._00-about_&lev._&dimoutc;
                    do j=1 to dim(avar_&lev); 
                         about_&lev(j)=avar_&lev.(j);
                    end;
                        keep _sample_ about_&lev._00-about_&lev._&dimoutc  ;
                run;
            %end;
            data outc;
            merge %do lev = 1 %to %eval(&outclev - 1); outc_&lev %end; ;
            by _sample_;
            run;

			proc datasets library = work nolist ;
    		delete param2 %do lev = 1 %to %eval(&outclev - 1); outc_&lev %end; ;
    		quit;

		%end;


    %* Censoring;    
    %if %bquote(&compevent)^= AND &compevent_cens = 0  %then %do;
        proc logistic descending data=param outest=compevent ;
            &ods_logit ;
             %if %bquote(&compeventwherem )^= %then where  &compeventwherem    ;;
            model &compevent = &compeventpred %if &uselabelo = 1 %then / parmlabel ;;		
            run;
        data compevent;
            set compevent;
            if _type_='PARMS';
            _sample_ = &bsample;
            array avar intercept &compeventpred;
            array abeta bcompevent00-bcompevent&dimcompevent;
            do i=1 to dim(avar);
                abeta(i)=avar(i);
                end;
            keep _sample_ bcompevent00-bcompevent&dimcompevent;
        run;

    %end;


	 %if %bquote(&censor)^=  AND &check_cov_models = 1  %then %do;

		%local covlisttmp ;
		%do i = 1 %to &ncov ;
			%let covlisttmp = &covlisttmp &&cov&i ;
			%if &&usevisitp&i = 1 %then %let covlisttmp = &covlisttmp &&cov&i.randomvisitp ;
		 %end;


        proc logistic descending data=param  ;
            &ods_logit ;
             %if %bquote(&censorwherem )^= %then where  &censorwherem    ;;
            model &censor = &censorpred %if &uselabelo = 1 %then / parmlabel ;;
			output out=_censor_ (keep = newid &time &censor &compevent  &outc pC_d  &covlisttmp ) pred = pC_d ;
            run;
			%if &compevent_cens = 1  %then %do;

			 	proc logistic descending data=param  ;
            	&ods_logit ;
             	%if %bquote(&censorwherem )^= %then where  &censorwherem    ;;
            	model &compevent = &censorpred %if &uselabelo = 1 %then / parmlabel ;;
				output out=_&compevent._ (keep = newid &time &compevent &outc pCE_d ) pred = pCE_d ;
            	run;
			%end;

			data _censor_ ;
			set _censor_ ;
			by newid ;
			retain _wtC_ ;
			if first.newid then _wtC_ = 1 ;
			_wtCl1_ = _wtC_ ; * the retain statement is the lagged value before it is updated below ;
			if &censor = 0 then _wtC_ = _wtC_ * (1/(1-pC_d));
			else _wtC_ = _wtC_ * 0; /* censor = 1 or missing */
			/* when there is no compevent or there is a compevent but it is NOT a second censoring variable **/
		     %if  &compevent_cens = 0  %then  rename  _wtC_ = _wt_  _wtCl1_ = _wtl1_ ;;
			run;

			%if &compevent_cens = 1 %then %do;
		    	data _&compevent._ ;
				set _&compevent._ ;
				by newid ;
				retain _wtCE_ ;
				if first.newid then _wtCE_ = 1 ;
				_wtCEl1_ = _wtCE_;
				
				if &compevent = 0 then _wtCE_ = _wtCE_ * (1/(1-pCE_d));
				else _wtCE_ = _wtCE_ * 0; /* compevent = 1 or missing */
				keep  newid &time &compevent _wtCE_ _wtCEl1_ ; 
				run;
	
				data _censor_ ;
        	    merge _censor_ _&compevent._ ;
				by newid &time ;
				_wt_ = _wtC_ * _wtCE_ ;
				_wtl1_ = _wtCl1_ * _wtCEl1_ ;
				keep &time _wt_ _wtl1_ &covlisttmp &outc &compevent ;
				run;
			%end;

			proc univariate data = _censor_ ;
			var _wt_ ;
			run;

			%if %upcase(%substr(&maxipw ,1,1)) = P %then %do; 

				proc means data = _censor_ &maxipw ;
				var _wt_ _wtl1_;
				output out = _p99_ (keep = _maxipw_ _maxipwl1_) &maxipw (_wt_ _wtl1_ ) = _maxipw_ _maxipwl1_   ;
				run;

				data _null_;
				set _p99_ ;
				call symput('maxipw',trim(left(_maxipw_)));
				call symput('maxipwl1',trim(left(_maxipwl1_)));
				run;
			%end;
			%else %do;
				%let maxipwl1 = &maxipw ;
			%end;

			data _censor_;
			set _censor_ ;
			if _wt_ > &maxipw then _wt_ = &maxipw;
			if _wtl1_ > &maxipwl1 then _wtl1_ = &maxipwl1 ;
			run;



            %if &outctype ne cateofu %then %do;
				proc means data = _censor_   ;
				class &time ;
				var &outc  %if %bquote(&compevent) ^= AND &compevent_cens = 0 %then &compevent ;;
				types &time ;
				weight _wt_ ;
				%if &outctype=binsurv %then %do ;
					output out = forwtY (keep = &time meanoutc %if %bquote(&compevent) ^= AND &compevent_cens = 0 %then meancompevent ;) 
	                        mean( &outc %if %bquote(&compevent) ^= AND &compevent_cens = 0 %then &compevent ;) = 
	       				meanoutc %if %bquote(&compevent) ^= AND &compevent_cens = 0 %then meancompevent ; ;
				%end;
				%else %do ;
					output out = forwtY (keep = &time &outc )   mean( &outc ) = &outc  ;
				%end;
				run;
			%end;
			%else %do;
				proc freq data = _censor_ (where = (&time = %eval(&timepoints - 1))) noprint ;
				tables &outc / out= proportions (keep = &outc percent rename = (percent = proportion) ) nocum ;
				weight _wt_ ;
				run;

				
				data proportions ;
				set proportions ;
				label &outc="&outc level" proportion="Proportion (%)";
				proportion = round(proportion,0.01);
				sample = &bsample ;
				run;
			 
			%end;
			*calculate E[L_k * (1-Y_k) &W_{k-1} ] / E[(1-Y_k)*W_{k-1}]
			 If there is a L_k then it must be that Y_k = 0 , Same as E[L_k * W_{k-1}]/E[W_{k-1}] ;
			proc means data = _censor_    ;
			class &time ;
			var &covlisttmp ;
			types &time ;
			weight _wtl1_ ;
			%if &outctype=binsurv %then %do ;
				output out = forwtCov (keep = &time &covlisttmp) mean( &covlisttmp) = &covlisttmp ;
			%end;
			%else %do ;
				output out = forwtCov (keep = &time &covlisttmp) mean( &covlisttmp) = &covlisttmp ;
			%end;
			run;


			%if &bsample = 0 %then %do;
				%if &outctype = binsurv  %then %do;  
	        		data cuminc;
	        		set forwtY;
	        		keep cuminc;
	        		by &time;
	        		retain cumsurv 1 cuminc 0;   
					%if (%bquote(&compevent) = ) OR ( &compevent_cens = 1 ) %then meancompevent = 0 ;;
	        		inc = cumsurv * meanoutc * (1.0 - meancompevent) ;         
	        		cuminc = cuminc + inc;
	        		surv = (1.0 - meanoutc) * (1.0 - meancompevent) ;         
	        		cumsurv = cumsurv * surv;
	        		if _N_ = &timepoints then call symput('obsp',trim(left(cuminc))); /* this is based on row number and not time value so there is no -1 here */
	        		run;
			   %end ;
			   %else %do;
			       %if &outctype ne cateofu %then %do;
					data cuminc ;
					set forwtY ;
					if _N_ = &timepoints then call symput('obsp',trim(left(&outc)));
					run;

				  %end;
				  %else %do ;
				  	data proportions0 ;
					set proportions ;
				    run;
				  %end;
			   %end;
          %end;
    %end;
    
     
        
    %* Looping over covariates for regressions;

    %do i = 0 %to &ncov;
        %if &outctype=conteofu %then %do;
            %let outcond=%nrstr((&outc ne .));
        %end;
       

       %if &&usevisitp&i = 1 %then %do;


            proc logistic descending
                data=param(keep =  _weight_ &outc &compevent   &time &&cov&i.randomvisitp &&cov&i.array 
                &wherevars  /* ts_last_&&cov&i.._l1  this should be included in &&cov&i.array */ )
                outest=&&cov&i.randomvisitp
                %if (&outputs ^= yes or %eval(&bsample) ^= 0) %then noprint ; ;
            where &time > 0  and &time not in ( &&cov&i.skip) and  %unquote(&&cov&i.visitpwherem)  and ts_last_&&cov&i.._l1 ne &&cov&i.visitpmaxgap  ;;
            model &&cov&i.randomvisitp = &&cov&i.array %if &uselabelc = 1 %then / parmlabel ;;
            weight _weight_ ;
            run;



            data &&cov&i.randomvisitp;
            set &&cov&i.randomvisitp;
            if _type_='PARMS';
            _sample_ = &bsample;
            array avar intercept &&cov&i.array;          
            array abeta bvar&i.visitp_00-bvar&i.visitp_&&dimvar&i;
            do j=1 to dim(avar); 
                abeta(j)=avar(j); 
            end;
            keep _sample_ bvar&i.visitp_00-bvar&i.visitp_&&dimvar&i;          
            run;

       %end;
        %if &&cov&i.otype=1 %then %do;
            proc logistic descending
                data=param(keep = _weight_ &outc &compevent   &time &&cov&i  &&cov&i.randomvisitp  &&cov&i.array &wherevars )
                outest=&&cov&i ;
                &ods_logit ;
                where &time > 0  and %unquote(&&cov&i.wherem) and &time not in (&&cov&i.skip) %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ;  ;
                model &&cov&i = &&cov&i.array %if &uselabelc = 1 %then / parmlabel ;;
                weight _weight_ ;
            run;
        %end;
        %else %if &&cov&i.otype=2 %then %do;
            proc logistic descending
                data=param(keep = _weight_ &outc &compevent   &time &&cov&i  &&cov&i.randomvisitp &&cov&i.array &&cov&i.._l1 &wherevars )
                outest=&&cov&i ;
                &ods_logit;
                where &time > 0  and %unquote(&&cov&i.wherem) and &&cov&i.._l1=0 and &time not in (&&cov&i.skip)   %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ; ;
                model &&cov&i = &&cov&i.array %if &uselabelc = 1 %then / parmlabel ;;
                weight _weight_ ;
            run;
        %end;

        %else %if &&cov&i.otype=3 %then %do;
            /* continuous out type will generate values within recorded range : attempt at truncated normal */
            proc reg
                data=param(keep =  _weight_ &outc &compevent   &time &&cov&i  &&cov&i.randomvisitp &&cov&i.array &wherevars )
                outest=&&cov&i ;
                &ods_reg ;
                where &time > 0  and %unquote(&&cov&i.wherem) and &time not in (&&cov&i.skip)  %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ; ;
                model &&cov&i = &&cov&i.array;
                weight _weight_ ;
            run;
            %end; 

        %else %if &&cov&i.otype=4 %then %do;
            %if %bquote(&&cov&i.class)^= %then %do;
                proc logistic descending 
                    data=param(keep =  _weight_ &outc &compevent   &time z&&cov&i   &&cov&i.randomvisitp &&cov&i.array &&cov&i.class &wherevars )
                    outest=c0_z&&cov&i ;
                    &ods_logit ;
                    where &time > 0  and %unquote(&&cov&i.wherem) and
                        &time not in (&&cov&i.skip)
                        and &&cov&i.class = 0   %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ; ;
                    model z&&cov&i = &&cov&i.array %if &uselabelc = 1 %then / parmlabel ;;
                    weight _weight_ ;
                run;
                proc logistic descending 
                    data=param(keep =  _weight_ &outc &compevent   &time z&&cov&i  &&cov&i.randomvisitp &&cov&i.array &&cov&i.class &wherevars )
                    outest=c1_z&&cov&i ;
                    &ods_logit ;
                    where &time > 0  and %unquote(&&cov&i.wherem) and
                        &time not in (&&cov&i.skip)
                        and &&cov&i.class ^= 0   %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ; ;
                    model z&&cov&i = &&cov&i.array  %if &uselabelc = 1 %then / parmlabel ; ;
                    weight _weight_ ;
                run;
             %end;
             %else %do;
                proc logistic descending 
                    data=param(keep =  _weight_ &outc &compevent   &time z&&cov&i  &&cov&i.randomvisitp &&cov&i.array &wherevars )
                    outest=z&&cov&i ;
                    &ods_logit ;
                    where &time > 0  and %unquote(&&cov&i.wherem) and &time not in (&&cov&i.skip)   %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ; ;
                    model z&&cov&i = &&cov&i.array %if &uselabelc = 1 %then / parmlabel ;;
                    weight _weight_ ;
                run;
             %end;
        
            proc reg
                data=param(keep =  _weight_ &outc &compevent   &time z&&cov&i l&&cov&i  &&cov&i.randomvisitp  &&cov&i.array &wherevars )
                outest=&&cov&i ;
                &ods_reg ;
                where &time > 0  and %unquote(&&cov&i.wherem) and z&&cov&i=1 and &time not in (&&cov&i.skip)  %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ; ;
                model l&&cov&i = &&cov&i.array;
                weight _weight_ ;
            run;
        %end;

        %else %if &&cov&i.otype=5 %then %do;
            %do lev = 1 %to %eval(&&cov&i.lev-1);
                proc logistic descending data=param(keep =  _weight_ &outc &compevent   &time &wherevars   &&cov&i.randomvisitp
                    &&cov&i.._&lev &&cov&i.array
                    %do lowerlev = 1 %to %eval(&lev-1);
                        &&cov&i.._&lowerlev
                        %end;
                        )
                    outest=&&cov&i.._&lev ;
                        &ods_logit ;
                        where  &time > 0  and %unquote(&&cov&i.wherem) and &time not in (&&cov&i.skip)
                            %do lowerlev = 1 %to %eval(&lev-1);
                                and &&cov&i.._&lowerlev ne 1
                                %end;
                               %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ; 
                                ;
                        model &&cov&i.._&lev = &&cov&i.array %if &uselabelc = 1 %then / parmlabel ;;
                        weight _weight_ ;
                run;
                %end;

            %end;

         %else %if  &&cov&i.otype = 6 %then %do ;
               /* true truncated normal using new proc qlim and bounds shifted by small ammount  */
 
            data _null_ ;
            cmin  = &&cov&i.min ;
            cmax = &&cov&i.max ;
            if abs(cmin) > 0 then minoffset = 10.0 ** (log10(abs(cmin)) - 2.0) ;                                    
            if abs(cmax) > 0 then maxoffset = 10.0 ** (log10(abs(cmax)) - 2.0) ;

            if cmin = 0 and abs(cmax) > 0 then minoffset = maxoffset ;
            if cmax = 0 and abs(cmin) > 0 then maxoffset = minoffset ;

            call symput('minoffset',trim(left(minoffset))) ;
            call symput('maxoffset',trim(left(maxoffset))) ;
            run;

            %if &printlogstats = 1 %then %put  ;
            %if &printlogstats = 1 %then %put  using truncated normal for &&cov&i : starting bounds are = [&&cov&i.min , &&cov&i.max] ;
           
            %let cov&i.min = %sysevalf(&&cov&i.min - &minoffset ) ;
            %let cov&i.max = %sysevalf(&&cov&i.max + &maxoffset ) ;
            %if &printlogstats = 1 %then %put  enlarge bounds by 1 percent of previous bounds ;
            %if &printlogstats = 1 %then %put  new bounds will be [&&cov&i.min , &&cov&i.max] using offset of &minoffset , &maxoffset  ;
            %if &printlogstats = 1 %then %put  ;

           data _covbounds_;
           set _covbounds_ ;
           &&cov&i.._min = &&cov&i.min ;
           &&cov&i.._max = &&cov&i.max ;
           run;
            proc qlim   data=param(keep =  _weight_ &outc &compevent   &time &&cov&i  &&cov&i.randomvisitp &&cov&i.array &wherevars )
                outest=&&cov&i  ;
                &ods_qlim ; 
                where &time > 0  and %unquote(&&cov&i.wherem) and &time not in (&&cov&i.skip)  %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ; ;
                model &&cov&i = &&cov&i.array  /  truncated (lb=  %sysfunc(compress(&&cov&i.min))  ub= %sysfunc(compress(&&cov&i.max)) );
                weight _weight_ ;
            run; ;
            

         %end ;
      
          %else %if  &&cov&i.otype = 7 %then %do ;
               /*  tobit density, similar to what is being done using otype = 3   */
        
              proc qlim   data=param(keep =  _weight_ &outc &compevent   &time &&cov&i  &&cov&i.randomvisitp &&cov&i.array &wherevars )
                outest=&&cov&i  ;
                &ods_qlim ;
                where &time > 0  and %unquote(&&cov&i.wherem) and &time not in (&&cov&i.skip)  %if &&usevisitp&i = 1 %then and &&cov&i.randomvisitp = 1 ; ;
                model &&cov&i = &&cov&i.array  /  censored (lb=  %sysfunc(compress(&&cov&i.min))  ub= %sysfunc(compress(&&cov&i.max)) );
                weight _weight_ ;
            run; 
            
        
         %end ;
          
         %else %if &&cov&i.otype = -1 %then %do; 
          /* user defined otype */
                  %&&cov&i.modusermacro ;

         %end;
    

    %end; /* end covariate loop. */
       

    sasfile param close;

    
    %*Looping over covariates to create the data sets;
    %do i = 0 %to &ncov;
        %if &&cov&i.otype=1 or &&cov&i.otype=2 or &&cov&i.otype=3 %then %do;
            data &&cov&i;
                set &&cov&i;
                if _type_='PARMS';
                _sample_ = &bsample;
                array avar intercept &&cov&i.array;       
                array abeta bvar&i._00-bvar&i._&&dimvar&i;
                do j=1 to dim(avar); abeta(j)=avar(j); end;
                    %if &&cov&i.otype=3 %then %do;
                        se&&cov&i=_rmse_;
                        %end;
                    keep _sample_ bvar&i._00-bvar&i._&&dimvar&i 
                        %if &&cov&i.otype=3 %then %do;
                        se&&cov&i
                            %end;
                        ;
            run;

            %if &&usevisitp&i = 1 %then %do;            

               data &&cov&i;
               merge &&cov&i.randomvisitp &&cov&i;
               by _sample_;
               run;

            %end; 
    
        %end;


        %else %if &&cov&i.otype=4 %then %do;
            %if %bquote(&&cov&i.class)^= %then %do;
                data c0_z&&cov&i;
                    set c0_z&&cov&i;
                    if _type_='PARMS';
                    _sample_ = &bsample;
                    array avar intercept &&cov&i.array;
                    array abeta bvar&i.z0_00-bvar&i.z0_&&dimvar&i;
                    do j=1 to dim(avar); abeta(j)=avar(j); end;
                        keep _sample_ bvar&i.z0_00-bvar&i.z0_&&dimvar&i;
                run;

                data c1_z&&cov&i;
                    set c1_z&&cov&i;
                    if _type_='PARMS';
                    _sample_ = &bsample;
                    array avar intercept &&cov&i.array;
                    array abeta bvar&i.z1_00-bvar&i.z1_&&dimvar&i;
                    do j=1 to dim(avar); abeta(j)=avar(j); end;
                        keep _sample_ bvar&i.z1_00-bvar&i.z1_&&dimvar&i;
                run;

                %end;
            %else %do;
                data z&&cov&i;
                    set z&&cov&i;
                    if _type_='PARMS';
                    _sample_ = &bsample;
                    array avar intercept &&cov&i.array;
                    array abeta bvar&i.z_00-bvar&i.z_&&dimvar&i;
                    do j=1 to dim(avar); abeta(j)=avar(j); end;
                        keep _sample_ bvar&i.z_00-bvar&i.z_&&dimvar&i;
                run;
                %end;
            data &&cov&i;
                set &&cov&i;
                if _type_='PARMS';
                _sample_ = &bsample;
                array avar intercept &&cov&i.array;
                array abeta bvar&i._00-bvar&i._&&dimvar&i;
                do j=1 to dim(avar); abeta(j)=avar(j); end;
                    se&&cov&i=_rmse_;
                    keep _sample_ bvar&i._00-bvar&i._&&dimvar&i se&&cov&i;
            run;
            data &&cov&i;
                %if %bquote(&&cov&i.class)^= %then %do;
                    merge &&cov&i c0_z&&cov&i c1_z&&cov&i;
                    %end;
                %else %do;
                    merge &&cov&i z&&cov&i;
                    %end;
                by _sample_;
            run;

          
           %if &&usevisitp&i = 1 %then %do;            

               data &&cov&i;
               merge &&cov&i.randomvisitp &&cov&i;
               by _sample_;
               run;

            %end;      
            
        %end;

        %else %if &&cov&i.otype=5 %then %do;
            %do l = 1 %to %eval(&&cov&i.lev - 1);
                data &&cov&i.._&l;
                    set &&cov&i.._&l;
                    if _type_='PARMS';
                    _sample_ = &bsample;
                    array avar_&l intercept &&cov&i.array;
                    array abeta_&l bvar&i._&l._00-bvar&i._&l._&&dimvar&i;
                    do j=1 to dim(avar_&l); abeta_&l.(j)=avar_&l.(j); end;
                        keep _sample_  bvar&i._&l._00-bvar&i._&l._&&dimvar&i;
                run;
                %end;
            data &&cov&i;
                merge %do l = 1 %to %eval(&&cov&i.lev - 1); &&cov&i.._&l %end; ;
                by _sample_;
            run;



            %if &&usevisitp&i = 1 %then %do;            

               data &&cov&i;
               merge &&cov&i.randomvisitp &&cov&i;
               by _sample_;
               run;

            %end; 

        %end;
        %else %if &&cov&i.otype=6 or &&cov&i.otype=7   %then %do;
            data &&cov&i;
                set &&cov&i;
                if _type_='PARM';
                _sample_ = &bsample;
                array avar intercept &&cov&i.array;       
                array abeta bvar&i._00-bvar&i._&&dimvar&i;
                do j=1 to dim(avar); abeta(j)=avar(j); end;
                    
                        se&&cov&i=_Sigma;
                        
                    keep _sample_ bvar&i._00-bvar&i._&&dimvar&i 
                        
                        se&&cov&i
                           
                        ;
            run;


           %if &&usevisitp&i = 1 %then %do;            

               data &&cov&i;
               merge &&cov&i.randomvisitp &&cov&i;
               by _sample_;
               run;

            %end; 


       %end;

      
       
       %else %if &&cov&i.otype=-1 %then %do; 

                 %&&cov&i.moddatausermacro ;


        %end;
    



        %end; /*end of i loop*/



    %* Merging all the parameter estimates into _beta_ and then into &betadata;


    data samplebetas;
        merge outc
            %if %bquote(&compevent)^= AND &compevent_cens = 0 %then %do; compevent %end;
             
                %do i = 0 %to &ncov;
                    %if &&cov&i.otype ne 0 %then %do;
                        &&cov&i
                        %end;
                    %end;
                    ;
          by _sample_;
                _sample_ = &bsample;
    run;

    proc append base=_beta_ data=samplebetas;
    run;
       
    proc sort data=_beta_; by _sample_;
    run;
    

    %if &bsample =   &sample_end   %then %do;            
         
        
        %if %bquote(&betadata)^= %then %do;
            %if &printlogstats = 1 %then %put ;
            %if &printlogstats = 1 %then %put  Outputting parameter estimates to &betadata;
            data &betadata;
                set _beta_;
            run;
            %end;
        %if &printlogstats = 1 %then %put ;
        
        %end;


/******** move code for covbetas to end of submacro *****/


     %if &check_cov_models = 1 OR &rungraphs = 1   %then %do;
          %local dataholder ;
		  %if &bsample = &sample_start %then %let dataholder = &covmeanname ;  
          %else %let dataholder = _cov_mean_tmp ; 

          proc sql ;
            create table &dataholder  as
            select &bsample as _sample_ , 
                &time as _time_ , 
                %if &outctype = binsurv %then %do;
                      max(sum(&outc),0) as e,
                      mean(&outc) as meanoutc ,
                      %if %bquote(&compevent)^= AND &compevent_cens = 0 %then max(sum(&compevent),0) as r , ;
                      %else 0 as r , ;
                      %if %bquote(&compevent)^= AND &compevent_cens = 0  %then mean(&compevent) as meancompevent , ;
                      %else 0 as meancompevent , ;                       
                      count(&time) as n ,
                %end;
                %else %do;
                     mean(&outc * _weight_ ) as &outc , 
                %end;
                %do i = 0 %to %eval( &ncov - 1) ;
                     mean(&&cov&i %if &&cov&i.otype > 0 %then  * _weight_ ; ) as &&cov&i ,
                     %if &&usevisitp&i = 1 %then mean(&&cov&i.randomvisitp * _weight_ ) as &&cov&i.randomvisitp  ,; 
                %end;
                %if &&usevisitp&ncov = 1 %then mean(&&cov&ncov.randomvisitp * _weight_ ) as &&cov&ncov.randomvisitp  ,; 
                mean(&&cov&ncov %if &&cov&i.otype > 0 %then * _weight_ ; ) as &&cov&ncov 
              
           from param(keep = &time  _weight_  %do i = 0 %to &ncov;  &&cov&i &&cov&i.randomvisitp  %end; &outc &compevent   )
        
           group by &time 
           order by &time 
           ;

      
      quit;

	  %if %bquote(&censor) ^= %then %do;
			data  &dataholder ;
			%if &outctype = binsurv %then %do;
 				merge  &dataholder 
                       forwtY (keep = &time meanoutc %if %bquote(&compevent) ^= AND &compevent_cens = 0 %then meancompevent ; ) 
                       forwtCov (keep = &time &covlisttmp ) ;
			%end;
			%else %do;
				merge  &dataholder 
                     %if &outctype ne cateofu %then  forwtY (keep = &time &outc  )   ;
                       forwtCov (keep = &time &covlisttmp )  ;
			%end;
			by &time ;
            run;		

			proc datasets library = work nolist ;
			delete _censor_ forwtY forwtCov;
			quit;
	  %end;
        
        %if &bsample > &sample_start  %then %do;
           proc append base = &covmeanname data = _cov_mean_tmp  ;
           run;

           proc datasets library = work nolist ;
                delete _cov_mean_tmp;
          quit;

        %end;  
        %if &chunked = 1 %then %do;
            
            proc copy in = work out = &savelib   ;
                select &covmeanname  ;
            run;
        %end;
        
     %end;


/********************************************************/





    proc datasets library=work nolist; 
        delete samplebetas outc %if &bsample > 0 %then  _paramsample_ ; 
            %if %bquote(&compevent)^= AND &compevent_cens = 0 %then %do; compevent %end;
            
                %do i = 0 %to &ncov;
                    %if &&cov&i.otype ne 0 %then %do; &&cov&i %end; 
                    %if &&cov&i.otype=4 %then %do; 
                        %if %bquote(&&cov&i.class)^= %then %do;
                            c0_z&&cov&i c1_z&&cov&i 
                        %end;
                        %else %do;
                            z&&cov&i 
                        %end;           
                    %end;
                    %if &&cov&i.otype=5 %then %do; 
                        %do lev = 1 %to %eval(&&cov&i.lev-1); &&cov&i.._&lev %end; 
                    %end;          
               %end;
                ;
    quit;



%mend parameters;

%macro mynull ;
/* do nothing here, just a place holder */
%mend;


%macro interv (intno= , intlabel=, nintvar= 0, intvisittype=1, intcond = 1=1, intsetup= mynull ,
          intvar1= ,
          inttype1 =, 
                inttimes1 = (-1),
          intpr1 = 1,
                intcov1 = ,
                intmean1 = ,
          intsd1 = ,
                intvalue1 = ,
          intmin1 =. ,
          intmax1 =. ,
          intchg1 = ,          
          intusermacro1= , 
        

          intvar2 = , inttype2 = , inttimes2 = (-1), intpr2 = 1, intcov2 = ,
          intmean2 = , intsd2 = , intvalue2 = , intmin2 =. , intmax2 =. , intchg2 = ,
          intusermacro2=,
         
          intvar3 = , inttype3 = , inttimes3 = (-1), intpr3 = 1, intcov3 = ,
                intmean3 = , intsd3 = , intvalue3 = , intmin3 =. , intmax3 =. , intchg3 = ,
          intusermacro3=,
          
                intvar4 = , inttype4 = , inttimes4 = (-1), intpr4 = 1, intcov4 = ,
                intmean4 = , intsd4 = , intvalue4 = , intmin4 =. , intmax4 =. , intchg4 = ,
          intusermacro4=,
        
          intvar5 = , inttype5 = , inttimes5 = (-1), intpr5 = 1, intcov5 = ,
                intmean5 = , intsd5 = , intvalue5 = , intmin5 =. , intmax5 =. , intchg5 = ,
          intusermacro5=,
         
          intvar6 = , inttype6 = , inttimes6 = (-1), intpr6 = 1, intcov6 = ,
                intmean6 = , intsd6 = , intvalue6 = , intmin6 =. , intmax6 =. , intchg6 = ,
          intusermacro6=,
         
          intvar7 = , inttype7 = , inttimes7 = (-1), intpr7 = 1, intcov7 = ,
                intmean7 = , intsd7 = , intvalue7 = , intmin7 =. , intmax7 =. , intchg7 = ,
          intusermacro7=,
         
          intvar8 = , inttype8 = , inttimes8 = (-1), intpr8 = 1, intcov8 = ,
                intmean8 = , intsd8 = , intvalue8 = , intmin8 =. , intmax8 =. , intchg8 = ,
      intusermacro8 = ,
      

      timesnotelig = -1  

    );

    %************ SIMULATION TO GET CUMULATIVE INCIDENCE UNDER INTERVENTION;

    %local n j i ;
    %if &printlogstats = 1 %then %put  Computing intervention &intno;

    %*Sampling observed distribution for intervention type 5 and 6;    
    %do n = 1 %to &nintvar;
        %if &&inttype&n = 5 or &&inttype&n = 6 %then %do;
 

            data distrib;
                set param;
                
                %if %bquote(&&intmax&n)^=. %then %do;
                    if &&intvar&n <= &&intmax&n;
                %end;
                %if %bquote(&&intmin&n)^= . %then %do;
                    if &&intvar&n >= &&intmin&n;
                %end;

                r_&&intvar&n.._&n = &&intvar&n; /* need to allow for the same variable to appear more than once in an intervention */
            run;

            
            data distrib;
                set distrib end=_end_;
                keep distid  r_&&intvar&n.._&n;
                distid = _N_;                    
                if _end_ then do;
                    call symput('dsize',trim(left(put(distid,8.))) );
                    end;
            run;
            
            proc sort data=distrib;
                by distid;
            run;

            data simul;
                set simul  ;
                call streaminit(%eval(4*&seed));
                U = rand('uniform');
                distid=int(U*&dsize)+1;
            run;

            proc sort data=simul;
                by distid;
            run;

 


            data simul;
                merge simul(in=s ) distrib;
                by distid;
                if s;
                drop U distid;
            run;

            proc datasets nolist library=work; 
                delete distrib;
            quit;
            
        %end;
   %end;            
    
    %*Sort data;
    proc sort data = simul;
        by _sample_;
    run;

  
    %*Outputting the mean of covariates and probability of event;
	%if &outctype = cateofu %then %do;

		 %put BSAMPLE = &bsample  ;

		proc sql ;
		create table mytest as 
		select &outc as &outc._level , count(*) as freq ,  sum (intervened) as intervened , sum(averinterv) as  averinterv  
		                    
        %if &minimalistic = no %then %do;
            %do i = 1 %to &ncov;
               %do j = 1 %to %eval(&timepoints );
                  %if &&usevisitp&i = 1 %then , sum(s&&cov&i.randomvisitp.&j ) as s&&cov&i.randomvisitp.&j ;
                  , sum(s&&cov&i..&j ) as  s&&cov&i..&j
				  %if &intno = 0   %then  %do;
					,	sum(ncs&&cov&i..&j) as 	ncs&&cov&i..&j  
                        %if &&usevisitp&i = 1 %then , sum( ncs&&cov&i.randomvisitp.&j) as  ncs&&cov&i.randomvisitp.&j ;
				   %end;
              %end;
           %end; 
       %end; 
 
       /**/ 			            
     , min( intervened ) as intervened_min, min( averinterv ) as averinterv_min  
            
            %if &minimalistic = no %then %do;
                %do i = 1 %to &ncov;
                   %do j = 1 %to %eval(&timepoints);
                      %if &&usevisitp&i = 1 %then , min(s&&cov&i.randomvisitp.&j) as s&&cov&i.randomvisitp.&j._min ; 
                     , min(s&&cov&i..&j) as  s&&cov&i..&j._min 
					  %if &intno = 0  %then %do;
					,	min(ncs&&cov&i..&j ) as ncs&&cov&i..&j._min  
                           %if &&usevisitp&i = 1 %then , min( ncs&&cov&i.randomvisitp.&j) as ncs&&cov&i.randomvisitp.&j._min ;
					   %end;
                  %end;
               %end; 
            %end;                     
    /**/
		              			               			                              
    /***/
    , max ( intervened ) as intervened_max, max( averinterv ) as averinterv_max  
            
            %if &minimalistic = no %then %do;
                %do i = 1 %to &ncov;
                   %do j = 1 %to %eval(&timepoints);
                      %if &&usevisitp&i = 1 %then , min(s&&cov&i.randomvisitp.&j) as s&&cov&i.randomvisitp.&j._min ; 
                     , max(s&&cov&i..&j) as  s&&cov&i..&j._max 
					  %if &intno = 0  %then %do;
					 ,	max(ncs&&cov&i..&j ) as ncs&&cov&i..&j._max  %if &&usevisitp&i = 1 %then , max( ncs&&cov&i.randomvisitp.&j) as ncs&&cov&i.randomvisitp.&j._max ;
					   %end;
                  %end;
               %end; 
            %end;                      
    /***/
		    from simulated&intno group by &outc ;

     create table mytest2 as 
	 select  sum (intervened)/&nsimul  as intervened , sum(averinterv) / &nsimul  as  averinterv  
		                    
        %if &minimalistic = no %then %do;
            %do i = 1 %to &ncov;
               %do j = 1 %to %eval(&timepoints );
                  %if &&usevisitp&i = 1 %then , sum(s&&cov&i.randomvisitp.&j ) / &nsimul as s&&cov&i.randomvisitp.&j ;
                  , sum(s&&cov&i..&j ) / &nsimul as  s&&cov&i..&j
				  %if &intno = 0   %then  %do;
					,	sum(ncs&&cov&i..&j)/&nsimul as 	ncs&&cov&i..&j  
                        %if &&usevisitp&i = 1 %then , sum( ncs&&cov&i.randomvisitp.&j)/nsimul as  ncs&&cov&i.randomvisitp.&j ;
				   %end;
              %end;
           %end; 
       %end; 
 
       /**/ 			            
     , min( intervened_min ) as intervened_min, min( averinterv_min ) as averinterv_min  
            
            %if &minimalistic = no %then %do;
                %do i = 1 %to &ncov;
                   %do j = 1 %to %eval(&timepoints);
                      %if &&usevisitp&i = 1 %then , min(s&&cov&i.randomvisitp.&j._min) as s&&cov&i.randomvisitp.&j._min ; 
                     , min(s&&cov&i..&j._min) as  s&&cov&i..&j._min 
					  %if &intno = 0  %then %do;
					,	min(ncs&&cov&i..&j._min ) as ncs&&cov&i..&j._min  
                           %if &&usevisitp&i = 1 %then , min( ncs&&cov&i.randomvisitp.&j._min) as ncs&&cov&i.randomvisitp.&j._min ;
					   %end;
                  %end;
               %end; 
            %end;                     
    /**/
		              			               			                              
    /***/
    , max ( intervened_max ) as intervened_max, max( averinterv_max ) as averinterv_max  
            
            %if &minimalistic = no %then %do;
                %do i = 1 %to &ncov;
                   %do j = 1 %to %eval(&timepoints);
                      %if &&usevisitp&i = 1 %then , min(s&&cov&i.randomvisitp.&j._max) as s&&cov&i.randomvisitp.&j._max ; 
                     , max(s&&cov&i..&j._max) as  s&&cov&i..&j._max 
					  %if &intno = 0  %then %do;
					 ,	max(ncs&&cov&i..&j._max ) as ncs&&cov&i..&j._max  %if &&usevisitp&i = 1 %then , max( ncs&&cov&i.randomvisitp.&j._max) as ncs&&cov&i.randomvisitp.&j._max ;
					   %end;
                  %end;
               %end; 
            %end;  
      from mytest ; 
      quit;

	  data mytest3 ;
	  set mytest (keep = &outc._level freq ) ;
	  s&outc = round(freq *100/&nsimul ,0.01) ;
	  run ;

	  proc transpose data = mytest3 out = mytest4 (drop = _name_ ) prefix = s&outc._ ;
	  var s&outc ;
	  id &outc._level ;
	  run;

	 
	  data stats_holder  ;
	  merge mytest4 mytest2 ;
	  run;

	  proc datasets library = work nolist ;
	  delete mytest mytest2 mytest3 mytest4 ;
	  quit;
			               
			             
			            





	%end;
    %else %do;


     proc means data=simulated&intno  noprint /*simulated&intno is a function*/; 
     
      output out=stats_holder
             mean (%if &outctype=binsurv %then %do;
                        cuminc 
                        %do j = 1 %to  %eval(&timepoints );
                            s&outc.&j
                        %end;
                   %end;
                   %else &outc ;
                  
                   intervened averinterv  
                    
                    %if &minimalistic = no %then %do;
                        %do i = 1 %to &ncov;
                           %do j = 1 %to %eval(&timepoints );
                              %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j ;
                              s&&cov&i..&j 
							  %if &intno = 0 /*** AND %bquote(&censor)^=  ****/ %then  %do;
									ncs&&cov&i..&j %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
							   %end;
                          %end;
                       %end; 
                   %end;)  = %if &outctype=binsurv %then %do;
                                pd  
                                %do j = 1 %to %eval(&timepoints);
                                    s&outc.&j

                                %end;
                             %end;
                             %else s&outc ;
                            intervened averinterv 
                       
                       %if &minimalistic = no %then %do;
                           %do i = 1 %to &ncov;
                               %do j = 1 %to %eval(&timepoints);
                                  %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j ;
                                  s&&cov&i..&j
								  %if &intno = 0 /*** AND %bquote(&censor)^= ***/ %then %do;
										ncs&&cov&i..&j %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
								   %end;
                               %end;
                           %end;  
                       %end;
                     
                      
          
            
             min ( %if &outctype=binsurv %then %do;
                        cuminc 
                        %do j = 1 %to %eval(&timepoints);
                            s&outc.&j
                        %end;
                    %end;
                    %else &outc ;
                    intervened averinterv  
                    
                    %if &minimalistic = no %then %do;
                        %do i = 1 %to &ncov;
                           %do j = 1 %to %eval(&timepoints);
                              %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j ; 
                              s&&cov&i..&j
							  %if &intno = 0 /*** AND %bquote(&censor)^=  ***/  %then %do;
									ncs&&cov&i..&j %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
							   %end;
                          %end;
                       %end; 
                    %end;) = %if &outctype=binsurv %then %do;
                                pd_min 
                                %do j = 1 %to %eval(&timepoints);
                                    s&outc.&j._min
                                %end;
                              %end;
                             %else s&outc._min ;
                           intervened_min averinterv_min 
                           
                           %if &minimalistic = no %then %do;
                               %do i = 1 %to &ncov;
                                   %do j = 1 %to %eval(&timepoints);
                                      %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j._min ; 
                                      s&&cov&i..&j._min
									  %if &intno = 0 /**** AND %bquote(&censor)^= ****/ %then  %do ;
											ncs&&cov&i..&j._min %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j._min ;
										%end;
                                   %end;
                               %end;
                           %end;                               
            
              max (%if &outctype=binsurv %then %do;
                        cuminc 
                        %do j = 1 %to %eval(&timepoints);
                          s&outc.&j
                        %end;
                    %end;
                    %else &outc ;
                    intervened averinterv  
                    
                    %if &minimalistic = no %then %do;
                        %do i = 1 %to &ncov;
                           %do j = 1 %to %eval(&timepoints);
                              %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j ; 
                              s&&cov&i..&j
							  %if &intno = 0 /*** AND %bquote(&censor)^= ***/ %then %do;
									ncs&&cov&i..&j %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
							   %end;
                          %end;
                       %end; 
                    %end; ) =%if &outctype=binsurv %then %do;
                                pd_max 
                                %do j = 1 %to %eval(&timepoints);
                                    s&outc.&j._max
                                %end;
                             %end;
                             %else s&outc._max;
                             intervened_max averinterv_max 
                             
                             %if &minimalistic = no %then %do;
                                 %do i = 1 %to &ncov;
                                    %do j = 1 %to %eval(&timepoints);
                                        %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j._max ; 
                                        s&&cov&i..&j._max
										%if &intno = 0 /*** AND %bquote(&censor)^=  ****/ %then %do;
											ncs&&cov&i..&j._max %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j._max ;
										 %end;
                                    %end;
                                %end; 
                             %end;
                        

              ;
               
             
             %if &outctype = binsurv %then %do;
                 output out=surv_holder mean( %do n = 1 %to %eval(&timepoints;
                                                cumincr&n cumsurvr&n cumcompeventr&n  
                                              %end; ) =     %do n = 1 %to %eval(&timepoints);
                                                                risk&n surv&n compevent&n   
                                                            %end;
                          ; 
            %end;
             
            run; 

%end;
            data interv&intno;
            set stats_holder;
            _sample_ = &bsample;
            length int2 $70 ;
            int=&intno;
            int2=&intlabel;
            n= &nsimul /*_FREQ_ */;
            keep int int2 _sample_ n %if &outctype=binsurv %then %do;
                                        pd 
                                        %do j = 1 %to %eval(&timepoints);
                                            s&outc.&j 										
                                        %end;
									
                                     %end;
                                    %else  %if &outctype = cateofu %then %do;
										%do j = 1 %to &outclev ;
											s&outc._&j
										 %end;
									%end ;
                                    %else s&outc ;
                                    intervened averinterv 
                
                %if &minimalistic = no %then %do;
                    %do j = 1 %to %eval(&timepoints);
                        %do i = 1 %to &ncov;
                             %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j ; 
                            s&&cov&i..&j
							%if &intno = 0 /***  AND %bquote(&censor)^=  ***/ %then  %do; 
								ncs&&cov&i..&j %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
							 %end;
                        %end;
                    %end; 
                %end;
                dataname ssize  
                        %if &outctype ne cateofu %then obsp ; 
                ;
           dataname = "&data";
           ssize = &ssize ;
           %if &outctype ne cateofu %then obsp = &obsp ;;
           run;
   
           %if &outctype = binsurv %then %do;

                data surv_tmp&intno;
                set surv_holder ;
                length int2 $70 ;
                _sample_ = &bsample ;
                int = &intno ;
                int2 = &intlabel ;
                n = _freq_ ;
                surv0 = 1 ;
                keep  int int2 _sample_ n surv0
                   %do n = 1 %to %eval(&timepoints);
                       risk&n surv&n compevent&n   
                   %end;
                 ;
               run;

               proc datasets library = work nolist ;
               delete surv_holder ;
               quit;
           %end;
		   %else %if &outctype = cateofu %then %do;
  				data surv_tmp&intno ;
              	set interv&intno ;
              
              	keep _sample_ %do j = 1 %to &outclev ; s&outc._&j %end ; int int2  ;
              run;
		   %end;
           %else %do;
              data surv_tmp&intno ;
              set interv&intno ;
              
              keep _sample_ s&outc int int2  ;
              run;
          %end; 



/*** temp code for censoring ****/
%if &intno = 0   %then %do;
	data interv&intno ;
	merge interv&intno surv_tmp&intno  %if &outctype = binsurv %then (keep = surv1 - surv&timepoints );;
	array surv{&timepoints } ;
	%do myi = 1 %to &ncov ;
		array ncs&&cov&myi {&timepoints } ;
		%if &&usevisitp&myi = 1 %then array ncs&&cov&myi.randomvisitp { &timepoints } ; ;
	%end;

    do j = 2 to &timepoints ; 
	    %if &outctype ne binsurv %then surv[j-1] = 1 ;;
	    %do myi = 1 %to &ncov ;
			ncs&&cov&myi [j ] = ncs&&cov&myi [ j ] / surv[j - 1 ] ; * should this be j or j-1;
			%if &&usevisitp&myi = 1 %then ncs&&cov&myi.randomvisitp [ j ] = ncs&&cov&myi.randomvisitp [ j ] / surv[j - 1] ;;
		%end;
     end;
	 drop j ;
	 %if &outctype = cateofu %then drop surv1 - surv&timepoints ;;
     run;
%end;


/**** end of temp code ********/



       
          %if &bsample = 0 & &print_stats = 1 %then %do;


              data mean_holder ;
              set stats_holder ;
              keep  %if &outctype=binsurv %then %do;
                        pd 
                        %do j = 1 %to %eval(&timepoints);
                            s&outc.&j
                        %end;
                     %end;
					 %else %if &outctype = cateofu %then %do;
					    %do j = 1 %to &outclev ;
							s&outc._&j
						 %end;
					 %end;
                    %else s&outc ;
                  intervened averinterv 
                  
               %if &minimalistic = no %then %do;
                   %do i = 1 %to &ncov;
                     %do j = 1 %to %eval(&timepoints);
                           %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j ; 
                          s&&cov&i..&j
						  %if &intno = 0 /*** AND %bquote(&censor)^= ***/ %then  %do;
								ncs&&cov&i..&j %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
						   %end;
                      %end;
                  %end;      
               %end; 
              ;
           run;
           data max_holder ;
           set stats_holder ;
           keep  %if &outctype=binsurv %then %do;
                    pd_max 
                    %do j = 1 %to %eval(&timepoints);
                        s&outc.&j._max
                    %end;
                  %end;
				  %else %if &outctype = cateofu %then %do;
				  	/* there is no max calculated, should be maximum level */
				  %end;
                 %else s&outc._max ;
                 intervened_max averinterv_max 
 
             %if &minimalistic = no %then %do;
                 %do i = 1 %to &ncov;
                     %do j = 1 %to %eval(&timepoints);
                          %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j._max ;
                         s&&cov&i..&j._max
						 %if &intno = 0 /**** AND %bquote(&censor)^= ***/ %then %do;
								ncs&&cov&i..&j._max %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j._max ;
						 %end;
                      %end;
                 %end;
              %end;
              ;
           run;

           data min_holder ;
           set stats_holder ;
           keep %if &outctype=binsurv %then %do;
                    pd_min 
                    %do j = 1 %to %eval(&timepoints);
                        s&outc.&j._min
                    %end;
                 %end;
				  %else %if &outctype = cateofu %then %do;
				  	/* no minimum level calculated for categorical outcome */
				  %end;
                %else s&outc._min ;
             intervened_min averinterv_min 
             %if &outctype=binsurv %then %do;
             
             %end;
             %if &minimalistic = no %then %do;
                 %do i = 1 %to &ncov;
                     %do j = 1 %to %eval(&timepoints);
                         %if &&usevisitp&i = 1 %then s&&cov&i.randomvisitp.&j._min ;
                         s&&cov&i..&j._min
						 %if &intno = 0 /****  AND %bquote(&censor)^= ***/ %then %do;
								ncs&&cov&i..&j._min %if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j._min ;
						  %end;
                      %end;
                 %end;
            %end;
            ;
           run;

           proc transpose data = mean_holder out= stats(rename = (col1= mean)) ;
           run;

           proc transpose data = min_holder out = stats_min (rename = (col1 = min) drop = _name_ );
           run;

           proc transpose data = max_holder out = stats_max (rename = (col1=max) drop=_name_ );
           run;


          data test ;
          merge stats( rename=( _name_ = variable))  stats_min stats_max ;
         
          label variable ="Variable"
                mean = "Mean"
                min = "Minimum"
                max = "Maximum";
          run;

          title " mean, min, max under intervention &intno ";
          proc print data = test noobs label  ;
          var variable mean min max ;
          run;
          title ;
          proc datasets library = work nolist ;
          delete stats test stats_min stats_max  max_holder min_holder mean_holder;
          quit;


    %end;
  
    data simul ;
    set simul (keep =  _sample_ &newsimulkeeplist);
    run;

    proc datasets library = work nolist ;
    delete stats_holder ;
    quit;

%mend interv;

%macro interv_init (intno= , intlabel=, nintvar= 0, intvisittype = 1 ,intcond = 1=1, intsetup = mynull ,
          intvar1= ,
          inttype1 =, 
          inttimes1 = (-1),
          intpr1 = 1,
          intcov1 = ,
          intmean1 = ,
          intsd1 = ,
          intvalue1 = ,
          intmin1 =. ,
          intmax1 =. ,
          intchg1 = ,  	
		  intcond1= 1=1 ,
          intgrace1= ,		  

intusermacro1=,  

      intvar2 = , inttype2 = , inttimes2 = (-1), intpr2 = 1, intcov2 = ,
          intmean2 = , intsd2 = , intvalue2 = , intmin2 =. , intmax2 =. , intchg2 = ,intcond2= 1=1 , intgrace2= ,
intusermacro2=,
 
          intvar3 = , inttype3 = , inttimes3 = (-1), intpr3 = 1, intcov3 = ,
                intmean3 = , intsd3 = , intvalue3 = , intmin3 =. , intmax3 =. , intchg3 = , intcond3= 1=1 , intgrace3= ,
intusermacro3=,
 
                intvar4 = , inttype4 = , inttimes4 = (-1), intpr4 = 1, intcov4 = ,
                intmean4 = , intsd4 = , intvalue4 = , intmin4 =. , intmax4 =. , intchg4 = , intcond4= 1=1 , intgrace4= ,
intusermacro4=,
 
          intvar5 = , inttype5 = , inttimes5 = (-1), intpr5 = 1, intcov5 = ,
                intmean5 = , intsd5 = , intvalue5 = , intmin5 =. , intmax5 =. , intchg5 = , intcond5= 1=1 , intgrace5= ,
intusermacro5=,
 
          intvar6 = , inttype6 = , inttimes6 = (-1), intpr6 = 1, intcov6 = ,
                intmean6 = , intsd6 = , intvalue6 = , intmin6 =. , intmax6 =. , intchg6 = , intcond6= 1=1 , intgrace6= ,
intusermacro6=,
 
          intvar7 = , inttype7 = , inttimes7 = (-1), intpr7 = 1, intcov7 = ,
                intmean7 = , intsd7 = , intvalue7 = , intmin7 =. , intmax7 =. , intchg7 = , intcond7= 1=1 , intgrace7= ,
intusermacro7=,
 
          intvar8 = , inttype8 = , inttimes8 = (-1), intpr8 = 1, intcov8 = ,
                intmean8 = , intsd8 = , intvalue8 = , intmin8 =. , intmax8 =. , intchg8 = , intcond8= 1=1 , intgrace8= ,
      intusermacro8 = ,
      
      timesnotelig = -1 

    );

    %************ SIMULATION TO GET CUMULATIVE INCIDENCE UNDER INTERVENTION;

    %if &printlogstats = 1 %then %put  Initializing intervention data view  &intno;

    %local f g h  i icov j  k l n intcovmap ;
     
    %do i = 1 %to &nintvar ;
	    %local intcovmap&i ;
		%let intcovmap&i = -1 ;
		%do covind = 1 %to &ncov ;
			%if &&intvar&i = &&cov&covind %then %let intcovmap&i = &covind ;
		%end ;
		/* %put FOR INTVAR &i : &&intvar&i intcovmap&i  = &&intcovmap&i ; */ 
		%if &&intcovmap&i = -1 %then %put ERROR WITH MAPPING INTVAR TO COV NUMBER ;
   %end;
 

    data simulated&intno   ( keep = &fixedcov
                                    intervened averinterv 
                                     
                                    %if &minimalistic = no %then %do;
                                        %do i = 1 %to &ncov;
                                            %do j = 1 %to &timepoints;
                                               s&&cov&i..&j
                                               %if &&usevisitp&i = 1 %then   s&&cov&i.randomvisitp.&j ;
											   %if (&intno = 0 /* AND %bquote(&censor) ^= */ ) %then  %do;
                                                        ncs&&cov&i..&j  
														%if &&usevisitp&i = 1 %then ncs&&cov&i.randomvisitp.&j ;
											   %end;
                                            %end;
                                        %end; 
                                    %end;
                                     
                                    %if &outctype = binsurv %then %do;
                                         cuminc 
                                         %do n = 1 %to &timepoints;
                                             s&outc.&n  cumincr&n cumsurvr&n cumcompeventr&n   
                                         %end;
                                    %end;   
                                    %else &outc ;
                                    %if &hazardratio = 1 %then _censor _newtime  ;
                                    
                                    
                                     )   /view = simulated&intno  ;
 

        %*Merge with regression coefficients;
        merge simul _betar_  _seedr_ _covbounds_ 
                    %if &hazardratio = 1 %then   _calchazard_   ; 
                      end=_end_;
        by _sample_;
        
        call streaminit(_seedr_);
        drop _seedr_ ; 

         
        %if &hazardratio = 1 %then %do;
            
            _censor = . ;
            _newtime = .;
            
        %end;
        %else %do;
           calchazard = 0 ;
        %end;
        %*Retain marker of out of range variables generated;
        %do i = 0 %to &ncov;
            %if &&cov&i.otype=3 or &&cov&i.otype=4   or  &&cov&i.otype=-1  %then %do;
                retain &&cov&i.._limit 0;
            %end;
        %end;

        %*Retain marker for intervention;
        retain totinterv 0;
        %do i = 1 %to &nintvar;
            retain &&intvar&i.._totinterv 0;

			array intervenedk_&&intvar&i {0:%eval(&timepoints - 1) } ;
        %end;
            
        _intno_ = &intno ;
        %* Arraying beta coefficients;
        uno=1;

        array aboutc boutc00-boutc&dimoutc;
		%if &outctype = conteofu4 %then %do ; 			                 
             array abzoutc boutcz_00-boutcz_&dimoutc;                 
        %end;
		%else %if &outctype = cateofu %then %do;
		     %do lev = 1 %to %eval(&outclev - 1);
                    array about_&lev about_&lev._00-about_&lev._&dimoutc;
             %end;
		%end;
        array s&outc {0:%eval(&timepoints - 1) }   ;
        %if %bquote(&compevent)^= AND &compevent_cens = 0  %then %do;
            array abcompevent bcompevent00-bcompevent&dimcompevent;
           
        %end;
        array scompevent {0:%eval(&timepoints - 1)} ;
       
        
        %do i = 0 %to &ncov;

            %if &&usevisitp&i = 1 %then %do;
                    array abvisitp&&cov&i bvar&i.visitp_00-bvar&i.visitp_&&dimvar&i;  /* * include the coefficients for visit process ; */
            %end;
            %if &&cov&i.otype=1 or &&cov&i.otype=2 or &&cov&i.otype=3 or &&cov&i.otype=4 or &&cov&i.otype=6 
                    or &&cov&i.otype=7 or &&cov&i.otype = -1 %then %do;
                array ab&&cov&i bvar&i._00-bvar&i._&&dimvar&i;
            %end;
            %if &&cov&i.otype=4 %then %do; 
                %if %bquote(&&cov&i.class)^= %then %do;
                    array abz0&&cov&i bvar&i.z0_00-bvar&i.z0_&&dimvar&i;
                    array abz1&&cov&i bvar&i.z1_00-bvar&i.z1_&&dimvar&i; 
                %end;
                %else %do;
                    array abz&&cov&i bvar&i.z_00-bvar&i.z_&&dimvar&i; 
                %end;
            %end;
            %if &&cov&i.otype=5 %then %do;
                %do lev = 1 %to %eval(&&cov&i.lev - 1);
                    array ab&&cov&i.._&lev bvar&i._&lev._00-bvar&i._&lev._&&dimvar&i;
                %end;
            %end;
           
            %if &&cov&i.otype=-1 %then %do; 
                   %&&cov&i.barray ;
           %end;                                    
       %end;  


   
        %*Arraying simulated variables   (X for XB);
        
            %let i = ;
            %let h=_l1 ;
            %let g=_l2 ;
            %let f=_l3 ;
            
    

            %do n=0 %to &ncov;
              
                %if &&cov&n.otype ne 0 %then %do;
                 
                   %if (&&cov&n.ptype=lag1cumavg | &&cov&n.ptype=lag2cumavg   ) & (%bquote(&&cov&n.cumint) ^= )%then %do;
                        array s&&cov&n..inter{0:%eval(&timepoints - 1)} ;
                   %end;
                   %if &&usevisitp&n = 1 %then array s&&cov&n.randomvisitp{0:%eval(&timepoints - 1)};;
                   array as&&cov&n
                        uno  &&cov&n.fixedcov  
                         %let changemtype = 0 ;
                         %if %upcase(&&cov&n.mtype)  = SOME %then %do;
                            %let changemtype = 1 ;
                            %let cov&n.mtype = all ;  
                        %end; 
                        %listpred(main,&n,0,&ncov,&f,&g,&h,&i) %listpred(contemp,,0,&n-1,&f,&g,&h,&i)  
                        %if &changemtype = 1 %then %let cov&n.mtype = some ;                                                   
                        &&cov&n.addvars             
                       %qcmpres(%interxarrays(main,cov&n)) ;             
                %end;
                                                                      

                %if &&cov&n.otype=-1 %then %do; 
                            %&&cov&n.sarray ;
                %end;

                array s&&cov&n {0:%eval(&timepoints - 1) }   ;
                s&&cov&n [0] = &&cov&n ;

				%if &intno = 0 /** AND %bquote(&censor) ^= ***/ %then %do;
					array ncs&&cov&n {0:%eval(&timepoints - 1) } ;
					%if &&usevisitp&n = 1 %then array ncs&&cov&n.randomvisitp{0:%eval(&timepoints - 1)};;
				%end;

          %end;  /* n */

 


 /*arraying time i outcome predictors */
		     %if &usehistory_eof = 0 %then %do;
              array as&outc uno &fixedcov %listpred(main,,0,&ncov,&g,&h,&i,&i)   
                     &eventaddvars 
                    %qcmpres(%interxarrays(main,outc))   ;
			  %end;
			  %else %if &usehistory_eof = 1 %then %do;
				array as&outc uno &fixedcov %listpred(eof ,,0,&ncov,&g,&h,&i,&i)   
                     &eventaddvars 
                    %qcmpres(%interxarrays(main,outc))   ;
			  %end;
                         
     
 /*same for censoring, if simulated*/
              %if %bquote(&compevent)^= AND &compevent_cens = 0  %then %do;
                    array ascompevent uno &fixedcov %listpred(main,,0,&ncov,&g,&h,&i,&i)                                                                           
                         &compeventaddvars 
                        %qcmpres(%interxarrays(main,compevent))   ;
              %end;

              

     
                       
            %&intsetup ; 
            
            %*Setting up censoring simulation;
            /***/
            %if &hazardratio= 1 %then %do;
               array Uoutc{0:%eval(&timepoints-1)} ;
               %if %bquote(&compevent)^= AND &compevent_cens = 0  %then array Ucompevent{0:%eval(&timepoints - 1)} ;;
            %end;
            
            /****/
            scensl = 0; /* * need to initialize simulated censoring for lag1 at baseline time ; */
            mygood = 1 ;
            array intervenedk{0:%eval(&timepoints - 1)} ;              
            suminterv = 0; /* count of any intervention for one person across all time-points */
            elig_persontime = 0;
 
            do &time  = 0 to %eval(&timepoints - 1) ;
                &time._l1 = &time - 1;
                &time._l2 = &time - 2 ;
                &time._l3 = &time - 3 ;
                _ind2 = 1;
                _ind3 = 1;
                _ind4 = 1;
                if &time < 2 then _ind2 = 0 ;
                if &time < 3 then _ind3 = 0 ;
                if &time < 4 then _ind4 = 0 ;

                
                %if &hazardratio = 1 %then %do;
                    Uoutc[&time] = rand('uniform');
                    %if %bquote(&compevent)^= AND &compevent_cens = 0  %then Ucompevent[&time] = rand('uniform');;
                %end;
                /************************************************** 
                need to re-create all the predictors since we have shifted the variables at the end of the loop and need to recreate any predictors
                like the splines for the lagged values ,
                NOTE : THIS WILL PRODUCE THE WRONG PREDICTORS FOR THE CURRENT TIME AND IS ONLY CALLED FOR UPDATING THE LAGGED VALUES, WHICH WERE UPDATED
                       AT THE END OF THE LOOP. THE CURRENT-TIME PREDICTORS FROM COVx WILL ONLY BE USED IN THE MODELS FOR COVy WHERE x<y. FOR y<=x, THE LAGGED
                       VALUES OF COVx WILL BE USED. WHEN A NEW VALUE OF COVx IS GENERATED IN SIMVAR, GENPRED WILL BE CALLED TO UPDATE THE PREDICTORS. DUE
                       TO HOW GENPRED WORKS, THE PREDICTORS BASED ON THE LAGGED VALUES WILL BE RECREATED.
                 *******************************************/
                *running genpred for lagged valuees ;
                %do i = 0 %to &ncov ;                     
                     %genpred(sim,lagtype=2);                    
                %end;       
                %do i = 1 %to &ncov ;
                       /* the interactions for cov&i may depend on variables derived from cov&j where i < j, need a separate loop */
                       %interactions(cov&i,,2,createvar); /* create any interactions. these are for the covariate models, and use lagged valued */
                %end; 

                %* Generating   covariates;
                %simvar;                
                if &time = 0  then intervened = 0;                 
                intervenedk[&time] = 0;                
                x = 0;


                %*Looping over intervention variables;
                %do i = 1 %to &nintvar;    

				


 					intervenedk_&&intvar&i [ &time ] = 0 ; 
                    %*Limiting to intervention times;
                    if &time in (&&inttimes&i) then do;                        
                        %*Limiting to when condition is met;
                            if &intcond then do;                                          
                                if &time in (&timesnotelig) then x = x ; 
                                else x = x + 1;
				
								
                                %*Intervention Type 1: Static deterministic intervention; 
                                %if &&inttype&i = 1 %then %do;                    
                                    if  &&intvar&i ^=&&intvalue&i and (rand('uniform')<=&&intpr&i) then do;
                                        &&intvar&i=&&intvalue&i;
                                        &&intvar&i.._totinterv = &&intvar&i.._totinterv + 1;
                                        s&&intvar&i [ &time] = &&intvar&i ;
                                        intervened = 1;
                                        intervenedk[&time] = 1;
									*	intervenedk_&&intvar&i [ &time ] = 1 ;
                                        totinterv = totinterv + 1;
                                    end;
                                %end;                            
                                %*Intervention Type 2: Threshold intervention; 
                                %else %if &&inttype&i = 2 %then %do;
                                    if %bquote(&&intmax&i)^=. then do;
                                        if &&intvar&i >&&intmax&i and (rand('uniform')<=&&intpr&i) then do;
                                            &&intvar&i=&&intmax&i;
                                            s&&intvar&i [ &time] = &&intvar&i ;
                                            &&intvar&i.._totinterv = &&intvar&i.._totinterv + 1;
                                            intervened = 1;
                                            intervenedk[&time] = 1;
										*	intervenedk_&&intvar&i [ &time ] = 1 ;
                                            totinterv = totinterv + 1;
                                        end;
                                    end;
                                    if %bquote(&&intmin&i)^=. then do;
                                        if &&intvar&i<&&intmin&i and (rand('uniform')<=&&intpr&i) then do;
                                            &&intvar&i=&&intmin&i;
                                            s&&intvar&i [ &time] = &&intvar&i ;
                                            &&intvar&i.._totinterv = &&intvar&i.._totinterv + 1;
                                            intervened = 1;
                                            intervenedk[&time] = 1; 
										*	intervenedk_&&intvar&i [ &time ] = 1 ;
                                            totinterv = totinterv + 1;
                                        end;
                                    end;
                                %end;
                                %*Intervention Type 3: Fixed Change; /*different int types were not working with else ifs so changed to if, still not sure why - change JGY*/
                                %else %if &&inttype&i = 3 %then %do;
                                    if (rand('uniform')<=&&intpr&i) then do;
                                        &&intvar&i = &&intvar&i + &&intchg&i;
                                        s&&intvar&i [ &time] = &&intvar&i ;
                                        &&intvar&i.._totinterv = &&intvar&i.._totinterv + 1;
                                        intervened = 1;
                                        intervenedk[&time] = 1; 
									*	intervenedk_&&intvar&i [ &time ] = 1 ;
                                        totinterv = totinterv + 1;
                                    end;      
                                %end; 

                                %*Intervention Type 4: Carry forward; 
                                %else %if &&inttype&i = 4 %then %do;
                                    if ( &&intvar&i  ne &&intvar&i.._l1) and (rand('uniform')<=&&intpr&i) then do;
                                        &&intvar&i  = &&intvar&i.._l1;
                                        s&&intvar&i [ &time] = &&intvar&i ;
                                        &&intvar&i.._totinterv = &&intvar&i.._totinterv + 1;
                                        intervened = 1;
                                        intervenedk[&time] = 1;
									*	intervenedk_&&intvar&i [ &time ] = 1 ;
                                        totinterv = totinterv + 1;
                                    end;      
                                %end;
							
								%else %if &&inttype&i = 5 %then %do; /* want to mimic the cutgrace macro */
                                   
                                   if &time = 0 then do ;
                                        _thresh&i = 0 ;
                                        _ts_thresh&i = 0 ;                                        
                                    end;
                                           
                                   * call streaminit(123);
                                    _testp_ = rand('uniform');
                                                                                           
                                    _tgrace&i = &time - &&intgrace&i;

                                    if _thresh&i = 0 and ( &&intcond&i ) then _thresh&i = 1;
                     
                                    if _thresh&i = 0 then &&intvar&i = 0 ;
										

									/* change this so that for grace period = 6, the _p_ are 1/6 then 1/ 5 ... 1/2 and then 1/1 for last period */
                                    if _thresh&i = 1 then do ;
                                       *** _p_ = (_ts_thresh&i+1) / &&intgrace&i ;
                                        /* 1/grace ,1/(grace - 1) ,...,1/(grace - 4), 1/(grace - 5) */
                                        _p_ = 1/(&&intgrace&i - _ts_trhesh&i) ;
                                        if  &&intvar&i.._l1 = 1 then &&intvar&i = 1 ;
											
                                        if _ts_thresh&i < &&intgrace&i and &&intvar&i = 0 then do ;                       
                                            if _testp_ < _p_ then	&&intvar&i = 1 ; 											
                                        end;
                                        else &&intvar&i = 1 ;
											

                                        _ts_thresh&i = _ts_thresh&i + 1; 
                                    end;
                

                                %end;      
								%else %if &&inttype&i = 6 %then %do; /* want to mimic the cutgrace macro */
                                                                      
                                   if &time = 0 then do ;
                                        _thresh&i = 0 ;
                                        _ts_thresh&i = 0 ;                                        
                                    end;
                                           
                                    *call streaminit(123);
                                    _testp_ = rand('uniform');
                                                                                           
                                    _tgrace&i = &time - &&intgrace&i;

                                    if _thresh&i = 0 and ( &&intcond&i ) then _thresh&i = 1;
                     
                                    if _thresh&i = 0 then do;
									    if &&intvar&i = 1 then do;
											&&intvar&i = 0 ;
											intervened = 1;
                                        	intervenedk[&time] = 1;
                                        	intervenedk[&time] = 1; 
                                        	totinterv = totinterv + 1;
										 end;
									 end;

									/* change this so that for grace period = 6, the _p_ are 1/6 then 1/ 5 ... 1/2 and then 1/1 for last period */
                                    if _thresh&i = 1 then do ;
                                       *** _p_ = (_ts_thresh&i+1) / &&intgrace&i ;
                                        /* 1/grace ,1/(grace - 1) ,...,1/(grace - 4), 1/(grace - 5) */
                                        _p_ = 1/(&&intgrace&i - _ts_trhesh&i) ;
                                        if  &&intvar&i.._l1 = 1 then do;
 											if &&intvar&i = 0 then do;
												&&intvar&i = 1 ;
												intervened = 1;
                                        		intervenedk[&time] = 1;
                                        		intervenedk[&time] = 1; 
                                        		totinterv = totinterv + 1;
											 end;
										 end;
                                        if _ts_thresh&i < &&intgrace&i and &&intvar&i = 0 then do ;                       
                                            if _testp_ < _p_ then do;
												&&intvar&i = 1 ; 
												intervened = 1;
                                        		intervenedk[&time] = 1;
                                        		intervenedk[&time] = 1; 
                                        		totinterv = totinterv + 1;
											 end;
                                        end;
                                        else  do;
											if &&intvar&i = 0 then  do ;
												&&intvar&i = 1 ;
												intervened = 1;
                                        		intervenedk[&time] = 1;
                                        		intervenedk[&time] = 1; 
                                        		totinterv = totinterv + 1;
                                             end;
										 end;

                                        _ts_thresh&i = _ts_thresh&i + 1; 
                                    end;
                

                                %end;                    
/*********/
                                %*Intervention Type 7: Randomly Assign for Observed Dist;
                                %else %if &&inttype&i = 7 %then %do;
                                    if (rand('uniform')<=&&intpr&i) then do;
                                        &&intvar&i = r_&&intvar&i.._&i;
                                        s&&intvar&i [ &time] = &&intvar&i ;
                                        &&intvar&i.._totinterv = &&intvar&i.._totinterv + 1;
                                        intervened = 1;
                                        intervenedk[&time] = 1;                                       
									*	intervenedk_&&intvar&i [ &time ] = 1 ;
                                        totinterv = totinterv + 1;
                                    end;      
                                %end; 

                                %*Intervention Type 8: Randomly Assign for Observed Dist, Conditional;
                                %else %if &&inttype&i = 8 %then %do;
                                    if %bquote(&&intmax&i)^=. then do;
                                        if &&intvar&i >&&intmax&i and (rand('uniform')<=&&intpr&i) then do;
                                            &&intvar&i..&k = r_&&intvar&i.._&i;
                                            s&&intvar&i [ &time] = &&intvar&i ;
                                            &&intvar&i.._totinterv = &&intvar&i.._totinterv + 1;
                                            intervened = 1;
                                            intervenedk[&time] = 1;
										*	intervenedk_&&intvar&i [ &time ] = 1 ;
                                            totinterv = totinterv + 1;
                                        end;
                                    end;
                                    if %bquote(&&intmin&i)^=. then do;
                                        if &&intvar&i<&&intmin&i and (rand('uniform')<=&&intpr&i) then do;
                                            &&intvar&i = r_&&intvar&i.._&i;
                                            s&&intvar&i [ &time] = &&intvar&i ;
                                            &&intvar&i.._totinterv = &&intvar&i.._totinterv + 1;
                                            intervened = 1;
                                            intervenedk[&time] = 1;
										*	intervenedk_&&intvar&i [ &time ] = 1 ;
                                            totinterv = totinterv + 1;
                                        end;
                                    end;
                                %end;   
/*********/	                                    								
                                %*Intervention Type -1: User defined intervention  ;
                                %else %if &&inttype&i = -1 %then %do;                                               
                                    %&&intusermacro&i ;  								 
                                %end; 


                            end;  /* intcond */ 							
                        end; /* inttimes */  
					 
                    %end;  /* nintvar */                    
                    if x >= 1 then elig_persontime = elig_persontime + 1;
                    suminterv = suminterv + intervenedk[&time] ;

                    if mygood  = 1 then do ;
                        %*Generating derived variables (categories, etc);
                        /* careful use of i, used in genpred, inherit from outside definition */
                        /* 10-24-2016 due to ptype can now involve more than one modeled variable, we need to call genpred for all variables, not
                            only those we have intervened on, example x has ptype = lag1cumavg and a cumint treatment variable. If we intervened on 
                            the treatment variable, then the predictors involving x will need to be re-calculated */
						
                        %do i = 1 %to &ncov ;
                              %genpred(sim,lagtype=1);
                        %end;
                        

                        /* outcome interactions,  for outcome we use factors at current time. */
                        %interactions(outc,,1,createvar);
                        %if %bquote(&compevent)^= AND &compevent_cens = 0  %then %do;
                            %interactions(compevent,,1,createvar);
                        %end;
                         
                      
						%if &usehistory_eof = 1 %then %do ;
							if &time = %eval(&timepoints - 1) then do ;
                                %do i = 1 %to &ncov;
 									array a&&cov&i{0:%eval(&timepoints - 1) }  %do ii = 0 %to %eval(&timepoints - 1) ;
                                                                                     &&cov&i.._&ii._eof 
																				%end; ;
									do _time_ = 0 to %eval(&timepoints - 1) ;
         	   							a&&cov&i [_time_ ] = s&&cov&i[ _time_ ]  ;
									end;
									%genpred(type=sim,useeof = 1 ) ;

                                %end; 

							end;
					

						%end ;
                       

                        %if &outctype=binsurv  %then %do;

						 	m&outc = 0 ;
                        	do i=1 to dim(aboutc);
                            	m&outc =sum(m&outc,aboutc(i)*as&outc(i));
                        	end;

                            if &outcwherenosim then do; * at no sim time ; 
                                /* evaluate user defined macro  because of some condition other than its too early to start simulating*/
                                %if %bquote(&outcnosimelsemacro)^= %then  %&outcnosimelsemacro ;                                                        
                            end;
                            else do ;
                                p&outc = 1/(1+exp(-m&outc));                            
                            end;
                            s&outc[&time] = p&outc ;
                            %if &hazardratio = 1 %then %do;
                                if calchazard = 1 then do ;  
                                    if Uoutc[&time] <= p&outc then _censor = 1;
                                    
                                    if  _censor = 1 then do ; 
                                        _newtime = &time ;                                         
                                        mygood = 0 ;  /* do not want to simulate any further covariate history    */
                                    end;                    
                                end; 

                            %end;
                        %end;
                        %if &outctype = bineofu %then %do;


 							m&outc = 0 ;
                        	do i=1 to dim(aboutc);
                            	m&outc =sum(m&outc,aboutc(i)*as&outc(i));
                        	end;

                            if &time = %eval(&timepoints -1) then do;
                                if &outcwherenosim then do; * at no sim time ; 
                                    /* evaluate user defined macro  because of some condition other than its too early to start simulating*/
                                    %if %bquote(&outcnosimelsemacro)^= %then %&outcnosimelsemacro ;                                        
                                end;
                                else do ;
                                    &outc = 1/(1+exp(-m&outc));
                                end;

                                 
                            end ;
                            else   &outc = . ;
                            s&outc [ &time] = &outc ;
                        %end;                      
                        /* FOR A CONTINUOUS OUTCOME MEASURED AT EOF  */ 
                        %if &outctype=conteofu or &outctype = conteofu2 or &outctype=conteofu3 %then %do;                                              
                            if &time = %eval(&timepoints - 1) then do ;

 								m&outc = 0 ;
                        		do _i=1 to dim(aboutc);
                            		m&outc =sum(m&outc,aboutc(_i)*as&outc(_i));
                        		end;

                                if &outcwherenosim then do; * at no sim time ; 
                                    /* evaluate user defined macro  because of some condition other than its too early to start simulating*/
                                    %if %bquote(&outcnosimelsemacro)^= %then %&outcnosimelsemacro ;                              
                                end;
                                else do ;
                                    /* outc_min and outc_max are read into the data from the _covbounds_ data set and can change based on the bootstrap sample */
                                    &outc = m&outc ;                                                 
                                    if  &outc < &outc._min then do;
                                        &outc = &outc._min; 
                                        &outc._limit = &outc._limit +1;                  
                                    end;                       
                                    if &outc > &outc._max then do;
                                        &outc = &outc._max; 
                                        &outc._limit = &outc._limit +1;
                                    end;
                                end;
                                 
                            end;                        
                            else do;
                                &outc. = . ;
                            end;

                            s&outc [ &time] = &outc ;
                        %end;
						%else %if &outctype=conteofu4 %then %do;                                              
                            if &time = (&timepoints - 1) then do ;

                                if &outcwherenosim then do; * at no sim time ; 
                                    /* evaluate user defined macro  because of some condition other than its too early to start simulating*/
                                    %if %bquote(&outcnosimelsemacro)^= %then %&outcnosimelsemacro ;                              
                                end;
                                else do ;
                                    /* outc_min and outc_max are read into the data from the _covbounds_ data set and can change based on the bootstrap sample */

										z&outc = 0 ;
                        				do i=1 to dim(abzoutc);
                         	   				z&outc =sum(z&outc,abzoutc(i)*as&outc(i));
                        				end;
										p&outc = 1/(1.0+exp(-1*z&outc ) );
										if U&outc <= p&outc then &outc = 1 ;
										else &outc = 0 ;

									    m&outc = 0 ;
			                        	do i=1 to dim(aboutc);
			                            	m&outc =sum(m&outc,aboutc(i)*as&outc(i));
			                        	end;

										if &outc = 1 then &outc = exp(m&outc) ;
									
                                                                                      
	                                    if  0 < &outc < &outc._min then do;
	                                        &outc = &outc._min; 
	                                        &outc._limit = &outc._limit +1;                  
	                                    end;                       
	                                    if &outc > &outc._max then do;
	                                        &outc = &outc._max; 
	                                        &outc._limit = &outc._limit +1;
	                                    end;
                                end;
                                 
                            end;                        
                            else do;
                                &outc. = . ;
                            end;

                            s&outc [ &time] = &outc ;
                        %end;
						%else %if &outctype = cateofu %then %do;
 
                            &outc = . ;
                            %do lev = 1 %to %eval(&outclev - 1);
                                if &outc = . then do;
                                    moutc_lev = 0 ;
                                    do j=1 to dim(about_&lev);
                                        moutc_lev=sum(moutc_lev,about_&lev(j)*as&outc(j));
                                    end;
                                    poutc_&lev=1/(1+exp(-moutc_lev));
                                    if Uoutc_&lev<=poutc_&lev  then &outc=&lev;
                                end;
                           %end;

                           if &outc = . then &outc = &outclev;                        


						%end;



                        %*Censoring at time k;
                        %if %bquote(&compevent) = OR &compevent_cens = 1 %then %do;                         
                            scompevent[&time] = 0 ;
                        %end;                   
                        %else %if %bquote(&compevent)^= AND &compevent_cens = 0  %then %do;
                             
                            mcompevent = 0 ;
                            do i=1 to dim(abcompevent);
                                mcompevent=sum(mcompevent,abcompevent(i)*ascompevent(i));
                            end;

                            if &compeventwherenosim then do; * at no sim time ; 
                                /* evaluate user defined macro  because of some condition other than its too early to start simulating*/
                                %if %bquote(&compeventnosimelsemacro) ^= %then %&compeventnosimelsemacro ;                                                               
                            end;
                            else do ;
                                pcompevent = 1/(1+exp(-mcompevent));
                            end;
                            scompevent[&time] = pcompevent ;
                            %if &hazardratio = 1 %then %do;

                                if calchazard = 1 then do;
                                    if Ucompevent[&time] <= pcompevent then _censor = 2;
                                    

                                    if _censor = 2 then do ;
                                        _newtime = &time ;                                        
                                        mygood = 0 ;  /* do not want to simulate any further covariate history   */
                                    end;
                                end;
                            %end;
                        %end;
                        
                            pcensl=0;                                   
                        
                    end; /* mygood = 1 */
                    else do ;
                        p&outc = . ;
                    end;
                 

                    %if &hazardratio = 1   %then %do;
                        /*  made it to end-of-followup without event or censor. need to censor due to eof */
                        if calchazard = 1 and mygood = 1 and &time = ( &timepoints -1) then do ; 
                            _censor = 0 ;                            
                            _newtime = &time ;
                        end;
                    %end;
                    /* increment lagged variables for the next time point */
                    %do icov = 1 %to &ncov ;
                        &&cov&icov.._l3 = &&cov&icov.._l2 ;
                        &&cov&icov.._l2 = &&cov&icov.._l1 ;
                        &&cov&icov.._l1 = &&cov&icov ;
                        %if &&usevisitp&icov = 1  %then %do;
                            ts_last_&&cov&icov.._l1 = ts_last_&&cov&icov ;
                        %end;
                        %if &&cov&icov.ptype = tsswitch1   %then ts&&cov&icov.._l1_inter = ts&&cov&icov.._inter ;;
                        %if &&cov&icov.ptype=cumavg |  &&cov&icov.ptype=lag1cumavg | &&cov&icov.ptype=lag2cumavg  |
                            &&cov&icov.ptype=cumavgcat |  &&cov&icov.ptype=lag1cumavgcat | &&cov&icov.ptype=lag2cumavgcat
                            %then %do;
                                &&cov&icov.._cumavg_l3 = &&cov&icov.._cumavg_l2 ;
                                &&cov&icov.._cumavg_l2 = &&cov&icov.._cumavg_l1 ;
                                &&cov&icov.._cumavg_l1 = &&cov&icov.._cumavg ;
                                %if %bquote(&&cov&icov.cumint)^= %then %do;
                                      if &&cov&icov.cumint = 1 then &&cov&icov.._normalization_factor  = &&cov&icov.._normalization_factor +1;
                                %end;
                        %end ;
                        %if &&cov&icov.ptype = rcumavg   %then %do;
                               &&cov&icov.._rcumavg_l1 =  &&cov&icov.._rcumavg ;
                        %end;
                    %end ;


            end ; /* time loop */
       

            
             if (&intno ne 0 & elig_persontime ne 0) then averinterv = suminterv / elig_persontime; 
             else averinterv = 0;           
            %if  &outctype=binsurv    %then %do ; 
				 
                   %simcuminc ;                 
            %end;
			%else %do;
				%if &intno = 0 /**** AND %bquote(&censor) ^= ****/ %then %do;
				   /* %put inside creating ncs variables &timepoints , &ncov ; */
				 	do time = 0 to %eval(&timepoints - 1) ;
						%do myi = 1 %to &ncov ;				   
							ncs&&cov&myi [ time] = s&&cov&myi [ time] ;
							%if &&usevisitp&myi = 1 %then ncs&&cov&myi.randomvisitp [ time ] = s&&cov&myi.randomvisitp [ time]  ;;
						%end;
				 	end;
				%end;
			%end;
            
                 
            %*Outputting some summary measures;
            if _end_ then do;

                %*Out of limit generations;
                %do i = 0 %to &ncov;
                    %if &&cov&i.otype=3 or &&cov&i.otype=4   or &&cov&i.otype=-1 %then %do;
                       %if &printlogstats = 1 %then put  %unquote(%str(%'&&cov&i predicted out of range = %')) &&cov&i.._limit;;
                    %end;            
                %end;

                %if &outctype=conteofu %then %do;
                   %if &printlogstats = 1 %then put  %unquote(%str(%'&outc predicted out of range = %')) &outc._limit;;
                %end;
            
                %*Total interventions;
                %if &intno ne 0 %then %do;
                    %do i = 1 %to &nintvar;
                       %if &printlogstats = 1 %then put  %unquote(%str(%'&&intvar&i intervened on = %')) &&intvar&i.._totinterv;;
                    %end;
                %end;
           
            end;
  
    run;




%mend interv_init;



%macro results;

     %************ SUMMARIZING AND OUTPUTTING RESULTS; 
     %local i j k useboot prec ;

     %if (&check_cov_models = 1 OR &rungraphs  ) AND &runnc = 1 %then %do;


          %if &save_raw_covmean = 1  %then %do;

               data &covmeanname._raw ;
               set &covmeanname ;
               run;

          %end;

          %if &outctype = binsurv %then %do;

               data rr ;              
               set &covmeanname (keep = _sample_ _time_   meanoutc   meancompevent ) ;
               by _sample_ _time_ ;
               retain cumsurv  cuminc newn ;
               if first._sample_ then do ;
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
               keep _sample_  _time_  cuminc cumsurv ;
               run;
 

               proc transpose data = rr  out=_stmp_ prefix= obsurv ;
               var cumsurv ;
               id _time_ ;
               by _sample_; 
               run;
               proc transpose data = rr out = _itmp_ prefix= obrisk ;
               var cuminc ;
               id _time_ ;
               by _sample_ ;
               run;

               data &observed_surv ;
               merge _stmp_ _itmp_ ;
               by _sample_ ;
               _int_ = 0 ;
               drop _NAME_ ;
               run;

               proc datasets library = work nolist ;
               delete rr _stmp_ _itmp_ ;
               quit;

          %end; /* binsurv */
          %else %do; /* for all other outcome types, must be eof type outcome  */
               data _obsoutcmean_ &observed_surv ;
               set &covmeanname (keep = _sample_ &time &outc ); *** change for cateofu ;;
               if &time = &timepoints -1  ;  
               run;
 
          %end;

          %do i = 1 %to &ncov ;

               proc transpose data = &covmeanname   out = _mean_&i(drop = _NAME_ ) prefix = &&cov&i ;
               var &&cov&i ;
               by _sample_ ;
               run;

               %if &&usevisitp&i = 1 %then %do;
                    proc transpose data = &covmeanname out = _mean_vp&i (drop = _NAME_) prefix=&&cov&i.randomvisitp ;
                    var &&cov&i.randomvisitp;
                    by _sample_;
                    run;

               %end;
          %end ;

          data &covmeanname ;
          merge %do i = 1 %to &ncov ;
                    _mean_&i 
                    %if &&usevisitp&i = 1 %then _mean_vp&i ;
               %end;
               ;
          by _sample_ ;
          run;

          data _diff_mean ;
          merge &covmeanname  (keep = _sample_ 
               %do i = 1 %to &ncov ;
                    %do k = 1 %to &timepoints ;
                         &&cov&i..&k
                         %if &&usevisitp&i = 1 %then &&cov&i.randomvisitp.&k ;
                    %end;
               %end;)                
               interv0_all(keep = _sample_ 
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
			);

           %do i = 1 %to &ncov ;
                    %do k = 1 %to &timepoints ;
                         d&&cov&i..&k  = &&cov&i..&k  - s&&cov&i..&k ;
                         %if &&usevisitp&i = 1 %then d&&cov&i.randomvisitp.&k = &&cov&i.randomvisitp.&k - s&&cov&i.randomvisitp.&k;;
                    %end;
           %end;
           run;

          proc univariate data = _diff_mean(where = (_sample_ > 0))  noprint  ;
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





          data &covmeanname ;
          merge _diff_mean (where = (_sample_ = 0)) _cov_std2 ; 
          drop  _sample_ ;
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



            %if &print_cov_means = 1 %then %do ;
                    %do i = 1 %to &ncov ;
                         proc print data = &covmeanname noobs label;
                         title "Comparison of means of observed &&cov&i and simulated &&cov&i over &time ";
                         var &time   &&cov&i   s&&cov&i  &&cov&i.._diff   &&cov&i.._stddev  &&cov&i.._lb    &&cov&i.._ub  &&cov&i.._lbp &&cov&i.._ubp           ;               
                         run;    
                         title ;

                         %if &&usevisitp&i = 1 %then %do;
                              proc print data = &covmeanname noobs label;
                              title "Comparison of means of observed &&cov&i.randomvisitp and simulated &&cov&i.randomvisitp over &time ";
                              var &time   &&cov&i.randomvisitp   s&&cov&i.randomvisitp  &&cov&i.randomvisitp._diff   &&cov&i.randomvisitp._stddev 
                                   &&cov&i.randomvisitp._lb    &&cov&i.randomvisitp._ub  &&cov&i.randomvisitp._lbp &&cov&i.randomvisitp._ubp           ;               
                              run;    
                              title ;

                         %end;
                    %end;
           %end;


           %if &rungraphs = 1 %then %do;
                    %if &nsamples > 0 %then %let useboot = 1 ;
                    %else %let useboot = 0 ;
                    %construct_graphs(
                         time=&time ,
                         outcome=&outc ,
                         compevent = &compevent ,
                         outctype = &outctype,
                         covmean=&covmeandata ,
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
     %end;




     data fin;

     %*Generating reference cumulative incidence;        
     data _ref_;
     set interv&refint._all; /*change jgy*/
     %if &outctype=binsurv %then   pDref=pD;
     %else %if &outctype = cateofu %then %do ;
          %do i = 1 %to &outclev ;
		      pDref_&i = s&outc._&i ;
          %end;
     %end ;
     %else pDref=s&outc ;;
     keep _sample_ pDref: ;
     run;



     %do i=&intstart %to &numint;

          %*Outputting summary of intervention;
          %*proc means data = interv&i._all;
          %*run;

          %*Comparing intervention to reference;
                 
          %let pd = pd ;
          %if &outctype ^= binsurv and &outctype ne cateofu %then %let pd = s&outc ;

          data interv&i; 
          merge interv&i._all _ref_;
          by _sample_; 
		  %if &outctype ^= cateofu %then %do;
          	if pDref^=0 then rr=&pd /pDref;
          	if pDref^=0 then rd=&pd - pDref;
          	if rd^=. and rd>0 then nnt = 1/rd;
			else nnt = . ; * rwl added 09Sep2024 ;
		  %end ;
		  %else %do ;
			 %do ii = 1 %to &outclev ;
				if pDref_&ii^=0 then rr_&ii=s&outc._&ii /pDref_&ii;
          	    if pDref_&ii^=0 then rd_&ii=s&outc._&ii - pDref_&ii;
          	    if rd_&ii^=. and rd_&ii>0 then nnt_&ii = 1/rd_&ii;
				else nnt_&ii = . ; * rwl added 09Sep2024;
			 %end;
		  %end;
          *logrr=log(rr); /* commented out since this was the only occurrance of this variable */
          run;

 


          %*Appending intervention datasets;
          data fin; 
          set fin interv&i; 
          if _sample_=0;
          if int=. then int=&i;
          run;



          %*Calculating bootstrap mean, variance and confidence intervals;
          proc univariate data=interv&i noprint;
          where _sample_ ne 0;
		  %if &outctype ne cateofu %then %do;
          	var &pd rr rd  ; * nnt;

          	output out = temp&i
          	mean = &pd._mean RR_mean RD_mean /* NNT_mean */ 
          	std =  &pd._std  RR_std  RD_std  /* NNT_std */
           	pctlpre = &pd._  RR_     RD_    /* NNT_ */

          	pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
		  %end;
		  %else %do;
		    	var %do ii = 1 %to &outclev ;
			        s&outc._&ii rr_&ii rd_&ii /* nnt_&ii */ 
                %end; ;
          	output out = temp&i

          		mean =  %do ii = 1 %to &outclev ;
                     s&outc._&ii._mean RR_&ii._mean RD_&ii._mean /* NNT_&ii._mean  */
				  %end;
          		std =   %do ii = 1 %to &outclev ;
                    s&outc._&ii._std  RR_&ii._std  RD_&ii._std /* NNT_&ii._std */
				  %end;
          		pctlpre =  %do ii = 1 %to &outclev ;
                        s&outc._&ii._  RR_&ii._     RD_&ii._    /* NNT_&ii._ */
                     %end;
          		pctlname = llim95 ulim95  pctlpts = 2.5 97.5;

		  %end;

          run;

          data temp&i;
          set temp&i;
          int = &i;		  
          run;

          data fin;
          merge fin temp&i;
          by int;
		  %if &outctype ne cateofu %then %do;
		  	if nnt ne . then do;
			   nnt_mean = 1/rd_mean ;
			   nnt_llim95 = 1/rd_llim95 ;
			   nnt_ulim95 = 1/rd_ulim95 ;
			end;
			else do ;
			   nnt_mean = . ;
			   nnt_llim95 = . ;
			   nnt_ulim95 = . ;
			end;
			if nnt_llim95 < 0 then nnt_llim95 = . ;
			if nnt_ulim95 < 0 then nnt_ulim95 = . ;
		  %end;
		  %else %do;
		  	%do ii = 1 %to &outclev ;
				if nnt_&ii ne . then do;
					nnt_&ii._mean = 1/rd_&ii._mean ;
					nnt_&ii._llim95 = 1/rd_&ii._ulim95 ;
					nnt_&ii._ulim95 = 1/rd_&ii._llim95 ;
			    end;
				else do ;
					nnt_&ii.mean = . ;
					nnt_&ii._llim95 = . ;
					nnt_&ii.ulim95 = . ;
				end;
				if nnt_&ii._llim95 < 0 then nnt_&ii._llim95 = . ;
				if nnt_&ii._ulim95 < 0 then nnt_&ii._ulim95 = . ;

			 %end;
		   %end;
          run;
;


          %*Deleting no longer needed datasets;
          proc datasets library=work nolist; 
          delete  interv&i interv&i._all temp&i;
          quit;

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


     %*Cleaning up results to print nicely; 

     data finfin;
     set fin;      
     %rescaleround; /* RESCALE AND ROUND OFF THE OUTPUT */
	 %if &outctype ne cateofu %then %do;
     	%labels;       /* LABEL THE OUTPUT */
	 %end
     run;

	 %if &outctype = cateofu %then %do;

	 	%let pd = s&outc ;

		data finfin(keep = int int2 &outc &pd &pd._llim95 &pd._ulim95 rd rd_llim95 rd_ulim95 nnt nnt_llim95 nnt_ulim95 
	                     rr rr_llim95 rr_ulim95 &pd._mean &pd._std intervened averinterv)  ;
		set finfin;
		label &pd ="Proportion (%)"  &outc = "&outc level"
	         &pd._llim95= "Lower limit 95% CI"
			 &pd._ulim95= "Upper limit 95% CI"
			 &pd._mean="Proportion mean"
			 &pd._std= "Proportion SE"

		;
    	 
   		*label s&outc._std    = 'Bootstrap Porportion SE';
   		*label s&outc._mean   = 'Bootstrap Proportion Mean';
   		label rr        = 'Ratio of Proportions';
   		label RD        = 'Difference of Proportions';
        label rr_llim95= "Lower limit 95% CI"
			  rr_ulim95= "Upper limit 95% CI"
			  rd_llim95= "Lower limit 95% CI"
			  rd_ulim95= "Upper limit 95% CI" ;
		label NNT       = '# Needed to Treat';
   		label NNT_ulim95 = 'Upper limit 95% CI';
   		label NNT_llim95 = 'Lower limit 95% CI';

		%do ii = 1 %to &outclev ;
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
	 %end;


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

     %if &outctype=binsurv or &outctype=bineofu or &outctype = cateofu %then %do;    /*** EM edited ***/
          title4 'PREDICTED PROPORTION UNDER SEVERAL INTERVENTIONS';
     %end;
     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4  %then %do;
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

     %if &outctype=binsurv  %then %do;
	     %if %bquote(&censor) ^= %then %do ;
		     title6 "IP-weighted natural course risk= %sysevalf(&obspm) % &additional_text "; 	
		 %end;
		 %else %do;
             title6 "Observed risk= %sysevalf(&obspm) % &additional_text ";
		 %end;
     %end;    
	%else %if  &outctype=bineofu %then %do; /*** EM edited ***/
	     %if %bquote(&censor) ^= %then %do ;
		     title6 "IP-weighted natural course proportion= %sysevalf(&obspm) % &additional_text "; 	
		 %end;
		 %else %do;
             title6 "Observed proportion= %sysevalf(&obspm) % &additional_text ";
		 %end;
     %end;  
     %else %if  &outctype =  cateofu  %then %do; /*** EM edited ***/
	     %if %bquote(&censor) ^= %then %do ;
		     title6 "IP-weighted natural course proportion  &additional_text "; 
             proc print data = proportions0 noobs label ;
			 var &outc proportion ;
			 run ;
		 
	
		 %end;
		 %else %do;
             title6 "Observed proportion  &additional_text ";
			  proc print data = proportions0 noobs label ;
			 var &outc proportion ;
			 run ;
		 %end;
     %end;   
     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 %then %do;
          title6 "Observed mean= %sysevalf(&obspm) ";
     %end;
     title7 "Data= &data, Sample size= &ssize, Monte Carlo sample size= &nsimul";
     title8 "Number of bootstrap samples= &nsamples";
     title9 "Reference intervention is &refint";

	  
	     proc print data=finfin noobs label double; 
	     var int %if &outctype = cateofu %then &outc ; &pd &pd._llim95 &pd._ulim95 rr rr_llim95 rr_ulim95 &pd._mean &pd._std intervened averinterv; 
	     run;

	     proc print data=finfin noobs label double; 
	     %if &outctype=binsurv or &outctype=bineofu %then %do;
	          var int &pd &pd._llim95 &pd._ulim95 rd rd_llim95 rd_ulim95 nnt nnt_llim95 nnt_ulim95;
	     %end;
	     %else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4 or &outctype = cateofu %then %do;
	          var int   %if &outctype = cateofu %then &outc ; s&outc s&outc._llim95 s&outc._ulim95 rd rd_llim95 rd_ulim95 
                         %if &outctype = cateofu %then nnt nnt_llim95 nnt_ulim95 ; ;    
	     %end;
	     run;

	  

     %if &outctype=binsurv %then %do ;
          data ausc (keep = _sample_ int ausc&timepoints ) ;
          set &survdata (keep = int _sample_ surv0 - surv&timepoints ) ;
          array surv{*} surv1 - surv&timepoints ;
          array ausc{*} ausc1 -  ausc&timepoints ;
          ausc1 = surv0 ;  /* = 1 */
          do i = 2 to  &timepoints ;
               ausc[i]=ausc[i-1]+surv[i];   /* integral of survival curve */
          end;
          run;
          proc transpose data=ausc out=drmst (drop= _NAME_) prefix=rmst;
          var ausc&timepoints;
          id int ;
          by _sample_;
          run;

          data ausc ;
          set ausc (where = (_sample_ = 0));
          run;

          data drmst ;
          set drmst ;
          array drmst{*} drmst&intstart - drmst&numint ;
          array rmst{*} rmst&intstart - rmst&numint ;

          do i = 1 to %eval(&numint + 1 - &intstart);
               drmst[i] = rmst[i] - rmst&refint ;
          end;
          drop i ;
          run;

          proc univariate data=drmst noprint;
          where _sample_ ne 0;
          var drmst&intstart - drmst&numint ;
          output out = drmst_stat
          mean = %do myint=&intstart %to &numint ; drmst&myint._mean %end;
          std =   %do myint=&intstart %to &numint ; drmst&myint._std %end;
          pctlpre =  %do myint=&intstart %to &numint ; drmst&myint._ %end;
          pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
          run;

          proc univariate data=drmst noprint;
          where _sample_ ne 0;
          var rmst&intstart - rmst&numint ;
          output out = rmst_stat
          mean = %do myint=&intstart %to &numint ;rmst&myint._mean %end;
          std =   %do myint=&intstart %to &numint ; rmst&myint._std %end;
          pctlpre =  %do myint=&intstart %to &numint ; rmst&myint._ %end;
          pctlname = llim95 ulim95  pctlpts = 2.5 97.5;
          run;

          %let prec = 0.001 ;

          data drmst_out ;
          merge drmst (where = (_sample_ = 0))  drmst_stat rmst_stat ;
          label rmst ='Restricted mean survival time';
          label rmst_ub = 'Upper limit 95% CI';
          label rmst_lb = 'Lower limit 95% CI';

          label drmst ='Restricted mean survival time difference';
          label drmst_ub = 'Upper limit 95% CI';
          label drmst_lb = 'Lower limit 95% CI';


          %do myint=&intstart %to &numint ;
               int = &myint ;
               rmst = round(rmst&myint,&prec) ;
               rmst_mean = round(rmst&myint._mean, &prec ) ;
               rmst_std = round(rmst&myint._std,&prec) ;
               rmst_lb = round(rmst&myint._llim95,&prec) ;
               rmst_ub = round(rmst&myint._ulim95,&prec) ;

               drmst = round(drmst&myint,&prec) ;
               drmst_mean = round(drmst&myint._mean,&prec) ;
               drmst_std = round(drmst&myint._std,&prec) ;
               drmst_lb = round(drmst&myint._llim95,&prec) ;
               drmst_ub = round(drmst&myint._ulim95,&prec) ;

               output ;
          %end;

          keep int rmst rmst_mean rmst_std rmst_lb rmst_ub  drmst drmst_mean drmst_std drmst_lb drmst_ub ;
          run;


          title4 "RESTRICTED MEAN SURVIVAL TIME AFTER &timepoints TIME POINTS";
          title7 "Data= &data, Sample size= &ssize, Monte Carlo sample size= &nsimul";
          title8 "Number of bootstrap samples= &nsamples";
          title9 "Reference intervention is &refint";
          proc print data = drmst_out noobs label double ;
          var int rmst rmst_lb rmst_ub  drmst drmst_lb drmst_ub ;
          run;

     %end;

     title;

     %* Deleting no longer needed datasets;

     proc datasets library=work nolist; 
     delete    _paramdata_  _simuldata_ _inputd_ _beta_  _ref_ fin finfin
     %if &outctype = binsurv %then ausc drmst drmst_stat rmst_stat drmst_out  ;
	 %if &outctype = cateofu %then proportions0 ;
     %if &check_cov_models %then _diff_mean  _cov_std2 
     %do ii = 0 %to &ncov ; _mean_&ii %if &&usevisitp&ii = 1 %then _mean_vp&ii ; %end ;
     ;        
     ;
     quit; 

%mend results;

%macro makeknots5(howmany);
   %local lloopp ;
   %do lloopp=1 %to %eval(&howmany-1); &lloopp..5  %end; 
%mend makeknots5;

%macro bootstrap_results(
       bootlib = ,
       bootname= ,
       survdata= ,
       outc = ,
       compevent= ,
       censor = , 
       outctype = binsurv,
       combine_survdata = 0 ,
       check_cov_models = 0,
       covmeandata =  ,
       observed_surv = ,
       print_cov_means = 0 ,
       savecovmean = 0,   
       ncov = , 
       hazardratio = 0,
       bootstrap_hazard = 0 ,
       hazardname = ,
       intcomp = ,
       time = ,
       timepoints = ,
       numparts = 1,
       samplestart = 0,
       sampleend = 20,
       numboot = 20, 
       numint = 0,
       refint = 0 ,
       rungraphs = 0 ,
       graphfile = graph.pdf ,
       resultsdata = ,
       runnc = 1 ,
       usevisitp = ,
       covvisitp= ,        
       titledata=,
       title1= ,
       title2= ,
       title3= ,
       tsize = 1,
       printlogstats = 1
       );

    
     

     %let save_raw_covmean = 0 ; 
     %if &runnc = 1 %then %let intstart = 0 ;
     %else %let intstart = 1 ;

     %local datalist i j part bstart0 bstart1 bend1 betas ;
 
     /* stack intervention data sets */

     %if &numparts > 1 %then %do;    
         %let bstart1 = %scan(&samplestart,1);
        
         %let bend1 = %scan(&sampleend,1);
         %let nsample = %eval(&bend1 - &bstart1 + 1) ;
         %let datalist = &bootlib..&bootname._&bstart1._&bend1 ;          
         %do part = 2 %to &numparts ;
              %let bstart = %scan(&samplestart,%eval(&part),%str( ));
              %let bend = %scan(&sampleend,%eval(&part),%str( ));
              %let nsample = %eval(&nsample + &bend - &bstart + 1);
              %let datalist = &datalist &bootlib..&bootname._&bstart._&bend ;
          %end;
     %end;
     %else %do;
        %let datalist = &bootlib..&bootname._&samplestart._&sampleend ;
        %let nsample = %eval( &sampleend - &samplestart + 1) ;
     %end; 

     data %do int = 0 %to &numint ;
              interv&int._all 
          %end;
          ;
     set &datalist ;
     %do int = 0 %to &numint ;
         if int = &int then output interv&int._all ;
     %end;
     run;



     data _null_;
     set interv&refint._all (obs = 1 keep = n dataname ssize obsp );
     call symput('nsimul',compress(n));
     call symput('data',compress(dataname));
     call symput('ssize',compress(ssize));
     call symput('obsp',compress(obsp));
     run;


     %if &check_cov_models = 1 %then %do;

        %if &numparts > 1 %then %do;    
            %let bstart1 = %scan(&samplestart,1,%str( ));
        
            %let bend1 = %scan(&sampleend,1,%str( ));
            %let nsample = %eval(&bend1 - &bstart1 + 1) ;
            %let datalist = &bootlib..&covmeandata._&bstart1._&bend1 ;          
            %do part = 2 %to &numparts ;
                %let bstart = %scan(&samplestart,%eval(&part),%str( ));
                %let bend = %scan(&sampleend,%eval(&part),%str( ));
                %let nsample = %eval(&nsample + &bend - &bstart + 1);
                %let datalist = &datalist &bootlib..&covmeandata._&bstart._&bend ;
            %end;
        %end;
        %else %do;
          %let datalist = &bootlib..&covmeandata._&samplestart._&sampleend ;
        %end; 


        %let covmeanname = &covmeandata ; 

        data &covmeanname ;
        set &datalist ;
        by _sample_  ;
        run;


        proc contents data = &covmeanname out = _aa (keep = name varnum) noprint  ;
        run;

        proc sort data = _aa ;
        by varnum ;
        run;

        %local varcutoff ;
        %if %upcase(&outctype) = BINSURV %then %let varcutoff = 8 ;
        %else %let varcutoff = 4 ;
 


        proc sql noprint  ;
          select name into :namelist separated by ' ' from _aa (where = (varnum > &varcutoff)) ;
        quit; 
 
      
       %let cov0 = &time ;
       %let usevisitp0 = 0 ;
       %let cov0randomvisitp = ;
       %let numvisitp = %numargs(&usevisitp ) ;
  
       %if &numvisitp = 0 %then %do;
  
           %do i = 1 %to &ncov ;
                 %let word = %scan(&namelist,&i);
                 %let cov&i = &word ;
                 %put cov&i = &&cov&i ;
           %end;
      %end;
      %else %do;
 
 
         %let visitcounter = 1 ;
         %let icov = 1 ;
         %let totalvarnum = %eval(&ncov + &numvisitp + 1 ) ; /* add 1 for time */
         %do i = 2 %to &totalvarnum ;

                 %let word = %scan(&namelist,&i );
                 %if  &visitcounter <= &numvisitp  %then %let visitindex = %scan(&usevisitp,&visitcounter);               
                 %if  %eval( &icov - &visitindex) = 0   %then %do;
                       %let cov&icov = &word ;
                       %let i2 = %eval(&i + 1) ;
                       %let word2 = %scan(&namelist,&i2);
                       %let usevisitp&icov = 1;
                       %let cov&icov.randomvisitp = &word2 ;
                       %let i = &i2 ;

                       
                       %let visitcounter = %eval(&visitcounter + 1);

                %end;
                %else %do;
                    %let cov&icov = &word ;
                    %let usevisitp&icov = 0 ;
                    %let cov&icov.randomvisitp = ;                                                           
                %end;

                 %put cov&icov = &&cov&icov ;
                 %put usevisitp&icov = &&usevisitp&icov ;
                 %put cov&icov.randomvisitp = &&cov&icov.randomvisitp ;
                 %let icov = %eval(&icov + 1); 
          
       %end;   

   

   %end; 

         proc datasets library = work nolist ;
         delete _aa ;
         quit ;
         
 %end;

      %let nsamples = %eval(&nsample - 1);

      %if not(%symexist(usevisitp1) ) %then %do;
             %do i = 0 %to &ncov ;
                  %let usevisitp&i = 0 ;
             %end;
             %if %bquote(&usevisitp)^= %then %do;
                 %let numvisitp = %numargs(&usevisitp) ;
                 %do j = 1 %to &numvisitp ;
                     %let covindex = %scan(&usevisitp,&j);
                     %let covindex = %sysfunc(compress(&covindex)) ;
                     %let visitpname = %scan(&covvisitp,&j);
                     %let usevisitp&covindex = 1;
                     %let cov&covindex.randomvisitp = &visitpname ;
                 %end;
             %end;
      %end;

      %if &check_cov_models = 1 AND &savecovmean = 1 %then %do;
          proc copy in = work out = &bootlib ;
          select &covmeandata &observed_surv  ;                   
          run;          
      %end;

      %if &combine_survdata = 1 AND  %bquote(&survdata)^= %then %do;
           %if &numparts > 1 %then %do;    
               %let bstart1 = %scan(&samplestart,1);
        
               %let bend1 = %scan(&sampleend,1);
            
               %let datalist = &bootlib..&survdata._&bstart1._&bend1 ;          
            %do part = 2 %to &numparts ;
                %let bstart = %scan(&samplestart,&part);
                %let bend = %scan(&sampleend,&part);
                %let nsample = %eval(&nsample + &bend - &bstart + 1);
                %let datalist = &datalist &bootlib..&survdata._&bstart._&bend ;
            %end;
        %end;
        %else %do;
          %let datalist = &bootlib..&survdata._&samplestart._&sampleend ;
        %end; 



        data &survdata ; /* &bootlib..&survdata ; */
        set &datalist ;
        run;
             

      %end;
 

      %results ;

%mend ;
/**********************/
%macro listpred(main,switchind,start,stop,t0,t1,t2,t3,pre,prep=0);
/* listpred is in general called two ways, with main and contemp. contemp is for how variables covX appear in the models 
   for covY with X < Y . for main this will be run for all variables and may change  based on switchind = &i when using main */
 
/* want to make a change for when the mtype of cov&index = some and check if cov&index is in cov&i.addvars. If so then need to change the value 
   of covYaddvars to include the correct variables */
    %local index lev knot mycov ;

 

    %if &main=main %then %do;
        %do index = &start %to &stop;
            %if %upcase(&&cov&index.mtype) = ALL %then %do;
                 %if &&cov&index.ptype = conbin %then %do;  
                        &pre.&&cov&index..&t3  
                 %end;
                 %if &&cov&index.ptype = conqdc  %then %do;  
                       &pre.&&cov&index..&t3 &pre.&&cov&index..s&t3  
                 %end;
                 %if &&cov&index.ptype = conzqdc %then %do; 
                       &pre.z&&cov&index..&t3 &pre.&&cov&index..&t3 &pre.&&cov&index..s&t3 
                 %end;
                 %if &&cov&index.ptype = concub %then %do;  
                        &pre.&&cov&index..&t3 &pre.&&cov&index..s&t3 &pre.&&cov&index..c&t3  
                 %end;
                 %if &&cov&index.ptype = conspl  %then %do; 
                     &pre.&&cov&index..&t3
                     %do knot = 1 %to %eval(&&cov&index.lev - 2); 
                          &pre.&&cov&index..&t3._spl&knot 
                     %end; 
                 %end;
                 %if &&cov&index.ptype = concat %then %do;
                     %do lev = 1 %to %eval(&&cov&index.lev - 1); 
                         &pre.&&cov&index..&t3._&lev 
                     %end; 
                 %end;

                 %if &&cov&index.ptype = lag1bin  %then %do; 
                       &pre.&&cov&index..&t2 
                 %end;
                 %if &&cov&index.ptype = lag1qdc  %then %do; 
                      &pre.&&cov&index..&t2 &pre.&&cov&index..s&t2 
                 %end;
                 %if &&cov&index.ptype = lag1zqdc %then %do; 
                      &pre.z&&cov&index..&t2 &pre.&&cov&index..&t2 &pre.&&cov&index..s&t2 
                 %end;
                 %if &&cov&index.ptype = lag1cub  %then %do; 
                        &pre.&&cov&index..&t2 &pre.&&cov&index..s&t2 &pre.&&cov&index..c&t2 
                 %end;
                 %if &&cov&index.ptype = lag1spl  %then %do; 
                       &pre.&&cov&index..&t2
                       %do knot = 1 %to %eval(&&cov&index.lev - 2); 
                           &pre.&&cov&index..&t2._spl&knot 
                       %end; 
                 %end;
                 %if &&cov&index.ptype = lag1cat %then %do;
                     %do lev = 1 %to %eval(&&cov&index.lev - 1); 
                           &pre.&&cov&index..&t2._&lev 
                     %end; 
                 %end;
                         
                 %if &&cov&index.ptype = lag2bin  %then %do; 
                        &pre.&&cov&index..&t1 &pre.&&cov&index..&t2 
                 %end;
                 %if &&cov&index.ptype = lag2qdc  %then %do; 
                       &pre.&&cov&index..&t1 &pre.&&cov&index..&t2 &pre.&&cov&index..s&t1 &pre.&&cov&index..s&t2 
                 %end;
                 %if &&cov&index.ptype = lag2zqdc %then %do; 
                       &pre.z&&cov&index..&t1 &pre.z&&cov&index..&t2 &pre.&&cov&index..&t1 &pre.&&cov&index..&t2 &pre.&&cov&index..s&t1 &pre.&&cov&index..s&t2 
                 %end;
                 %if &&cov&index.ptype = lag2cub  %then %do; 
                         &pre.&&cov&index..&t1 &pre.&&cov&index..&t2 &pre.&&cov&index..s&t1 &pre.&&cov&index..s&t2 &pre.&&cov&index..c&t1 &pre.&&cov&index..c&t2 
                 %end;
                 %if &&cov&index.ptype = lag2spl  %then %do; 
                        &pre.&&cov&index..&t1 &pre.&&cov&index..&t2
                        %do knot = 1 %to %eval(&&cov&index.lev - 2); 
                            &pre.&&cov&index..&t1._spl&knot &pre.&&cov&index..&t2._spl&knot 
                        %end; 
                 %end;
                 %if &&cov&index.ptype = lag2cat %then %do;
                     %do lev = 1 %to %eval(&&cov&index.lev - 1); 
                          &pre.&&cov&index..&t1._&lev &pre.&&cov&index..&t2._&lev 
                     %end; 
                 %end;
                 %if &&cov&index.ptype = lag3bin  %then %do; 
                       &pre.&&cov&index..&t0 &pre.&&cov&index..&t1 &pre.&&cov&index..&t2 
                 %end;
                 %if &&cov&index.ptype = lag3qdc  %then %do; 
                         &pre.&&cov&index..&t0 &pre.&&cov&index..&t1 &pre.&&cov&index..&t2 &pre.&&cov&index..s&t0 &pre.&&cov&index..s&t1 &pre.&&cov&index..s&t2 
                 %end;
                 %if &&cov&index.ptype = lag3zqdc  %then %do; 
                       &pre.z&&cov&index..&t0 &pre.z&&cov&index..&t1 &pre.z&&cov&index..&t2 
                 %end;
                 %if &&cov&index.ptype = lag3zqdc %then %do; 
                        &pre.&&cov&index..&t0 &pre.&&cov&index..&t1 &pre.&&cov&index..&t2 &pre.&&cov&index..s&t0 &pre.&&cov&index..s&t1 &pre.&&cov&index..s&t2 
                 %end;
                 %if &&cov&index.ptype = lag3cub  %then %do;
                     &pre.&&cov&index..&t0 &pre.&&cov&index..&t1 &pre.&&cov&index..&t2
                     &pre.&&cov&index..s&t0 &pre.&&cov&index..s&t1 &pre.&&cov&index..s&t2
                     &pre.&&cov&index..c&t0 &pre.&&cov&index..c&t1 &pre.&&cov&index..c&t2 
                 %end;
                 %if &&cov&index.ptype = lag3spl  %then %do; 
                      &pre.&&cov&index..&t0 &pre.&&cov&index..&t1 &pre.&&cov&index..&t2
                      %do knot = 1 %to %eval(&&cov&index.lev - 2); 
                             &pre.&&cov&index..&t0._spl&knot &pre.&&cov&index..&t1._spl&knot &pre.&&cov&index..&t2._spl&knot 
                      %end; 
                 %end;
                 %if &&cov&index.ptype = lag3cat %then %do;
                     %do lev = 1 %to %eval(&&cov&index.lev - 1); 
                         &pre.&&cov&index..&t0._&lev &pre.&&cov&index..&t1._&lev &pre.&&cov&index..&t2._&lev 
                     %end; 
                 %end;
                 %if &&cov&index.ptype = skpbin  %then %do; 
                       &pre.&&cov&index..&t2 &pre.&&cov&index..&t2._ti 
                 %end;
                 %if &&cov&index.ptype = skpqdc  %then %do; 
                       &pre.&&cov&index..&t2 &pre.&&cov&index..&t2._ti &pre.&&cov&index..s&t2 &pre.&&cov&index..s&t2._ti 
                 %end;
                 %if &&cov&index.ptype = skpzqdc  %then %do; 
                       &pre.z&&cov&index..&t2 &pre.z&&cov&index..&t2._ti &pre.&&cov&index..&t2 &pre.&&cov&index..&t2._ti &pre.&&cov&index..s&t2 &pre.&&cov&index..s&t2._ti  
                 %end;
                 %if &&cov&index.ptype = skpcub  %then %do; 
                         &pre.&&cov&index..&t2 &pre.&&cov&index..&t2._ti &pre.&&cov&index..s&t2 &pre.&&cov&index..s&t2._ti &pre.&&cov&index..c&t2 &pre.&&cov&index..c&t2._ti  
                 %end;
                 %if &&cov&index.ptype = skpcat %then %do;
                     %do lev = 1 %to %eval(&&cov&index.lev - 1); 
                          &pre.&&cov&index..&t2._&lev &pre.&&cov&index..&t2._&lev._ti  
                     %end; 
                 %end;
                 %if &&cov&index.ptype = skpspl %then %do; 
                        &pre.&&cov&index..&t2 &pre.&&cov&index..&t2._ti
                        %do knot = 1 %to %eval(&&cov&index.lev - 2); 
                               &pre.&&cov&index..&t2._spl&knot &pre.&&cov&index..&t2._spl&knot._ti  
                        %end; 
                 %end;

                 %if &&cov&index.ptype = tsswitch1 %then %do;

                     %if &switchind = %then %do; 
                          &pre.&&cov&index..&t2 &pre.ts&&cov&index..&t2._inter   
					          %do knot = 1 %to %eval(&&cov&index.lev - 2);
                                &pre.ts&&cov&index..&t2._inter_spl&knot    
							   %end;
																				  

                     %end;
                     %else %do;
                          %if %eval(&index)>%eval(&switchind) %then %do;
                                  &pre.&&cov&index..&t2 &pre.ts&&cov&index..&t2._inter  
									%do knot = 1 %to %eval(&&cov&index.lev - 2); 
										&pre.ts&&cov&index..&t2._inter_spl&knot 
									%end; 
																					     
                          %end;
                          %else %if %eval(&index)<%eval(&switchind) %then %do;
                                  &pre.&&cov&index..&t3 &pre.ts&&cov&index..&t3._inter  
									%do knot = 1 %to %eval(&&cov&index.lev - 2); 
										&pre.ts&&cov&index..&t3._inter_spl&knot
									%end;
																						 
                          %end;
                     %end;



                  %end; /*end tsswitch1*/
                  
                   %if &&cov&index.ptype=cumavg  %then %do;
                      /* prep = 1 option needed for cumavg type variables to keep the original variable and not the cumavg variable in the data used
                        for the simulated data sets . this option is only used in the dataprep macro */
                   /**
                        %if &index ge &switchind %then %do;
                           %if &prep = 0 %then   &pre.&&cov&index.._cumavg&t2;
                           %else %if &prep = 1 %then &pre.&&cov&index..&t2 ;
                         %end;
                    ***/

				   /* when listpred is called to create initial model lists cov&index.knots is either 0 or 3,4,5. and the first part of the code is run. After the 
				      various knots have been determined cov&index.knots is a list of at least 3 knots and the second part of the code is run */

/* for new version, allow for setting of the acutal knots to use. Here covX.knots is a list with at least 3 values */ 
                                                                  
						%if &prep = 0 %then   %do;
							%if &index ge &switchind %then %do;							   
								%if %bquote(&&cov&index.knots ) ^= %then %do; 								
									%if %numargs(&&cov&index.knots ) = 1 %then %do;

										%if &&cov&index.knots = 0 %then  &pre.&&cov&index.._cumavg&t2   ;
										%if &&cov&index.knots = 3 %then  
											&pre.&&cov&index.._cumavg&t2   &pre.&&cov&index.._cumavg&t2._spl1   ;
										%if &&cov&index.knots = 4 %then                                                                     
											&pre.&&cov&index.._cumavg&t2 &pre.&&cov&index.._cumavg&t2._spl1 &pre.&&cov&index.._cumavg&t2._spl2;
										%if &&cov&index.knots = 5 %then                                                 
											&pre.&&cov&index.._cumavg&t2  &pre.&&cov&index.._cumavg&t2._spl1 &pre.&&cov&index.._cumavg&t2._spl2 &pre.&&cov&index.._cumavg&t2._spl3;
									%end;                                
									%else %do;
										&pre.&&cov&index.._cumavg&t2
										%do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
											&pre.&&cov&index.._cumavg&t2._spl&knotcount  
										%end;
									%end; 
								%end;
							%end;
							%else %if &index < &switchind %then %do;
								%if %bquote(&&cov&index..knots ) ^= %then %do; 
									%if %numargs(&&cov&index.knots ) = 1 %then %do;

										%if &&cov&index.knots = 0 %then    &pre.&&cov&index.._cumavg&t3 ;
										%if &&cov&index.knots = 3 %then   
											&pre.&&cov&index.._cumavg&t3 &pre.&&cov&index.._cumavg&t3._spl1 ;
										%if &&cov&index.knots = 4 %then 						  
											&pre.&&cov&index.._cumavg&t3 &pre.&&cov&index.._cumavg&t3._spl1 &pre.&&cov&index.._cumavg&t3._spl2;
										%if &&cov&index.knots = 5 %then 
											&pre.&&cov&index.._cumavg&t3  &pre.&&cov&index.._cumavg&t3._spl1 &pre.&&cov&index.._cumavg&t3._spl2 &pre.&&cov&index.._cumavg&t3._spl3;
									%end;                                
									%else %do;

										&pre.&&cov&index.._cumavg&t3 
										%do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
											&pre.&&cov&index.._cumavg&t3._spl&knotcount 
										%end;
									%end; 
								%end;
							%end;
						%end;                              
						%else %if &prep = 1 %then &pre.&&cov&index..&t2  ;


                   %end;

                     
                   %if &&cov&index.ptype=cumavgcat  %then %do;
                      /* prep = 1 option needed for cumavg type variables to keep the original variable and not the cumavg variable in the data used
                        for the simulated data sets . this option is only used in the dataprep macro */

                        %if &index ge &switchind %then %do;
                            %do lev = 1 %to %eval(&&cov&index.lev - 1); 
                          
                     
                                %if &prep = 0 %then   &pre.&&cov&index.._cumavg&t2._&lev ;
                                %else %if &prep = 1 %then &pre.&&cov&index..&t2._&lev  ;
                             %end;
                        %end;

                   %end;

                     %if &&cov&index.ptype=lag1cumavgcat  %then %do;
                      /* prep = 1 option needed for cumavg type variables to keep the original variable and not the cumavg variable in the data used
                        for the simulated data sets . this option is only used in the dataprep macro */

                      
                           
                          
                     
                                %if &prep = 0 %then  %do;
                                      %if &index ge &switchind %then %do;
                                            %do lev = 1 %to %eval(&&cov&index.lev - 1);
                                                &pre.&&cov&index..&t2._&lev 
                                            %end;
                                            %do lev = 1 %to %eval(&&cov&index.lev - 1);
                                                &pre.&&cov&index.._cumavg&t1._&lev
                                            %end;
                                       %end;
                                       %else %if &index < &switchind %then %do;
                                             %do lev = 1 %to %eval(&&cov&index.lev - 1) ;
                                                &pre.&&cov&index..&t3._&lev 
                                             %end;
                                              
                                             
                                              %do lev = 1 %to %eval( &&cov&index.lev - 1) ;
                                                  &pre.&&cov&index.._cumavg&t2._&lev  
                                             %end;
                                        %end;
                                 %end;
                                 %else %if &prep = 1 %then &pre.&&cov&index..&t2  ;
                            
                        

                   %end;
                   

                      %if &&cov&index.ptype=lag1cumavg  %then %do;
                      /* prep = 1 option needed for cumavg type variables to keep the original variable and not the cumavg variable in the data used
                        for the simulated data sets . this option is only used in the dataprep macro */
 
                                                                                                       
                                 %if &prep = 0 %then   %do;
                                     %if &index ge &switchind %then %do;
                                     %if %bquote(&&cov&index.knots ) ^= %then %do; 
                                         %if %numargs(&&cov&index.knots ) = 1 %then %do;
                                          
                                             %if &&cov&index.knots = 0 %then  &pre.&&cov&index..&t2 &pre.&&cov&index.._cumavg&t1 ;
                                             %if &&cov&index.knots = 3 %then  &pre.&&cov&index..&t2   &pre.&&cov&index..&t2._spl1  
                                                                            &pre.&&cov&index.._cumavg&t1 &pre.&&cov&index.._cumavg&t1._spl1 ;
                                             %if &&cov&index.knots = 4 %then 
                                                                   &pre.&&cov&index..&t2  &pre.&&cov&index..&t2._spl1 &pre.&&cov&index..&t2._spl2
                                                              &pre.&&cov&index.._cumavg&t1 &pre.&&cov&index.._cumavg&t1._spl1 &pre.&&cov&index.._cumavg&t1._spl2;
                                             %if &&cov&index.knots = 5 %then 
                                                &pre.&&cov&index..&t2  &pre.&&cov&index..&t2._spl1 &pre.&&cov&index..&t2._spl2  &pre.&&cov&index..&t2._spl3
                                                &pre.&&cov&index.._cumavg&t1  &pre.&&cov&index.._cumavg&t1._spl1 &pre.&&cov&index.._cumavg&t1._spl2 &pre.&&cov&index.._cumavg&t1._spl3;
                                         %end;                                
                                         %else %do;
                                             &pre.&&cov&index..&t2  
                                             %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                 &pre.&&cov&index..&t2._spl&knotcount   
                                             %end;
                                             &pre.&&cov&index.._cumavg&t1
                                             %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                   &pre.&&cov&index.._cumavg&t1._spl&knotcount  
                                             %end;
                                         %end; 
                                     %end;
                                     %end;
                                     %else %if &index < &switchind %then %do;
                                        %if %bquote(&&cov&index..knots ) ^= %then %do; 
                                         %if %numargs(&&cov&index.knots ) = 1 %then %do;
                                          
                                             %if &&cov&index.knots = 0 %then  &pre.&&cov&index..&t3 &pre.&&cov&index.._cumavg&t2 ;
                                             %if &&cov&index.knots = 3 %then  &pre.&&cov&index..&t3   &pre.&&cov&index..&t3._spl1  
                                                                            &pre.&&cov&index.._cumavg&t2 &pre.&&cov&index.._cumavg&t2._spl1 ;
                                             %if &&cov&index.knots = 4 %then 
                                                                   &pre.&&cov&index..&t3  &pre.&&cov&index..&t3._spl1 &pre.&&cov&index..&t3._spl2
                                                              &pre.&&cov&index.._cumavg&t2 &pre.&&cov&index.._cumavg&t2._spl1 &pre.&&cov&index.._cumavg&t2._spl2;
                                             %if &&cov&index.knots = 5 %then 
                                                &pre.&&cov&index..&t3  &pre.&&cov&index..&t3._spl1 &pre.&&cov&index..&t3._spl2  &pre.&&cov&index..&t3._spl3
                                                &pre.&&cov&index.._cumavg&t2  &pre.&&cov&index.._cumavg&t2._spl1 &pre.&&cov&index.._cumavg&t2._spl2 &pre.&&cov&index.._cumavg&t2._spl3;
                                         %end;                                
                                         %else %do;
                                            &pre.&&cov&index..&t3   
                                             %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                 &pre.&&cov&index..&t3._spl&knotcount   
                                             %end;
                                             &pre.&&cov&index.._cumavg&t2 
                                              %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                  &pre.&&cov&index.._cumavg&t2._spl&knotcount 
                                             %end;
                                         %end; 
                                     %end;
                                     %end;
                                 %end;                              
                                 %else %if &prep = 1 %then &pre.&&cov&index..&t2  ;
                             
                              
                                             
                   %end;

                    %if &&cov&index.ptype=lag2cumavgcat  %then %do;
                      /* prep = 1 option needed for cumavg type variables to keep the original variable and not the cumavg variable in the data used
                        for the simulated data sets . this option is only used in the dataprep macro */

                      
                           
                          
                     
                                %if &prep = 0 %then  %do;
                                      %if &index ge &switchind %then %do;
                                            %do lev = 1 %to %eval(&&cov&index.lev - 1);
                                                &pre.&&cov&index..&t2._&lev  
                                            %end;
                                            %do lev = 1 %to %eval(&&cov&index.lev - 1);
                                                &pre.&&cov&index..&t1._&lev 
                                             %end;
                                             %do lev = 1 %to %eval(&&cov&index.lev - 1);
                                                 &pre.&&cov&index.._cumavg&t0._&lev 

                                             %end;

                                       %end;
                                       %else %if &index < &switchind %then %do;
                                            
                                             %do lev = 1 %to %eval(&&cov&index.lev - 1);
                                                &pre.&&cov&index..&t3._&lev  
                                            %end;
                                            %do lev = 1 %to %eval(&&cov&index.lev - 1);
                                                &pre.&&cov&index..&t2._&lev 
                                             %end;
                                             %do lev = 1 %to %eval(&&cov&index.lev - 1);
                                                 &pre.&&cov&index.._cumavg&t1._&lev 

                                             %end;
                                        %end;
                                 %end;
                                 %else %if &prep = 1 %then &pre.&&cov&index..&t2  ;
                            
                        

                   %end;

                      %if &&cov&index.ptype=lag2cumavg  %then %do;
                      /* prep = 1 option needed for cumavg type variables to keep the original variable and not the cumavg variable in the data used
                        for the simulated data sets . this option is only used in the dataprep macro */
                       
                        
                               
                                 %if &prep = 0 %then   %do;
                                     %if &index ge &switchind %then %do;
                                          
                                         %if %bquote(&&cov&index.knots ) ^= %then %do; 
                                             %if %numargs(&&cov&index.knots ) = 1 %then %do;

                                                 %if &&cov&index.knots = 0 %then  &pre.&&cov&index..&t2 &pre.&&cov&index..&t1 &pre.&&cov&index.._cumavg&t0 ;
                                                 %if &&cov&index.knots = 3 %then  &pre.&&cov&index..&t2 &pre.&&cov&index..&t2._spl1  
                                                                                    &pre.&&cov&index..&t1 &pre.&&cov&index..&t1._spl1
                                                                                    &pre.&&cov&index.._cumavg&t0 &pre.&&cov&index.._cumavg&t0._spl1 ;
                                                 %if &&cov&index.knots = 4 %then 
                                                                             &pre.&&cov&index..&t2  &pre.&&cov&index..&t2._spl1 &pre.&&cov&index..&t2._spl2 
                                                                             &pre.&&cov&index..&t1  &pre.&&cov&index..&t1._spl1 &pre.&&cov&index..&t1._spl2 
                                                                             &pre.&&cov&index.._cumavg&t0 &pre.&&cov&index.._cumavg&t0._spl1 &pre.&&cov&index.._cumavg&t0._spl2;
                                                 %if &&cov&index.knots = 5 %then 
                                                         &pre.&&cov&index..&t2  &pre.&&cov&index..&t2._spl1 &pre.&&cov&index..&t2._spl2  &pre.&&cov&index..&t2._spl3 
                                                         &pre.&&cov&index..&t1  &pre.&&cov&index..&t1._spl1 &pre.&&cov&index..&t1._spl2  &pre.&&cov&index..&t1._spl3 
                                                         &pre.&&cov&index.._cumavg&t0  &pre.&&cov&index.._cumavg&t0._spl1 &pre.&&cov&index.._cumavg&t0._spl2 &pre.&&cov&index.._cumavg&t0._spl3;
                                             %end;                                
                                             %else %do;
                                                 &pre.&&cov&index..&t2
                                                 %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                     &pre.&&cov&index..&t2._spl&knotcount  
                                                 %end;
                                                 &pre.&&cov&index..&t1
                                                 %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                     &pre.&&cov&index..&t1._spl&knotcount 
                                                 %end;
                                                 &pre.&&cov&index.._cumavg&t0
                                                 %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                     &pre.&&cov&index.._cumavg&t0._spl&knotcount 
                                                 %end;
                                             %end; 
                                         %end;
                                     %end;
                                     %else %if &index < &switchind %then %do;
                                        
                                         %if %bquote(&&cov&index.knots ) ^= %then %do; 
                                             %if %numargs(&&cov&index.knots ) = 1 %then %do;

                                                 %if &&cov&index.knots = 0 %then  &pre.&&cov&index..&t3 &pre.&&cov&index..&t2 &pre.&&cov&index.._cumavg&t1 ;
                                                 %if &&cov&index.knots = 3 %then  &pre.&&cov&index..&t3 &pre.&&cov&index..&t3._spl1  
                                                                                    &pre.&&cov&index..&t2 &pre.&&cov&index..&t2._spl1
                                                                                    &pre.&&cov&index.._cumavg&t1 &pre.&&cov&index.._cumavg&t1._spl1 ;
                                                 %if &&cov&index.knots = 4 %then 
                                                                         &pre.&&cov&index..&t3  &pre.&&cov&index..&t3._spl1 &pre.&&cov&index..&t3._spl2 
                                                                         &pre.&&cov&index..&t2  &pre.&&cov&index..&t2._spl1 &pre.&&cov&index..&t2._spl2 
                                                                         &pre.&&cov&index.._cumavg&t1 &pre.&&cov&index.._cumavg&t1._spl1 &pre.&&cov&index.._cumavg&t1._spl2;
                                                 %if &&cov&index.knots = 5 %then 
                                                         &pre.&&cov&index..&t3  &pre.&&cov&index..&t3._spl1 &pre.&&cov&index..&t3._spl2  &pre.&&cov&index..&t3._spl3 
                                                         &pre.&&cov&index..&t2  &pre.&&cov&index..&t2._spl1 &pre.&&cov&index..&t2._spl2  &pre.&&cov&index..&t2._spl3 
                                                         &pre.&&cov&index.._cumavg&t1  &pre.&&cov&index.._cumavg&t1._spl1 &pre.&&cov&index.._cumavg&t1._spl2 &pre.&&cov&index.._cumavg&t1._spl3;
                                             %end;                                
                                             %else %do;
                                                 &pre.&&cov&index..&t3
                                                 %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                     &pre.&&cov&index..&t3._spl&knotcount  
                                                 %end;
                                                 &pre.&&cov&index..&t2
                                                 %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                     &pre.&&cov&index..&t2._spl&knotcount  
                                                 %end;
                                                 &pre.&&cov&index.._cumavg&t1
                                                 %do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
                                                     &pre.&&cov&index.._cumavg&t1._spl&knotcount 
                                                 %end;
                                             %end; 
                                         %end;
                                     %end;           
                                 %end;                              
                                 %else %if &prep = 1 %then &pre.&&cov&index..&t2 ;
                             
                             
                   %end;
                   %if &&cov&index.ptype =rcumavg %then %do ;
                           /* prep = 1 option needed for cumavg type variables to keep the original variable and not the cumavg variable in the data used
                        for the simulated data sets this option is only used in the dataprep macro */
                          %if &index ge &switchind %then %do ;
                              %if &prep = 0 %then &pre.&&cov&index.._rcumavg&t2;
                              %else %if &prep = 1 %then &pre.&&cov&index..&t2 ;
                           %end;
                   %end;

                  /**** original code                  
                  %if &&usevisitp&index = 1  &    ( (&switchind = &index) OR &switchind = )  %then %do;      
                    &pre.ts_last_&&cov&index..&t2
                  %end;
                   *****/
                    %if &&usevisitp&index = 1   %then %do; 
                        /* this combines the main and contemp options for including the ts_last portion of the visit 
                           process into a single call using main = main and the switchind variable */

                        %let mycov = ;
                        %if &switchind = %then %let mycov= ts_last_&&cov&index ;
                        %else %if &index >= &switchind  %then %let mycov = ts_last_&&cov&index.._l1 ; 
                        %else %if &index < &switchind %then %let mycov = ts_last_&&cov&index ;
                     
                        &mycov                       
                  %end;
             %end;
          %end ; /* for mtype = all */
    %end;

    /* contemp is using start = 0 and stop = &i - 1. In this case cov&index is used for variable after &index so lags need to 
       be changed, lag1 is at the same time as the &i variable */

    %if &main=contemp %then %do;
        %do index = &start %to &stop;
            %if %upcase(&&cov&index.mtype) = ALL %then %do;
                 %if    &&cov&index.ptype = lag1bin
                     or &&cov&index.ptype = lag2bin
                     or &&cov&index.ptype = lag3bin %then %do; 
                          &pre.&&cov&index..&t3  
                 %end;
                 %if   &&cov&index.ptype = lag1qdc
                    or &&cov&index.ptype = lag2qdc
                    or &&cov&index.ptype = lag3qdc %then %do; 
                         &pre.&&cov&index..&t3 &pre.&&cov&index..s&t3  
                 %end;
                 %if &&cov&index.ptype = lag1zqdc
                     or &&cov&index.ptype = lag2zqdc
                     or &&cov&index.ptype = lag3zqdc
                     %then %do; 
                           &pre.z&&cov&index..&t3 &pre.&&cov&index..&t3 &pre.&&cov&index..s&t3 
                 %end;
                 %if   &&cov&index.ptype = lag1spl
                    or &&cov&index.ptype = lag2spl
                    or &&cov&index.ptype = lag3spl  %then %do; 
                        &pre.&&cov&index..&t3
                        %do knot = 1 %to %eval(&&cov&index.lev - 2); 
                             &pre.&&cov&index..&t3._spl&knot 
                        %end; 
                 %end;
                 %if &&cov&index.ptype = lag1cat
                    or &&cov&index.ptype = lag2cat
                    or &&cov&index.ptype = lag3cat %then %do;
                         %do lev = 1 %to %eval(&&cov&index.lev - 1); 
                           &pre.&&cov&index..&t3._&lev 
                         %end; 
                 %end;
                 %if &&cov&index.ptype = skpbin %then %do; 
                        &pre.&&cov&index..&t3 &pre.&&cov&index..&t3._ti  
                 %end;
                 %if &&cov&index.ptype = skpqdc  %then %do; 
                       &pre.&&cov&index..&t3 &pre.&&cov&index..&t3._ti &pre.&&cov&index..s&t3 &pre.&&cov&index..s&t3._ti 
                 %end;
                 %if &&cov&index.ptype = skpzqdc  %then %do; 
                        &pre.z&&cov&index..&t3 &pre.z&&cov&index..&t3._ti &pre.&&cov&index..&t3 &pre.&&cov&index..&t3._ti 
                        &pre.&&cov&index..s&t3 &pre.&&cov&index..s&t3._ti  
                  %end;
                 %if &&cov&index.ptype = skpcub  %then %do; 
                          &pre.&&cov&index..&t3 &pre.&&cov&index..&t3._ti &pre.&&cov&index..s&t3 &pre.&&cov&index..s&t3._ti 
                          &pre.&&cov&index..c&t3 &pre.&&cov&index..c&t3._ti  
                 %end;
                 %if &&cov&index.ptype = skpcat %then %do;
                     %do lev = 1 %to %eval(&&cov&index.lev - 1); 
                          &pre.&&cov&index..&t3._&lev &pre.&&cov&index..&t3._&lev._ti 
                     %end; 
                 %end;
                 %if &&cov&index.ptype = skpspl %then %do; 
                     &pre.&&cov&index..&t3 &pre.&&cov&index..&t3._ti
                     %do knot = 1 %to %eval(&&cov&index.lev - 2); 
                           &pre.&&cov&index..&t3._spl&knot &pre.&&cov&index..&t3._spl&knot._ti  
                     %end; 
                 %end; 
               
                /** 
                 %if &&usevisitp&index = 1   %then ts_last_&&cov&index ;
                **/
                 %if &&cov&index.ptype=cumavg  %then %do;
                        /* &&cov&index.._cumavg  */ /* is this taken into account in the main option like for the other cumavg ptypes ?? */       
                  %end;

                  %if &&cov&index.ptype = cumavgcat %then %do ;
                        %do lev = 1 %to %eval(&&cov&index.lev - 1); 
                            &&cov&index.._cumavg_&lev 
                        %end; 
                  %end;

                  %if &&cov&index.ptype=lag1cumavg  %then %do;
                      /* nothing here, taken into account in with the main option */
                   %end;
                   %if &&cov&index.ptype=lag2cumavg  %then %do;
                      /* nothing here, taken into account with the main option */
                        
                   %end;

                  %if &&cov&index.ptype=rcumavg %then %do; &pre.&&cov&index.._rcumavg&t3  %end;
             %end;
         %end;
   %end;
/* code for including complete history in eof-type outcome models */
	%local timeindex  startindex stopindex;
	%if &main=eof %then %do;
        %do index = 1 %to &ncov;
		    %if &&cov&index.etype_part2 = &timepoints %then %let  startindex = 0 ;
			%else %if &&cov&index.etype_part2 = 0 %then %let startindex = %eval(&timepoints);
			%else %let startindex = %eval(&timepoints - &&cov&index.etype_part2 ) ;

			%let stopindex = %eval(&timepoints - 1);
			%if &&cov&index.etype = cumavg or &&cov&index.etype = cumavgcat or &&cov&index.etype = cumsum or &&cov&index.etype = cumsumcat  %then %do;
	     		%if &&cov&index.etype_part2 = &timepoints %then %do;
					%let startindex = 0;
            		%let stopindex =  0;
			 	%end;
		 		%else %do;
		 			%let startindex = %eval(&timepoints - &&cov&index.etype_part2 ) ;
            		%let stopindex = %eval(&timepoints - &&cov&index.etype_part2 ) ;
		 		%end;
			%end;	 
			%if &&cov&index.etype = cumavgnew or &&cov&index.etype = cumavgcatnew or &&cov&index.etype = cumsumnew or &&cov&index.etype = cumsumcatnew %then %do;
	     		%if &&cov&index.etype_part2 = &timepoints %then %do;
					%let startindex = 0;
            		%let stopindex =  0;
			 	%end;
		 		%else %do;
		 			%let startindex = %eval(&timepoints - &&cov&index.etype_part2 - 1 ) ;
            		%let stopindex  = %eval(&timepoints - &&cov&index.etype_part2 - 1 ) ;
		 		%end;
			%end;	

		    %do timeindex = &startindex %to &stopindex ;
				%if %upcase(&&cov&index.mtype) = ALL %then %do;
					%if &&cov&index.etype = bin %then %do;  
                        &&cov&index.._&timeindex._eof  
					%end;
					%if &&cov&index.etype = qdc  %then %do;  
                       &&cov&index.._&timeindex._eof &&cov&index.._&timeindex._eofs  
					%end;
					%if &&cov&index.etype = zqdc %then %do; 
                       z&&cov&index.._&timeindex._eof &&cov&index.._&timeindex._eof &&cov&index.._&timeindex._eofs 
					%end;
					%if &&cov&index.etype = cub %then %do;  
                        &&cov&index.._&timeindex._eof &&cov&index.._&timeindex._eofs &&cov&index.._&timeindex._eofc  
					%end;
					%if &&cov&index.etype = spl  %then %do; 
						&&cov&index.._&timeindex._eof
						%do knot = 1 %to %eval(&&cov&index.elev - 2); 
                          &&cov&index.._&timeindex._eof_spl&knot 
						%end; 
					%end;
					%if &&cov&index.etype = cat %then %do;
						%do lev = 1 %to %eval(&&cov&index.elev  - 1); 
							&&cov&index.._&timeindex._eof_&lev 
						%end; 
					%end;

                 
					%if &&cov&index.etype = skpbin  %then %do; 
					   %if not( &timeindex in &&cov&index.skip ) %then &&cov&index.._&timeindex._eof /* &&cov&index.._&timeindex._eof_ti */ ;
					%end;
					%if &&cov&index.etype = skpqdc  %then %do; 
                        %if not( &timeindex in &&cov&index.skip ) %then &&cov&index.._&timeindex._eof  &&cov&index.._&timeindex._eofs ;
					%end;
					%if &&cov&index.etype = skpzqdc  %then %do; 
                      %if not( &timeindex in &&cov&index.skip ) %then  z&&cov&index.._&timeindex._eof  &&cov&index.._&timeindex._eof &&cov&index.._&timeindex._eofs ;
					%end;
					%if &&cov&index.etype = skpcub  %then %do; 
                         %if not( &timeindex in &&cov&index.skip ) %then &&cov&index.._&timeindex._eof  &&cov&index.._&timeindex._eofs   &&cov&index.._&timeindex._eofc ;  
					%end;
					%if &&cov&index.etype = skpcat %then %do;
						%do lev = 1 %to %eval(&&cov&index.elev - 1); 
							%if not( &timeindex in &&cov&index.skip ) %then  &&cov&index.._&timeindex._eof_&lev  ; 
						%end; 
					%end;
					%if &&cov&index.etype = skpspl %then %do; 
                        %if not( &timeindex in &&cov&index.skip ) %then  &&cov&index.._&timeindex._eof ; 
                        %do knot = 1 %to %eval(&&cov&index.elev - 2); 
                              %if not( &timeindex in &&cov&index.skip ) %then  &&cov&index.._&timeindex._eof_spl&knot  ;
                        %end; 
					%end;

/*******
GENERAL FORM FOR SPLINE VARIABLES WHEN COV&I..KNOTS HAS ONE ENTRY OR THREE ENTRIES  
							%if %bquote(&&cov&index.knots ) ^= %then %do; 								
									%if %numargs(&&cov&index.knots ) = 1 %then %do;

										%if &&cov&index.knots = 0 %then  &&cov&index.._cumavg   ;
										%if &&cov&index.knots = 3 %then  
											&&cov&index.._cumavg   &&cov&index.._cumavg_spl1   ;
										%if &&cov&index.knots = 4 %then                                                                     
											 &&cov&index.._cumavg &&cov&index.._cumavg_spl1 &&cov&index.._cumavg_spl2;
										%if &&cov&index.knots = 5 %then                                                 
											&&cov&index.._cumavg  &&cov&index.._cumavg_spl1 &&cov&index.._cumavg_spl2 &&cov&index.._cumavg_spl3;
									%end;                                
									%else %do;
										 &&cov&index.._cumavg&t1
										%do knotcount = 1 %to %eval(%numargs(&&cov&index.knots ) - 2) ;
											&&cov&index.._cumavg_spl&knotcount  
										%end;
									%end; 

                              %end;

********/

                  
					%if &&cov&index.etype=cumsum  %then %do; 						   								
									%if %numargs(&&cov&index.eknots ) = 1 %then %do; /* {0,3,4,5} */

										%if &&cov&index.eknots = 0 %then  &&cov&index.._cumsum_&timeindex._eof   ;
										%if &&cov&index.eknots = 3 %then  
											&&cov&index.._cumsum_&timeindex._eof   &&cov&index.._cumsum_&timeindex._eof_spl1   ;
										%if &&cov&index.eknots = 4 %then                                                                     
											 &&cov&index.._cumsum_&timeindex._eof &&cov&index.._cumsum_&timeindex._eof_spl1 &&cov&index.._cumsum_&timeindex._eof_spl2;
										%if &&cov&index.eknots = 5 %then                                                 
											&&cov&index.._cumsum_&timeindex._eof  &&cov&index.._cumsum_&timeindex._eof_spl1 &&cov&index.._cumsum_&timeindex._eof_spl2 &&cov&index.._cumsum_&timeindex._eof_spl3;
									%end;                                
									%else %do; /* list of 3 or more knots */
										 &&cov&index.._cumsum_&timeindex._eof
										%do knotcount = 1 %to %eval(%numargs(&&cov&index.eknots ) - 2) ;
											&&cov&index.._cumsum_&timeindex._eof_spl&knotcount  
										%end;
									%end;                             						                           				
					%end;
					%if &&cov&index.etype=cumsumnew  %then %do;   

						 								
									%if %numargs(&&cov&index.eknots ) = 1 %then %do;

										%if &&cov&index.eknots = 0 %then  &&cov&index.._cumsum_&timeindex._eof   ;
										%if &&cov&index.eknots = 3 %then  
											&&cov&index.._cumsum_&timeindex._eof   &&cov&index.._cumsum_&timeindex._eof_spl1   ;
										%if &&cov&index.eknots = 4 %then                                                                     
											 &&cov&index.._cumsum_&timeindex._eof &&cov&index.._cumsum_&timeindex._eof_spl1 &&cov&index.._cumsum_&timeindex._eof_spl2;
										%if &&cov&index.eknots = 5 %then                                                 
											&&cov&index.._cumsum_&timeindex._eof  &&cov&index.._cumsum_&timeindex._eof_spl1 &&cov&index.._cumsum_&timeindex._eof_spl2 &&cov&index.._cumsum_&timeindex._eof_spl3;

                                        %do iii = %eval(&timeindex + 1) %to %eval(&timepoints - 1) ;
									   		&&cov&index.._&iii._eof  
										%end;
									%end;                                
									%else %do;
									    &&cov&index.._cumsum_&timeindex._eof
										%do knotcount = 1 %to %eval(%numargs(&&cov&index.eknots ) - 2) ;
											&&cov&index.._cumsum_&timeindex._eof_spl&knotcount  
										%end;
									 	%do iii = %eval(&timeindex + 1) %to %eval(&timepoints - 1) ;
									   		&&cov&index.._&iii._eof  
										%end;
									%end;                            	
					%end;

					%if &&cov&index.etype=cumsumcat  %then %do;
                     
                            %do lev = 1 %to %eval(&&cov&index.elev - 1);                                                
                                 &&cov&index.._cumsum_&timeindex._eof_&lev                                  
                             %end;                        
					%end;

					%if &&cov&index.etype = cumsumcatnew %then %do;
					 						 						   							 
					 		%do lev = 1 %to %eval(&&cov&index.elev - 1);                                                
                                 &&cov&index.._cumsum_&timeindex._eof_&lev                                  
                             %end; 
					 
						    %do iii = %eval(&timeindex + 1) %to %eval(&timepoints - 1) ;
							   	&&cov&index.._&iii._eof  
							%end;							 							 						
					 
					%end;   

					%if &&cov&index.etype=cumavg  %then %do;  	
									%if %numargs(&&cov&index.eknots ) = 1 %then %do;

										%if &&cov&index.eknots = 0 %then  &&cov&index.._cumavg_&timeindex._eof   ;
										%if &&cov&index.eknots = 3 %then  
											&&cov&index.._cumavg_&timeindex._eof   &&cov&index.._cumavg_&timeindex._eof_spl1   ;
										%if &&cov&index.eknots = 4 %then                                                                     
											 &&cov&index.._cumavg_&timeindex._eof &&cov&index.._cumavg_&timeindex._eof_spl1 &&cov&index.._cumavg_&timeindex._eof_spl2;
										%if &&cov&index.eknots = 5 %then                                                 
											&&cov&index.._cumavg_&timeindex._eof  &&cov&index.._cumavg_&timeindex._eof_spl1 &&cov&index.._cumavg_&timeindex._eof_spl2 &&cov&index.._cumavg_&timeindex._eof_spl3;
									%end;                                
									%else %do;
										 &&cov&index.._cumavg_&timeindex._eof
										%do knotcount = 1 %to %eval(%numargs(&&cov&index.eknots ) - 2) ;
											&&cov&index.._cumavg_&timeindex._eof_spl&knotcount  
										%end;
									%end; 
					%end;
					%if &&cov&index.etype=cumavgnew  %then %do;                                                                                       												   
					 								
									%if %numargs(&&cov&index.eknots ) = 1 %then %do;

										%if &&cov&index.eknots = 0 %then  &&cov&index.._cumavg_&timeindex._eof   ;
										%if &&cov&index.eknots = 3 %then  
											&&cov&index.._cumavg_&timeindex._eof   &&cov&index.._cumavg_&timeindex._eof_spl1   ;
										%if &&cov&index.eknots = 4 %then                                                                     
											 &&cov&index.._cumavg_&timeindex._eof &&cov&index.._cumavg_&timeindex._eof_spl1 &&cov&index.._cumavg_&timeindex._eof_spl2;
										%if &&cov&index.eknots = 5 %then                                                 
											&&cov&index.._cumavg_&timeindex._eof  &&cov&index.._cumavg_&timeindex._eof_spl1 &&cov&index.._cumavg_&timeindex._eof_spl2 &&cov&index.._cumavg_&timeindex._eof_spl3;
                                        %do iii = %eval(&timeindex + 1) %to %eval(&timepoints - 1) ;
									   		&&cov&index.._&iii._eof  
										%end;
									%end;                                
									%else %do;
									    &&cov&index.._cumavg_&timeindex._eof
										%do knotcount = 1 %to %eval(%numargs(&&cov&index.eknots ) - 2) ;
											&&cov&index.._cumavg_&timeindex._eof_spl&knotcount  
										%end;
									 	%do iii = %eval(&timeindex + 1) %to %eval(&timepoints - 1) ;
									   		&&cov&index.._&iii._eof  
										%end;
									%end;                             																                         				
					%end;

                     
					%if &&cov&index.etype=cumavgcat  %then %do;
                     
                            %do lev = 1 %to %eval(&&cov&index.elev - 1);                                                
                                 &&cov&index.._cumavg_&timeindex._eof_&lev                                  
                             %end;                        
					%end;

					%if &&cov&index.etype = cumavgcatnew %then %do;
					 						 						   							 
					 		%do lev = 1 %to %eval(&&cov&index.elev - 1);                                                
                                 &&cov&index.._cumavg_&timeindex._eof_&lev                                  
                             %end; 
					 
						    %do iii = %eval(&timeindex + 1) %to %eval(&timepoints - 1) ;
							   	&&cov&index.._&iii._eof  
							%end;							 							 						
					 
					%end;                                    
				
                  /*****************
                    %if &&usevisitp&index = 1   %then %do; 
                         this combines the main and temp options for including the ts_last portion of the visit 
                           process into a single call using main = main and the switchind variable ;

                        %let mycov = ;
                        %if &switchind = %then %let mycov= ts_last_&&cov&index ;
                        %else %if &index >= &switchind  %then %let mycov = ts_last_&&cov&index.._l1 ; 
                        %else %if &index < &switchind %then %let mycov = ts_last_&&cov&index ;
                     
                        &mycov                       
                    %end;
				   *******************/
                %end ; /* for mtype = all */
			%end; /* end timeindex */
		%end ;/* end of ncov */	
	%end ; /* main = eof */	

 
 
 
 
 

%mend listpred;    


/********************************************************************************************************************************************************/

/*********************************************************************/
/** THE MACROS IN THIS FILE HANDLE INTERACTIONS BETWEEN DICHOTOMOUS **/
/** COVARIATES, OR ONE CONTINUOUS AND ONE DICHOTOMOUS. THEY CAN NOW **/
/** HANDLE FIXED BASELINE COVARIATES IN THE INTERACTION TERMS, PLUS **/
/** THANKS TO ROGER, THEY CAN NOW DEAL WITH CATEGORICAL VARIABLES.  **/
/*********************************************************************/

/* NOTE: IF THE INTERACTION TERM IS BETWEEN TWO BASELINE VARIABLES, 
   THEN IT SHOULD BE TREATED AS A FIXED VARIABLE TOO. THIS STILL 
   NEEDS TO BE ADJUSTED. */
/* ALSO, THE PRODUCT OF A DICHOTOMOUS AND A CONTINUOUS VARIABLE SHOULD
   PROBABLY INHERIT THE PREDICTOR-TYPE OF THE CONTINUOUS VARIABLE. AT
   PRESENT IT DOESN'T. */

/*FIRST A MINI-MACRO FUNCTION BY ROGER LOGAN 
FOR COUNTING THE NUMBER OF SPACE-SEPARATED WORDS IN A STRING*/
%macro numargs(arg);
 %local n word ;
   %let n = 1;
   %if %bquote(&arg)^= %then %do;
      %do %until (%qscan(&arg,%eval(&n),%str( ))=%str());
         %let word = %qscan(&arg,&n);
         %let n = %eval(&n+1);
         %end;
      %end;
   %eval(&n-1) /* there is no ; here since it will be used as %let a = %numargs(&b) ;
          and the ; is included at the end of this line  */
%mend numargs;
/*THIS MACRO GOES THROUGH THE LIST OF INTERACTION TERMS PROVIDED BY 
USER AND CREATES THE INTERACTION TERM VARIABLES.
--but only for predictors of outc, compevent and censl, 
not the time-varying covariates, which are handled below in %interactionsb*/

%macro interactions(vrbl,indexer,mytype ,fortype);
  
      %local baslineind first second firstind firstvarcat firstvarspl firstvarqdc firstvarcub firstvarzqdc firstvartss1 firstvarcum0 firstvarcum1 firstvarcum2 
             firstvarcum3 firstvarcon firstvarbin secondind secondvarcat secondvarspl secondvarqdc secondvarcub secondvarzqdc secondvartss1 secondvarcum0 secondvarcum1 
             secondvarcum2 secondvarcum3 secondvarcon secondvarbin ;

    %let &vrbl.interact=%qcmpres(&&&vrbl.interact);



    /* first count number of interaction terms for modeling this variable */
    %if %bquote(&vrbl.interact)^=%str() %then %do;
        %let &vrbl.ninterx=%eval(%numargs(&&&vrbl.interact));

        /* for each 2-way interaction term, split into first and second member of pair */
        %do iota=1 %to %eval(&&&vrbl.ninterx);
            %local interx&iota ;
            %let interx&iota=%qscan(&&&vrbl.interact,&iota,%str( ));
            %let first=%scan(&&interx&iota,1,*);
            %let second=%scan(&&interx&iota,2,*);             
            %let baselineind = 0;
                         
            %if %index(%upcase(&first),B)=1 %then %do;
                %let firstind=%substr(&first,2);
                %let firstvar=%scan(&fixedcov,&firstind);
                %let firstvarcat = 0 ;
                %let firstvarspl = 0 ;
                %let firstvarqdc = 0 ;
                %let firstvarcub = 0 ;
                %let firstvarzqdc = 0 ;
                %let firstvartss1 = 0 ;
                %let firstvarcum0 = 0 ;
                %let firstvarcum1 = 0 ;
                %let firstvarcum2 = 0 ;
                %let firstvarcum3 = 0 ;
                
                %let firstvarcon = 0 ;
                %let firstvarbin = 0 ;
                %let baselineind = 1;
            %end;
            %else %if %index(%upcase(&first),L)=1 %then %do; /* 2-25-2015, rwl create new interaction for internal indicator function */
                %let firstind = %substr(&first,2);
                %let firstvar = _ind&firstind ; 
                %let firstvarcat = 0 ;
                %let firstvarspl = 0 ;
                %let firstvarqdc = 0 ;
                %let firstvarcub = 0 ;
                %let firstvarzqdc = 0 ;
                %let firstvartss1 = 0 ;
                %let firstvarcum0 = 0 ;
                %let firstvarcum1 = 0 ;
                %let firstvarcum2 = 0 ;
                %let firstvarcum3 = 0 ;
                %let firstvarcon = 0 ;
                %let firstvarbin = 0 ;
                %let baselineind = 1;
            %end ;
            %else %do;
                %let firstvar=&&cov&first;                    
                %let firstind = &first ;
                %let firstvarcat = 0 ;
                %let firstvarspl = 0 ;
                %let firstvarqdc = 0 ;
                %let firstvarcub = 0 ;
                %let firstvarzqdc = 0 ;
                %let firstvartss1 = 0 ;
                %let firstvarcum0 = 0 ;
                %let firstvarcum1 = 0 ;
                %let firstvarcum2 = 0 ;
                %let firstvarcum3 = 0 ;                
                %let firstvarbin = 0 ;

                %if  &&cov&first.ptype=lag1bin or &&cov&first.ptype=lag2bin or
                &&cov&first.ptype=lag3bin      %then   %let firstvarbin=1;


                %else %if( &&cov&first.ptype=concat or &&cov&first.ptype=lag1cat or
                &&cov&first.ptype=lag2cat or &&cov&first.ptype=lag3cat or
                &&cov&first.ptype=skpcat)   %then   %let firstvarcat=1;

                %else  %if( &&cov&first.ptype=conspl or &&cov&first.ptype=lag1spl or
                &&cov&first.ptype=lag2spl or &&cov&first.ptype=lag3spl or
                &&cov&first.ptype=skpspl)    %then   %let firstvarspl=1;

                %else %if  &&cov&first.ptype=conqdc or &&cov&first.ptype=lag1qdc or
                &&cov&first.ptype=lag2qdc or &&cov&first.ptype=lag3qdc or
                &&cov&first.ptype=skpqdc    %then  %let firstvarqdc=1;

                %else %if  &&cov&first.ptype=concub or &&cov&first.ptype=lag1cub or
                &&cov&first.ptype=lag2cub or &&cov&first.ptype=lag3cub or
                &&cov&first.ptype=skpcub    %then  %let firstvarcub=1;

                %else %if  &&cov&first.ptype=conzqdc or &&cov&first.ptype=lag1zqdc or
                &&cov&first.ptype=lag2zqdc or &&cov&first.ptype=lag3zqdc or
                &&cov&first.ptype=skpzqdc    %then  %let firstvarzqdc=1;

                %else %if &&cov&first.ptype = tsswtich1  
                %then %let firstvartss1 = 1;
                   
                %else %if &&cov&first.ptype = cumavg   %then %let firstvarcum0 = 1 ;

                %else %if &&cov&first.ptype = lag1cumavg  %then %let firstvarcum1 = 1 ;

                %else %if &&cov&first.ptype = lag2cumavg    %then %let firstvarcum2 = 1;

                %else %if &&cov&first.ptype = rcumavg    %then %let firstvarcum3 = 1;
            %end;            

            %if %index(%upcase(&second),B)=1 %then %do;
                %let secondind=%substr(&second,2);
                %let secondvar=%scan(&fixedcov,&secondind);
                %let secondvarcat = 0 ;
                %let secondvarspl = 0 ;
                %let secondvarqdc = 0 ;
                %let secondvarcub = 0 ;
                %let secondvarzqdc = 0 ;
                %let secondvartss1 = 0 ;
                %let secondvarcum0 = 0 ;
                %let secondvarcum1 = 0 ;
                %let secondvarcum2 = 0 ;
                %let secondvarcum3 = 0 ;
                %let secondvarbin = 0 ;
                %let baselineind = 1;
            %end;
            %else %if %index(%upcase(&second),L)=1 %then %do; /* create new interaction for internal indicator function */
                %let secondind = %substr(&second,2);
                %let secondvar = _ind&secondind ; 
                %let secondvarcat = 0 ;
                %let secondvarspl = 0 ;
                %let secondvarqdc = 0 ;
                %let secondvarcub = 0 ;
                %let secondvarzqdc = 0 ;
                %let secondvartss1 = 0 ;
                %let secondvarcum0 = 0 ;
                %let secondvarcum1 = 0 ;
                %let secondvarcum2 = 0 ;
                %let secondvarcum3 = 0 ;         
                %let secondvarbin = 0 ;
                %let baselineind = 1;
            %end ;
            %else %do;
                %let secondvar=&&cov&second;                
                %let secondind = &second ;
                %let secondvarcat = 0;
                %let secondvarspl = 0 ;
                %let secondvarqdc = 0 ;
                %let secondvarcub = 0 ;
                %let secondvarzqdc = 0 ;
                %let secondvartss1 = 0;
                %let secondvarcum0 = 0 ;
                %let secondvarcum1 = 0 ;
                %let secondvarcum2 = 0 ;
                %let secondvarcum3 = 0 ;                
                %let secondvarbin = 0 ;


                %if &&cov&second.ptype=lag1bin or &&cov&second.ptype=lag2bin or
                &&cov&second.ptype=lag3bin   %then   %let secondvarbin=1;


                %else %if (&&cov&second.ptype=concat or &&cov&second.ptype=lag1cat or
                &&cov&second.ptype=lag2cat or &&cov&second.ptype=lag3cat or
                &&cov&second.ptype=skpcat )   %then   %let secondvarcat=1;
                %else %if (&&cov&second.ptype=conspl or &&cov&second.ptype=lag1spl or
                &&cov&second.ptype=lag2spl or &&cov&second.ptype=lag3spl or
                &&cov&second.ptype=skpspl )  %then     %let secondvarspl=1;

                %else %if  &&cov&second.ptype=conqdc or &&cov&second.ptype=lag1qdc or
                &&cov&second.ptype=lag2qdc or &&cov&second.ptype=lag3qdc or
                &&cov&second.ptype=skpqdc   %then  %let secondvarqdc=1;
                %else %if  &&cov&second.ptype=concub or &&cov&second.ptype=lag1cub or
                &&cov&second.ptype=lag2cub or &&cov&second.ptype=lag3cub or
                &&cov&second.ptype=skpcub    %then %let secondvarcub=1;


                %else %if  &&cov&second.ptype=conzqdc or &&cov&second.ptype=lag1zqdc or
                &&cov&second.ptype=lag2zqdc or &&cov&second.ptype=lag3zqdc or
                &&cov&second.ptype=skpzqdc    %then %let secondvarzqdc=1;

                %else %if &&cov&second.ptype = tsswitch1   %then %let secondvartss1 = 1 ;                   
                %else %if &&cov&second.ptype = cumavg  %then %let secondvarcum0 = 1 ;
                %else %if &&cov&second.ptype = lag1cumavg    %then %let secondvarcum1 = 1 ;
                %else %if &&cov&second.ptype = lag2cumavg  %then %let secondvarcum2 = 1;
                %else %if &&cov&second.ptype = rcumavg   %then %let secondvarcum3 = 1;

            %end;
         
            /* the sum can only be 0 or 1 */
            %local firsttype secondtype type1 type2 ;
            %let firsttype = %eval(&firstvarbin + &firstvarcat + &firstvarspl + &firstvarqdc + &firstvarcub +&firstvarzqdc + &firstvartss1 + &firstvarcum0 +  &firstvarcum1 + &firstvarcum2 + &firstvarcum3) ;
            %let secondtype = %eval(&secondvarbin + &secondvarcat + &secondvarspl + &secondvarqdc + &secondvarcub + &secondvarzqdc + &secondvartss1 +  &secondvarcum0 +  &secondvarcum1 + &secondvarcum2 + &secondvarcum3) ;
            
            %if       &firstvarbin = 1 %then %let type1 = 0 ;
            %else %if &firstvarcat = 1 %then %let type1 = 1 ;
            %else %if &firstvarspl = 1 %then %let type1 = 2 ;
            %else %if &firstvarqdc = 1 %then %let type1 = 3 ;
            %else %if &firstvarcub = 1 %then %let type1 = 4 ;
            %else %if &firstvarzqdc = 1 %then %let type1 = 5 ;
            %else %if &firstvartss1 = 1 %then %let type1 = 6 ;
            %else %if &firstvarcum0 = 1 %then %let type1 = 7 ;
            %else %if &firstvarcum1 = 1 %then %let type1 = 8 ;
            %else %if &firstvarcum2 = 1 %then %let type1 = 9 ;
            %else %if &firstvarcum3 = 1 %then %let type1 = 10;
            %else %let type1 = -1 ;

            %if &secondvarbin = 1       %then  %let type2 = 0 ;
            %else %if &secondvarcat = 1 %then  %let type2 = 1 ;
            %else %if &secondvarspl = 1 %then  %let type2 = 2 ;
            %else %if &secondvarqdc = 1 %then  %let type2 = 3 ;
            %else %if &secondvarcub = 1 %then  %let type2 = 4 ;
            %else %if &secondvarzqdc = 1 %then %let type2 = 5 ;
            %else %if &secondvartss1 = 1 %then %let type2 = 6 ;
            %else %if &secondvarcum0 = 1 %then %let type2 = 7 ;
            %else %if &secondvarcum1 = 1 %then %let type2 = 8 ;
            %else %if &secondvarcum2 = 1 %then %let type2 = 9 ;
            %else %if &secondvarcum3 = 1 %then %let type2 = 10;
            %else                              %let type2 = -1 ;
            /* when one factor is a baseline variable, both of the corresponding type variables will be 0. */


            %if &fortype = createvar %then %do;
               %if         &firsttype = 1 &  &secondtype = 0 %then %makecatint(&vrbl,&iota,&first,&secondind, &type1,&type2,&mytype,&baselineind);               
               %else %if   &firsttype = 1 &  &secondtype =  1%then %makecatint(&vrbl,&iota,&first,&second,    &type1,&type2,&mytype,&baselineind);               
               %else %if   &firsttype = 0 & &secondtype =  1 %then %makecatint(&vrbl,&iota,&second ,&firstind,&type2,&type1,&mytype,&baselineind) ;               
               %else %do;                  
                    label &&vrbl._I&iota._1_1="&firstvar*&secondvar";
                    &&vrbl._I&iota._1_1=&firstvar*&secondvar;             
               %end;
            %end;
            %else %if &fortype = createlist %then %do;
                              
               %if         &firsttype = 1 &  &secondtype = 0 %then  %let &vrbl._I&iota = %listcatint(&vrbl,&iota,&first,&second,&type1,&type2 ) ;
               %else %if   &firsttype = 1 &  &secondtype = 1 %then  %let &vrbl._I&iota = %listcatint(&vrbl,&iota,&first,&second,&type1,&type2 ) ;
               %else %if   &firsttype = 0 &  &secondtype = 1 %then  %let &vrbl._I&iota = %listcatint(&vrbl,&iota,&second,&first,&type2,&type1 ) ;      
               %else %if   &firsttype = 0 &  &secondtype = 0 %then  %let &vrbl._I&iota = &vrbl._I&iota._1_1;
                 
               %if &printlogstats = 1 %then %put  &vrbl._I&iota = &&&vrbl._I&iota;    
            %end;
        %end; /*iota*/

    %end; /*if there are interactions*/
    %else %do;
       %let &vrbl.ninterx = 0;
    %end;

 
%mend interactions;

/*THIS MACRO CREATES A LIST OF INTERACTION TERMS TO BE USED IN THE MODEL FOR OUTCOME OR CENSORING.*/

%macro interxarrays(type,vrbl,indexer);
   %local ii ;
   %if &&&vrbl.ninterx^=0  %then %do;

      %if &type=main %then %do;
         %do ii=1 %to %eval(&&&vrbl.ninterx);      
            &&&vrbl._I&ii
         %end;
      %end;   
      %else %if &type=sim %then %do;
 
         %do ii=1 %to %eval(&&&vrbl.ninterx);
            s&vrbl._I&ii._&indexer 
         %end;
      %end;

   %end;   
   %else %do; 
       %str() 
   %end;
%mend interxarrays;   

  
/***********/
 


%macro listcatint(vrbl,iota,first,second,type1,type2);

%local level1 level2 nlevels1 nlevels2;
/* these counts when at least one of the factors is a spline variable. There are lev-2 spline variables + original
   variable */
/* type1 and type2 = type of variable for counting number of levels in loop 
   1 = cat (lev - 1) , 2 = spl (lev - 2 + 1)  , 3 = qdc (lev = 2) , 4 = cub (lev = 3), 5 = zqdc (lev = 3) 6 = tsswitch1 (lev = 1) 
   7 = cumavg (lev = 1) 8 = lag1cumavg (lev = 2) 9 = lag2cumavg (lev = 3) 10= rcumavg (lev = 1) */

/* when usespline = 1, there can be additional spline variables for types 6,8 and 9.
   for type = 6 there is one spline for each term
       type = 7 can have spline cov.lev - 2 spline variables
       type = 8, 9 there are cov.lev - 2 spline variables for each variable */


%let nlevels1 = 1 ; /* default value for bin type and baseline variables */

%if &type1 = 1 or &type1 = 2 %then %let nlevels1 = %eval(&&cov&first.lev - 2 + 1); 
%else %if &type1 = 3 %then %let nlevels1 = 2 ;
%else %if &type1 = 4 %then %let nlevels1 = 3 ;
%else %if &type1 = 5 %then %let nlevels1 = 3 ;
%else %if &type1 = 6 %then %let nlevels1 = %eval(2 +  2* (&&cov&first.lev - 2)) ;
%else %if &type1 = 7 %then %let nlevels1=  %eval(1 +    (&&cov&first.lev - 2));
%else %if &type1 = 8 %then %let nlevels1 = %eval(2 + 2* (&&cov&first.lev - 2)) ;
%else %if &type1 = 9 %then %let nlevels1 = %eval(3 + 3* (&&cov&first.lev - 2)) ;
%else %if &type1 = 10 %then %let nlevels1 = 1;


%let nlevels2 = 1 ;  /* default value for bin type and baseline variables */

%if &type2 = 1 or &type2 = 2 %then %let nlevels2 = %eval(&&cov&second.lev-1); 
%else %if &type2 = 3 %then  %let nlevels2 = 2 ;
%else %if &type2 = 4 %then  %let nlevels2 = 3 ;
%else %if &type2 = 5 %then  %let nlevels2 = 3 ;
%else %if &type2 = 6 %then  %let nlevels2 = %eval(2 +  2* (&&cov&first.lev - 2))  ;
%else %if &type2 = 7 %then  %let nlevels2 = %eval(1 +     (&&cov&second.lev - 2)) ;
%else %if &type2 = 8 %then  %let nlevels2 = %eval(2 +  2 *(&&cov&second.lev - 2)) ;
%else %if &type2 = 9 %then  %let nlevels2 = %eval(3 +  3 *(&&cov&second.lev - 2)) ;
%else %if &type2 = 10 %then %let nlevels2 = 1;

%do level1=1 %to &nlevels1 ;
     %do level2 = 1 %to &nlevels2 ; 
         &&vrbl._I&iota._&level1._&level2
     %end; 
%end;
%mend listcatint;



%macro makecatint(vrbl,iota,first,second,type1,type2,vartype,baseline);
  /* when first and second are spline variables, the same number of factors is used, except the variable
     names need to be changed. */

  /* vartype will indicate an outcome-type variable or covariate-type variable */

   %local level1 level2 nlevels1 nlevels2 var1 var2 splinelevel1 splinelevel2 mylag1 mylag2 lag1 lag2 lag3 spl1 spl2 varlevel1 varlevel2 ;

   /* type is indicator of ptype 1 = cat, 2 = spl , 3 = qdc, 4 = cub  , 5 = zqdc
      nlevels changes depding of type, 1 and 2 have lev - 1 , 3 has 2, 4 has 3  , 5 has 3 */

  %let nlevels1 = 1 ; /* for -bin,  rcumavg, and baseline variables, type1 = -1 (?? baseline conbin )  */


  %if &type1 = 1 or &type1 = 2 %then %let nlevels1 = %eval(&&cov&first.lev-1); 
  %else %if &type1 = 3 %then %let nlevels1 = 2 ;
  %else %if &type1 = 4 %then %let nlevels1 = 3 ;
  %else %if &type1 = 5 %then %let nlevels1 = 3 ;
  %else %if &type1 = 6 %then %let nlevels1 = %eval(2+  2 * %sysfunc(max(0,(&&cov&first.lev - 2 )))) ;
  %else %if &type1 = 7 %then %let nlevels1 = %eval(1+  1 * %sysfunc(max(0,(&&cov&first.lev - 2 )))) ;
  %else %if &type1 = 8 %then %let nlevels1 = %eval(2+  2 * %sysfunc(max(0,(&&cov&first.lev - 2 )))) ;
  %else %if &type1 = 9 %then %let nlevels1 = %eval(3+  3 * %sysfunc(max(0,(&&cov&first.lev - 2 )))) ;
  


  %let nlevels2 = 1 ;

  %if &type2 = 1 or &type2 = 2 %then %let nlevels2 = %eval(&&cov&second.lev-1); 
  %else %if &type2 = 3 %then %let nlevels2 = 2 ;
  %else %if &type2 = 4 %then %let nlevels2 = 3 ;
  %else %if &type2 = 5 %then %let nlevels2 = 3 ;
  %else %if &type2 = 6 %then %let nlevels2 = %eval(2+  2 * %sysfunc(max (0,(&&cov&second.lev - 2)))) ;
  %else %if &type1 = 7 %then %let nlevels1 = %eval(1+ %sysfunc(max(0,(&&cov&second.lev - 2 )))) ;
  %else %if &type2 = 8 %then %let nlevels2 = %eval(2+  2 * %sysfunc(max (0,(&&cov&second.lev - 2)))) ;
  %else %if &type2 = 9 %then %let nlevels2 = %eval(3+  3 * %sysfunc(max(0, (&&cov&second.lev - 2)))) ;
 

 %if &vartype = 1 %then %do ;
      /* no lag values for outcome-type model interactions */
      %let mylag1 = ;
      %let mylag2 = ;
 %end;
 %else %if &vartype = 2 %then %do;
    /* only use lagged values for covariate models */
    %let mylag1 = _l1 ;
    %if %upcase(%substr(&&cov&first.ptype , 1,3)) = CON %then %let mylag1 = ;
    %let mylag2 = _l1 ;
    %if %upcase(%substr(&&cov&second.ptype , 1,3)) = CON %then %let mylag2 = ;     
 %end ;

%let var1 = &&cov&first..&mylag1 ;
%let var2 = &&cov&second..&mylag2 ;


/* when there is a single baseline variable, this is contained in the second variable. there will be no lag */
%if &baseline = 1 %then %let var2 = %scan(&fixedcov,&second); 
  
/* when there is a baseline variable it will be in the second variable and &second will indicate the position in fixedcov  */

   %let splinelevel1 = 0 ;
   %do level1=1 %to &nlevels1 ;
       %if &type1 > 1 %then %let var1 =  ; 
      

      %if &type1 = 1 %then %do;
                %let var1 =  &&cov&first..&mylag1._&level1 ;
      %end;
      %else %if &type1 = 2 %then %do;
                %if &level1 = 1 %then %let var1 = &&cov&first..&mylag1 ;
                %else %let var1 = &&cov&first..&mylag1._spl%eval(&level1-1) ;
      %end;
      %else %if &type1 = 3 %then %do;
               %if &level1 = 1 %then %let var1 = &&cov&first.&mylag1;
               %else %if &level1 = 2 %then %let var1 = &&cov&first..&mylag1.s ;
      %end;
      %else %if &type1 = 4 %then %do;
                 %if &level1 = 1 %then %let var1 = &&cov&first..&mylag1;
                 %else %if &level1 = 2 %then %let var1 = &&cov&first..&mylag1.s ;
                 %else %if &level1 = 3 %then %let var1 = &&cov&first..&mylag1.c ;
      %end;
      %else %if &type1 = 5 %then %do;
                 %if &level1 = 1 %then %let var1 = z&&cov&first.&mylag1;
                 %else %if &level1 = 2 %then %let var1 = &&cov&first.&mylag1 ;
                 %else %if &level1 = 3 %then %let var1 = &&cov&first..&mylag1.s ;
      %end;
      
      %else %if &type1 = 6 %then %do;
             %if       &level1 = 1 %then %let var1 = &&cov&first..&mylag1;
             %else %if &level1 = 2 %then %let var1 = ts&&cov&first..&mylag1._inter ;           
             %else %if &level1 = 3 %then %let var1 = ts&&cov&first..&mylag1._inter_spl1 ;             
             %else %if &level1 = 4 %then %let var1 = ts&&cov&first..&mylag1._inter_spl2 ;
      %end; 
      %else %if &type1 = 7 %then %let var1 = &&cov&first.._cumavg&mylag1 ;

      %else %if &type1 = 8 %then %do;
              %if &vartype = 1 %then %do;
                %let lag1 = ;
                %let lag2 = _l1 ;
                 
           %end;
           %else %if &vartype = 2 %then %do ;
               %let lag1 = _l1 ;
               %let lag2 = _l2 ;              
          %end; 
      
          %let varlevel1 = %sysfunc(mod(&level1,2)) ;
          %if &splinelevel1 = 0 %then %let spl1 = ;
          %else %let spl1 = _spl&splinelevel1 ;
          
          
  
               %if &varlevel1 = 1 %then %let var1= &&cov&first..&lag1.&spl1    ;
         %else %if &varlevel1 = 0 %then %let var1 =  &&cov&first.._cumavg&lag2.&spl1     ;

         %if &varlevel1 = 0 %then %let splinelevel1 = %eval(&splinelevel1 + 1) ;
      %end;
      %else %if &type1 = 9 %then %do;
            %if &vartype = 1 %then %do;
                %let lag1 = ;
                %let lag2 = _l1 ;
                %let lag3 = _l2 ;
           %end;
           %else %if &vartype = 2 %then %do ;
               %let lag1 = _l1 ;
               %let lag2 = _l2 ;
               %let lag3 = _l3 ;
          %end; 


          %let varlevel1 = %sysfunc(mod(&level1,3)) ;
          %if &splinelevel1 = 0 %then %let spl1 = ;
          %else %let spl1 = _spl&splinelevel1 ;
          
          
  
               %if &varlevel1 = 1 %then %let var1= &&cov&first..&lag1.&spl1    ;
         %else %if &varlevel1 = 2 %then %let  var1= &&cov&first..&lag2.&spl1   ;   
         %else %if &varlevel1 = 0 %then %let var1 =  &&cov&first.._cumavg&lag3.&spl1     ;

         %if &varlevel1 = 0 %then %let splinelevel1 = %eval(&splinelevel1 + 1) ;

 
      %end;

      %else %if &type1 = 10 %then %let var1 = &&cov&first.._rcumavg&mylag1 ;
 %let splinelevel2 = 0 ;
      %do level2 = 1 %to &nlevels2;
           %if &type2 > 1 %then %let var2 =  ;
           
           %if &type2 = 1 %then %let var2 = &&cov&second..&mylag2._&level2 ;
           %else %if &type2 = 2 %then %do;
                 %if &level2 = 1 %then %let var2 = &&cov&second.&mylag2 ;
                %else %let var2 = &&cov&second..&mylag2._spl%eval(&level2-1) ;
          %end;
           %else %if &type2 = 3 %then %do;
               %if &level2 = 1 %then %let var2 = &&cov&second.&mylag2;
               %else %if &level2 = 2 %then %let var2 = &&cov&second..&mylag2.s ;
           %end;
           %else %if &type2 = 4 %then %do;
                 %if &level2 = 1 %then %let var2 = &&cov&second.&mylag2;
                 %else %if &level2 = 2 %then %let var2 = &&cov&second..&mylag2.s ;
                 %else %if &level2 = 3 %then %let var2 = &&cov&second..&mylag2.c ;
           %end;
            %else %if &type2 = 5 %then %do;
                 %if &level2 = 1 %then %let var2 = z&&cov&second.&mylag2;
                 %else %if &level2 = 2 %then %let var2 = &&cov&second.&mylag2 ;
                 %else %if &level2 = 3 %then %let var2 = &&cov&second..&mylag2.s ;
           %end;
           %else %if &type2 = 6 %then %do;
             %if &level2 = 1        %then %let var2 = &&cov&second..&mylag1 ;
             %else %if &level2 = 2  %then %let var2 = ts&&cov&second..&mylag1._inter ;              
             %else %if &level2 = 3  %then %let var2 = ts&&cov&second..&mylag1._inter_spl1 ;              
             %else %if &level2 = 4  %then %let var2 = ts&&cov&second..&mylag1._inter_spl2 ;
           %end; 
           %else %if &type2 = 7 %then %let var2 = &&cov&second.._cumavg&mylag1 ;
           %else %if &type2 = 8 %then %do;
                 %if &vartype = 1 %then %do;
                      %let lag1 = ;
                      %let lag2 = _l1 ;                 
                %end;
                %else %if &vartype = 2 %then %do ;
                    %let lag1 = _l1 ;
                    %let lag2 = _l2 ;              
               %end; 

         
              %let varlevel2 = %sysfunc(mod(&level2,2)) ;
              %if &splinelevel2 = 0 %then %let spl2 = ;
              %else %let spl2 = _spl&splinelevel2 ;
          
          
  
               %if &varlevel2 = 1 %then %let var2= &&cov&second..&lag1.&spl2    ;
               %else %if &varlevel2 = 0 %then %let var2 =  &&cov&second.._cumavg&lag2.&spl2     ;

               %if &varlevel2 = 0 %then %let splinelevel2 = %eval(&splinelevel2 + 1) ;
         
          %end;
          %else %if &type2 = 9 %then %do;
            %if &vartype = 1 %then %do;
                %let lag1 = ;
                %let lag2 = _l1 ;
                %let lag3 = _l2 ;
           %end;
           %else %if &vartype = 2 %then %do ;
               %let lag1 = _l1 ;
               %let lag2 = _l2 ;
               %let lag3 = _l3 ;
          %end; 
  
          %let varlevel2 = %sysfunc(mod(&level2,3)) ;
          %if &splinelevel2 = 0 %then %let spl2 = ;
          %else %let spl2 = _spl&splinelevel2 ;
                     
               %if &varlevel2 = 1 %then %let var2= &&cov&second..&lag1.&spl2    ;
         %else %if &varlevel2 = 2 %then %let  var2= &&cov&second..&lag2.&spl2      ;
         %else %if &varlevel2 = 0 %then %let var2 =  &&cov&second.._cumavg&lag3.&spl2     ;

         %if &varlevel2 = 0 %then %let splinelevel2 = %eval(&splinelevel2 + 1) ;          
      %end;    
      %else %if &type2 = 10 %then %let var2 = &&cov&second.._rcumavg&mylag2 ;   

          label &&vrbl._I&iota._&level1._&level2 = "&var1 * &var2" ;
          &&vrbl._I&iota._&level1._&level2 = &var1 * &var2 ;
         
       %end; /* level2 */
   %end;   /* level1 */
   
   
%mend makecatint;
/*****/


%macro remove_created_global;       

            /* remove the global macro variables that were created in the macros used to create
               the interaction terms used in the models */

            %local  ii jj  rw nvars n_removed test word _var_removed removed ;
            %let ii = 1;
            %let _var_removed= ;
            %let n_removed = 0;
            %let nvars = %numargs(&created_global) ;
            %do ii = 1 %to %eval(&nvars) ;
           

                %let test = %scan(&created_global,%eval(&ii),%str( ));
                %let removed = 0;
                %do jj = 1 %to %eval(&n_removed);
                    %let rw = %scan(&_var_removed,%eval(&jj),%str( ));
                    %if &test = &rw %then %let removed = 1 ;
                %end;
                %if &removed = 0 %then %do;
                       %symdel &test ;
                       %let _var_removed = &_var_removed &test ;
                       %let n_removed = %eval(&n_removed + 1);
                %end;              
            %end; 
  
%mend remove_created_global;
/**********************************************************************************************************************************/
/**************************************************************/
/**************************************************************/
/**** THIS MACRO IS CALLED BY DATAPREP, INTERV, AND SIMVAR ****/
/**** IT CREATES THE VARIABLES NEEDED FOR THE PREDICTORS, *****/
/**** DEPENDING ON THEIR TYPE. ********************************/
/**************************************************************/
/**************************************************************/

%macro genpred(type,lagtype=0,useeof = 0);

        /* IMPORTATNT : THE VARIABLE I IS NOT A LOCAL VARIABLE, BUT IS USED IN THE CALLING PROGRAM FOR SETTING THE INDEX OF THE VARIABLE. COULD 
           CHANGE THIS SO IT IS AN ARGUMENT OF THE MACRO */
 
        %local current lagged  lev knot ;
        %let current = 1 ;
        %let lagged = 1 ;
        %if &lagtype = 1 %then %let lagged  = 0 ;
        %else %if &lagtype = 2 %then %let current = 0 ;      
       

/* ALGEBRAIC MANIPULATIONS, ISH */
      
%if &useeof = 0 %then %do; 


        %if &&cov&i.ptype = conqdc or &&cov&i.ptype = conzqdc    %then %do;
            &&cov&i..s = &&cov&i*&&cov&i;             /* SQUARE */
        %end;

        %if &&cov&i.ptype = concub  %then %do;
            &&cov&i..s = &&cov&i*&&cov&i;             /* SQUARE */
            &&cov&i..c = &&cov&i*&&cov&i*&&cov&i;           /* CUBE */
        %end;
            
        %if &&cov&i.ptype = conspl AND %numargs(&&cov&i.knots) > 1    %then %do;
            %rcspline(&&cov&i,&&cov&i.knots);           /* SPLINE, SEE BELOW! */
        %end;
            
        %if &&cov&i.ptype = concat   %then %do;    
            %if &&cov&i.otype ne 5 AND %numargs(&&cov&i.knots) > 1  %then %do;             /* MAKE CATEGORIES, SEE BELOW! */
                %makecat(&&cov&i, &&cov&i.knots, &&cov&i.lev);
            %end;
            %if &&cov&i.otype = 5 %then %do;
                %do lev = 1 %to %eval(&&cov&i.lev - 1);
                     &&cov&i.._&lev = (&&cov&i = &lev);
                %end;
           %end;  
        %end;

/* MAKE LAGS */
    /* LAG 1 */


        %if (&&cov&i.ptype = lag1qdc  or &&cov&i.ptype = lag1zqdc )  %then %do;
            %if &current = 1 %then &&cov&i..s = &&cov&i*&&cov&i;;
            %if &lagged = 1  %then &&cov&i..s_l1 = &&cov&i.._l1*&&cov&i.._l1;;
        %end;

        %if (&&cov&i.ptype = lag1cub )   %then %do;
            %if &current = 1 %then  &&cov&i..s = &&cov&i*&&cov&i;;
            %if &lagged = 1 %then &&cov&i..s_l1 = &&cov&i.._l1*&&cov&i.._l1;;
            %if &current = 1 %then &&cov&i..c = &&cov&i*&&cov&i*&&cov&i;;
            %if &lagged = 1 %then &&cov&i..c_l1 = &&cov&i.._l1*&&cov&i.._l1*&&cov&i.._l1;;
       %end;
            
        %if &&cov&i.ptype = lag1spl AND %numargs(&&cov&i.knots) > 1   %then %do;
            %if &current = 1 %then %rcspline(&&cov&i,&&cov&i.knots);
            %if &lagged = 1 %then %rcspline(&&cov&i.._l1,&&cov&i.knots);
        %end;
            
        %if &&cov&i.ptype = lag1cat AND %numargs(&&cov&i.knots) > 1    %then %do;
           %if &current = 1 %then  %makecat(&&cov&i, &&cov&i.knots, &&cov&i.lev);
           %if &lagged = 1  %then  %makecat(&&cov&i.._l1, &&cov&i.knots, &&cov&i.lev);
        %end;
    
    /* LAG 2 */ 
        %if &&cov&i.ptype = lag2qdc or &&cov&i.ptype = lag2zqdc   %then %do;
            %if &current = 1 %then &&cov&i..s = &&cov&i*&&cov&i;;
            %if &lagged = 1 %then  &&cov&i..s_l1 = &&cov&i.._l1*&&cov&i.._l1;;
            %if &lagged = 1 %then &&cov&i..s_l2 = &&cov&i.._l2*&&cov&i.._l2;;
        %end;

        %if &&cov&i.ptype = lag2cub   %then %do;
             %if &current = 1 %then &&cov&i..s = &&cov&i*&&cov&i;;
             %if &lagged = 1 %then &&cov&i..s_l1 = &&cov&i.._l1*&&cov&i.._l1;;
             %if &lagged = 1 %then &&cov&i..s_l2 = &&cov&i.._l2*&&cov&i.._l2;;
             %if &current = 1 %then &&cov&i..c = &&cov&i*&&cov&i*&&cov&i;;
             %if &lagged = 1 %then &&cov&i..c_l1 = &&cov&i.._l1*&&cov&i.._l1*&&cov&i.._l1;;
             %if &lagged = 1 %then &&cov&i..c_l2 = &&cov&i.._l2*&&cov&i.._l2*&&cov&i.._l2;;
        %end;
            
        %if &&cov&i.ptype = lag2spl AND %numargs(&&cov&i.knots) > 1  %then %do;
             %if &current = 1 %then %rcspline(&&cov&i,&&cov&i.knots);
             %if &lagged = 1 %then %rcspline(&&cov&i.._l1,&&cov&i.knots);
             %if &lagged = 1 %then %rcspline(&&cov&i.._l2,&&cov&i.knots);
        %end;
            
        %if &&cov&i.ptype = lag2cat AND %numargs(&&cov&i.knots) > 1  %then %do;
            %if &current = 1 %then %makecat(&&cov&i, &&cov&i.knots, &&cov&i.lev);
            %if &lagged = 1 %then  %makecat(&&cov&i.._l1, &&cov&i.knots, &&cov&i.lev);
            %if &lagged = 1 %then %makecat(&&cov&i.._l2, &&cov&i.knots, &&cov&i.lev);
       %end;

    /* LAG 3 */
        %if &&cov&i.ptype = lag3qdc or &&cov&i.ptype = lag3zqdc  
              %then %do;
            %if &current = 1 %then &&cov&i..s = &&cov&i*&&cov&i;;
            %if &lagged = 1 %then  &&cov&i..s_l1 = &&cov&i.._l1*&&cov&i.._l1;;
            %if &lagged = 1 %then  &&cov&i..s_l2 = &&cov&i.._l2*&&cov&i.._l2;;
             %if &lagged = 1 %then &&cov&i..s_l3 = &&cov&i.._l3*&&cov&i.._l3;;
        %end;

        %if &&cov&i.ptype = lag3cub    %then %do;
            %if &current = 1 %then &&cov&i..s = &&cov&i*&&cov&i;;
            %if &lagged = 1 %then  &&cov&i..s_l1 = &&cov&i.._l1*&&cov&i.._l1;;
            %if &lagged = 1 %then  &&cov&i..s_l2 = &&cov&i.._l2*&&cov&i.._l2;;
            %if &lagged = 1 %then  &&cov&i..s_l3 = &&cov&i.._l3*&&cov&i.._l3;;
            %if &current = 1 %then &&cov&i..c = &&cov&i*&&cov&i*&&cov&i;;
            %if &lagged = 1 %then  &&cov&i..c_l1 = &&cov&i.._l1*&&cov&i.._l1*&&cov&i.._l1;;
             %if &lagged = 1 %then &&cov&i..c_l2 = &&cov&i.._l2*&&cov&i.._l2*&&cov&i.._l2;;
            %if &lagged = 1 %then  &&cov&i..c_l3 = &&cov&i.._l3*&&cov&i.._l3*&&cov&i.._l3;;
       %end;

        
        %if &&cov&i.ptype = lag3spl AND %numargs(&&cov&i.knots) > 1   %then %do;
            %if &current = 1 %then %rcspline(&&cov&i,&&cov&i.knots);
            %if &lagged = 1 %then %rcspline(&&cov&i.._l1,&&cov&i.knots);
            %if &lagged = 1 %then %rcspline(&&cov&i.._l2,&&cov&i.knots);
            %if &lagged = 1 %then  %rcspline(&&cov&i.._l3,&&cov&i.knots);
        %end;
            
        %if &&cov&i.ptype = lag3cat AND %numargs(&&cov&i.knots) > 1   %then %do;
            %if &current = 1 %then %makecat(&&cov&i, &&cov&i.knots, &&cov&i.lev);
            %if &lagged = 1 %then  %makecat(&&cov&i.._l1, &&cov&i.knots, &&cov&i.lev);
            %if &lagged = 1 %then  %makecat(&&cov&i.._l2, &&cov&i.knots, &&cov&i.lev);
            %if &lagged = 1 %then  %makecat(&&cov&i.._l3, &&cov&i.knots, &&cov&i.lev);
       %end;

/* IF TIMES WERE SKIPPED FOR THIS VARIABLE, WE NEED THE INTERACTION WITH TIME. */
 
  
      %if &&cov&i.ptype = skpbin    %then %do;
           %if &current = 1 %then  %maketi(&&cov&i, &time, &time._l1, &&cov&i.skip, &interval);
           %if &lagged = 1 %then  %maketi(&&cov&i.._l1, &time._l1, &time._l2, &&cov&i.skip, &interval);            
      %end;

        %if &&cov&i.ptype = skpqdc or &&cov&i.ptype = skpzqdc   
               %then %do;
            %if &current = 1 %then %maketi(&&cov&i, &time, &time._l1, &&cov&i.skip, &interval);
            %if &lagged = 1 %then  %maketi(&&cov&i.._l1, &time._l1, &time._l2, &&cov&i.skip, &interval);
            
            %if &current = 1 %then &&cov&i..s = &&cov&i*&&cov&i;;
            %if &lagged = 1 %then  &&cov&i..s_l1 = &&cov&i.._l1*&&cov&i.._l1;;
            
            %if &current = 1 %then %maketi(&&cov&i..s, &time, &time._l1, &&cov&i.skip, &interval);
            %if &lagged = 1 %then  %maketi(&&cov&i..s_l1, &time._l1, &time._l2, &&cov&i.skip, &interval);            
       %end;


        %if &&cov&i.ptype = skpcub   %then %do;
            %if &current = 1 %then %maketi(&&cov&i, &time, &time._l1, &&cov&i.skip, &interval);
            %if &lagged = 1 %then  %maketi(&&cov&i.._l1, &time._l1, &time._l2, &&cov&i.skip, &interval);
            
            %if &current = 1 %then  &&cov&i..s = &&cov&i*&&cov&i;;
            %if &lagged = 1 %then  &&cov&i..s_l1 = &&cov&i.._l1*&&cov&i.._l1;;
            %if &current = 1 %then  &&cov&i..c = &&cov&i*&&cov&i*&&cov&i;;
            %if &lagged = 1 %then  &&cov&i..c_l1 = &&cov&i.._l1*&&cov&i.._l1*&&cov&i.._l1;;

            %if &current = 1 %then %maketi(&&cov&i..s, &time, &time._l1, &&cov&i.skip, &interval);
            %if &lagged = 1 %then  %maketi(&&cov&i..s_l1, &time._l1, &time._l2, &&cov&i.skip, &interval);
            %if &current = 1 %then  %maketi(&&cov&i..c, &time, &time._l1, &&cov&i.skip, &interval);
            %if &lagged = 1 %then   %maketi(&&cov&i..c_l1, &time._l1, &time._l2, &&cov&i.skip, &interval);
       %end;

                        
        %if &&cov&i.ptype = skpspl AND %numargs(&&cov&i.knots) > 1   %then %do;
            %if &current = 1 %then %maketi(&&cov&i,&time,&time._l1,&&cov&i.skip, &interval);
            %if &lagged = 1 %then  %maketi(&&cov&i.._l1,&time._l1,&time._l2,&&cov&i.skip, &interval);
            
            %if &current = 1 %then  %rcspline(&&cov&i,&&cov&i.knots);
            %if &lagged = 1 %then  %rcspline(&&cov&i.._l1,&&cov&i.knots);
            
            %do knot = 1 %to %eval(&&cov&i.lev - 2);
               %if &current = 1 %then  %maketi(&&cov&i.._spl&knot,&time,&time._l1,&&cov&i.skip, &interval);
               %if &lagged = 1 %then    %maketi(&&cov&i.._l1_spl&knot,&time._l1,&time._l2,&&cov&i.skip, &interval);
            %end;
        %end;
            
        %if &&cov&i.ptype = skpcat AND %numargs(&&cov&i.knots) > 1   %then %do;

            %if &current = 1 %then  %makecat(&&cov&i, &&cov&i.knots, &&cov&i.lev);
            %if &lagged = 1 %then  %makecat(&&cov&i.._l1, &&cov&i.knots, &&cov&i.lev);
            
            %do lev = 1 %to %eval(&&cov&i.lev - 1);
               %if &current = 1 %then  %maketi(&&cov&i.._&lev,&time,&time._l1,&&cov&i.skip, &interval);
               %if &lagged = 1 %then   %maketi(&&cov&i.._l1_&lev,&time._l1,&time._l2,&&cov&i.skip, &interval);
            %end;
        %end; /* SKPCAT */


       
         %if &&cov&i.ptype = tsswitch1  %then %do; /*change JGY*/  
            %if &type = main %then %do;
               
               /*ts&&cov&i.._l1_inter = lag(ts&&cov&i.._inter) ; */
               retain  ts&&cov&i.._inter ;;
               if &time = 0 then do;                             
                   ts&&cov&i.._inter = 0 ;
                   ts&&cov&i.._l1_inter = 0 ;
               end;
               ts&&cov&i.._l1_inter = ts&&cov&i.._inter ;
               ts&&cov&i.._inter = ts&&cov&i.._l1_inter + &&cov&i ;

               *sumtemp&i = sumtemp&i + &&cov&i; 
               *sumtemp_l1&i = sumtemp_l1&i + &&cov&i.._l1;
               *ts&&cov&i.._inter = &&cov&i*sumtemp&i;
               *ts&&cov&i.._l1_inter = &&cov&i.._l1*sumtemp_l1&i;                   
            %end;
            %else %if &type = sim %then %do;                            
               %if &lagged = 1 %then %do;
                   if &time = 0 then do ;
                       ts&&cov&i.._l1_inter = 0 ;
                   end;
               %end;
               %if &current = 1 %then  ts&&cov&i.._inter = ts&&cov&i.._l1_inter + &&cov&i ;;                                        
               %if /* &usespline = 1 */ %numargs(&&cov&i.knots) >=  3 %then %do;
                     %if &current = 1 %then  %rcspline(ts&&cov&i.._inter ,&&cov&i.knots);
                     %if &lagged  = 1 %then  %rcspline(ts&&cov&i.._l1_inter ,&&cov&i.knots);
               %end;         
            %end;
        %end;

        %if ( &&cov&i.ptype=cumavg |  &&cov&i.ptype=lag1cumavg | &&cov&i.ptype=lag2cumavg | &&cov&i.ptype = cumavgcat |
               &&cov&i.ptype=lag1cumavgcat | &&cov&i.ptype = lag2cumavgcat ) 
           %then %do;
            %if &type = main %then %do;
                retain sum_&&cov&i  %if %bquote(&&cov&i.cumint) ^= %then &&cov&i.._timestart_ &&cov&i..sumtimestart ;;

                if first.newid then do ;               
                    sum_&&cov&i = 0 ;  
                    %if %bquote(&&cov&i.cumint) ^= %then &&cov&i.._timestart_ = -1 ; ;
                end;

                %if %bquote(&&cov&i.cumint) ^= %then %do;
                    if &&cov&i.._timestart_ = -1 and &&cov&i.cumint = 1 then &&cov&i.._timestart_ = &time ;
                %end;


                _counter = &time + 1 ;
                %if %bquote(&&cov&i.cumint) = %then %do;
                    sum_&&cov&i = sum_&&cov&i + &&cov&i ;           
                    &&cov&i.._cumavg = sum_&&cov&i / _counter ; 
                    &&cov&i.._cumavg_l1 =  lag(&&cov&i.._cumavg) ; 
                    &&cov&i.._cumavg_l2 =  lag2(&&cov&i.._cumavg) ; 
                    &&cov&i.._cumavg_l3 =  lag3(&&cov&i.._cumavg)  ; 

                    if &time < 1 then do ;
                        &&cov&i.._cumavg_l1 =  0 ; 
                        &&cov&i.._cumavg_l2 =  0 ; 
                        &&cov&i.._cumavg_l3 =  0  ; 
                    end ;
                    else if &time < 2 then do ;

                        &&cov&i.._cumavg_l2 =  0 ; 
                        &&cov&i.._cumavg_l3 =  0  ; 
                    end ;
                    else if &time < 3 then do ;
                        &&cov&i.._cumavg_l3 =  0  ;
                    end;
                    /**
                    %if &&cov&i.ptype = cumavgcat %then %do ;
                        %if &current = 1 %then  %makecat(&&cov&i.._cumavg, &&cov&i.knots, &&cov&i.lev);
                        %if &lagged = 1  %then  %makecat(&&cov&i.._cumavg_l1, &&cov&i.knots, &&cov&i.lev);
                        
                    %end;
                    ***/
                %end;

                %if %bquote(&&cov&i.cumint) ^= %then %do;
                    &&cov&i.._newcum = 0 ;
                    if &&cov&i.._timestart_ ge 0 then do ;
                        if &&cov&i.._timestart_ = &time then do ;
                            &&cov&i..sumtimestart = sum_&&cov&i - &&cov&i ;   /* sum (m = 0 to m=timestart - 1 L(m) */             
                        end;   
                        &&cov&i.._newcum = (sum_&&cov&i - &&cov&i..sumtimestart) / (_counter - &&cov&i.._timestart_ ) ; 
                    end;
                    &&cov&i.._cumavg_l3 = lag3(&&cov&i.._newcum);
                    &&cov&i.._cumavg_l2 = lag2(&&cov&i.._newcum);
                    &&cov&i.._cumavg_l1 = lag1(&&cov&i.._newcum);
                    if &time < 1 then do ;
                        &&cov&i.._cumavg_l3 = 0;
                        &&cov&i.._cumavg_l2 =0;
                        &&cov&i.._cumavg_l1 = 0;
                    end;
                    if &time < 2 then do ;
                        &&cov&i.._cumavg_l3 = 0 ;
                        &&cov&i.._cumavg_l2 =  0;
                    end;
                    else if &time < 3 then do ;                
                        &&cov&i.._cumavg_l3 = 0;
                    end;                      
                    drop &&cov&i.._newcum ;
                %end;

                drop  _counter sum_&&cov&i %if %bquote(&&cov&i.cumint) ^= %then &&cov&i..sumtimestart &&cov&i.._timestart_ ; ;   
            %end;
            %else %if &type = sim %then %do ;
                %if %bquote(&&cov&i.cumint)= %then %do;
                    %if &lagged = 1 %then %do;
                        if &time = 0  then do;
                            /* set all negative-time averages to be 0, k = 3 is baseline  */
                            &&cov&i.._cumavg_l3 = 0 ;
                            &&cov&i.._cumavg_l2 = 0 ;
                            &&cov&i.._cumavg_l1 = 0 ;                                          
                        end;
                    %end;
                     /*due to a possible change of the current value of variable due to intervention, we need to 
                       recalculate the average. will use the mean function with the assumption that future times are missing */
                       %if &current = 1 %then &&cov&i.._cumavg = (&time * &&cov&i.._cumavg_l1 + &&cov&i )/(&time + 1)  ;;              
                   
                %end;
                %else %do ;
                    * unless the variable contained in covX.cumint is created before covX then current value of the cumavg term ;
                    *  will be wrong. This will be corrected in the call to genpred preceeding the creation of the outcome. ;
                    %if &lagged = 1 %then %do;
                        if &time = 0  then do;
                            &&cov&i.._timestart = -1 ;
                            &&cov&i.._cumavg_l3 = 0 ;
                            &&cov&i.._cumavg_l2 = 0 ;
                            &&cov&i.._cumavg_l1 = 0 ;
                        end;
                    %end;
                    %if &current = 1 %then %do;
                        if &&cov&i.cumint = 0 then &&cov&i.._cumavg = 0 ;
                        if &&cov&i.cumint = 1 then do;
                           if &&cov&i.._timestart = -1 then do;
                               &&cov&i.._timestart = &time ;
                               &&cov&i.._normalization_factor = 1 ;                           
                               &&cov&i.._cumavg = &&cov&i ;
                           end;
                           else if &&cov&i.._timestart > -1 then do ;
                            &&cov&i.._cumavg =  ( (&&cov&i.._normalization_factor -1) * &&cov&i.._cumavg_l1 + &&cov&i)/(&&cov&i.._normalization_factor ) ;
                             
                          end;
                       end;
                      
                    %end;      

                %end;
               /* need > 2 for splines, when covXknots = 0 there are no splines used */
                 
                %if &&cov&i.ptype = cumavg AND %numargs(&&cov&i.knots) > 2 %then %do;
                    %if &current = 1 %then %rcspline(&&cov&i.._cumavg, &&cov&i.knots);
                    %if &lagged = 1  %then %rcspline(&&cov&i.._cumavg_l1, &&cov&i.knots);
                %end;
                %else %if (&&cov&i.ptype = lag1cumavg or &&cov&i.ptype = lag2cumavg ) AND %numargs(&&cov&i.knots) > 2  %then %do;
                    %if &current = 1 %then %rcspline(&&cov&i ,&&cov&i.knots);  
                    %if &lagged = 1 %then %do;
                         %rcspline(&&cov&i.._l1 ,&&cov&i.knots);                                                                                                              
                         %rcspline(&&cov&i.._cumavg_l1 ,&&cov&i._cumavg_l1_knots);
                         %rcspline(&&cov&i.._cumavg_l2 ,&&cov&i._cumavg_l1_knots);
                         %if &&cov&i.ptype= lag2cumavg %then %do ;
                             %rcspline(&&cov&i.._l2 ,&&cov&i.knots); 
                             %rcspline(&&cov&i.._cumavg_l3 ,&&cov&i._cumavg_l1_knots);
                         %end;
                     %end;
                 %end ;
                 %else %if &&cov&i.ptype = cumavgcat %then %do ;                      
                    %if &current = 1 %then  %makecat(&&cov&i.._cumavg, &&cov&i.knots, &&cov&i.lev);
                    %if &lagged = 1  %then  %makecat(&&cov&i.._cumavg_l1, &&cov&i.knots, &&cov&i.lev);                                           
                 %end;
                  %else %if &&cov&i.ptype = lag1cumavgcat OR &&cov&i.ptype = lag2cumavgcat %then %do ;                      
                    %if &current = 1 %then  %makecat(&&cov&i, &&cov&i.knots, &&cov&i.lev);
                    %if &lagged = 1  %then  %do;
                        %makecat(&&cov&i.._l1, &&cov&i.knots , &&cov&i.lev) ;
                        %makecat(&&cov&i.._cumavg_l1, &&cov&i._cumavg_l1_knots, &&cov&i._cumavg_l1_lev);                                           
                        %makecat(&&cov&i.._cumavg_l2,&&cov&i._cumavg_l1_knots, &&cov&i._cumavg_l1_lev);
                        %if &&cov&i.ptype = lag2cumavgcat %then %do;
                            %makecat(&&cov&i.._l2, &&cov&i._cumavg_l1_knots, &&cov&i._cumavg_l1_lev);
                            %makecat(&&cov&i.._l3, &&cov&i._cumavg_l1_knots, &&cov&i._cumavg_l1_lev);
                        %end;
                    %end;
                 %end;
            %end;

            
        %end;

        %if &&cov&i.ptype=rcumavg   %then %do;
            %if &type = main %then %do;
                 retain  &&cov&i.._rcumavg ;
                 if first.newid then do ;                                  
                       &&cov&i.._rcumavg = 0 ;
                 end;                 
                 * at time = 1 the lagged average is x(0)  = x at baseline ;
                 * at time = 2 the lagged value would be  (x(0) + x(1) ) / 2 , different from baseline value  ; 
                 &&cov&i.._rcumavg_l1 = &&cov&i.._rcumavg ;
                 if &time = 0 then &&cov&i.._rcumavg =  &&cov&i ;
                 else &&cov&i.._rcumavg = ( &&cov&i.._rcumavg_l1 + &&cov&i ) / 2.0 ;
           %end;
           %else %if &type = sim %then %do;
               if &time = 0 then do ;
                  %if &lagged = 1 %then &&cov&i.._rcumavg_l1 = 0 ;;
                  %if &current = 1 %then  &&cov&i.._rcumavg = &&cov&i ;;
               end ;
               %if &current = 1 %then else &&cov&i.._rcumavg = ( &&cov&i.._rcumavg_l1 + &&cov&i ) / 2.0 ; ;
               
           %end;        

        %end;

       /* for user-defined functions of cov&i and any cov&j where j < i (or possibly lagged values of cov&h where h > i ) */
        %if &&cov&i.genmacro^=  %then %do;
            %&&cov&i.genmacro;
       %end;

%if &type = main %then %do;

/* GET THE Z VARIABLES FOR 0TYPE 4 COVARIATES, AND LOG-TRANSFORM THE NONZEROS. */
    %if &&cov&i.otype = 4 %then %do;
            /* IF COV&I NE 0, THEN ZCOV&I=1, IF COV&I=0, THEN ZCOV&I=0. */
            z&&cov&i = (&&cov&i ne 0);
            if z&&cov&i = 1 then l&&cov&i = log(&&cov&i);
            
            z&&cov&i.._l1 = (&&cov&i.._l1 ne 0);
            if z&&cov&i.._l1 = 1 then l&&cov&i.._l1 = log(&&cov&i.._l1);
            
            %if &&cov&i.skip = -1 %then %do;
                %if &&cov&i.ptype = lag2bin or &&cov&i.ptype = lag2qdc or &&cov&i.ptype = lag2zqdc
                    or &&cov&i.ptype = lag2cub or &&cov&i.ptype = lag2cat or &&cov&i.ptype = lag2spl 
                    or &&cov&i.ptype = lag3bin or &&cov&i.ptype = lag3qdc or &&cov&i.ptype = lag3zqdc
                    or &&cov&i.ptype = lag3cub or &&cov&i.ptype = lag3cat or &&cov&i.ptype = lag3spl
                    %then %do;
                    z&&cov&i.._l2 = (&&cov&i.._l2 ne 0);
                    if z&&cov&i.._l2 = 1 then l&&cov&i.._l2 = log(&&cov&i.._l2);
                    %end;

                %if &&cov&i.ptype = lag3bin or &&cov&i.ptype = lag3qdc or &&cov&i.ptype = lag3zqdc
                    or &&cov&i.ptype = lag3cub or &&cov&i.ptype = lag3cat or &&cov&i.ptype = lag3spl
                    %then %do;
                    z&&cov&i.._l3 = (&&cov&i.._l3 ne 0);
                    if z&&cov&i.._l3 = 1 then l&&cov&i.._l3 = log(&&cov&i.._l3);
                    %end;
                %end;
            
            %else %do;
                %maketi(z&&cov&i,&time,&time._l1,&&cov&i.skip, &interval);
                %maketi(z&&cov&i.._l1,&time._l1,&time._l2,&&cov&i.skip, &interval);
                %end;
                
            %end;


     
  %end;     
%end;
/* code for complete hitory variables used ineof models */

%local timeindex startindex stopindex  ;
%if &useeof = 1 %then %do ;
    %if &&cov&i.etype_part2 = &timepoints %then %let  startindex = 0 ;
	%else %if &&cov&i.etype_part2 = 0 %then %let startindex = %eval(&timepoints);
	%else %let startindex = %eval(&timepoints - &&cov&i.etype_part2 ) ;

	%let stopindex = %eval(&timepoints - 1);
    %if &&cov&i.etype = cumavg or &&cov&i.etype = cumavgcat or &&cov&i.etype = cumsum or &&cov&i.etype = cumsumcat %then %do;
	     %if &&cov&i.etype_part2 = &timepoints %then %do;
			%let startindex = 0;
            %let stopindex =  0;
		 %end;
		 %else %do;
		 	%let startindex = %eval(&timepoints - &&cov&i.etype_part2 ) ;
            %let stopindex = %eval(&timepoints - &&cov&i.etype_part2 ) ;
		 %end;
	%end;
	%if &&cov&i.etype = cumavgnew or &&cov&i.etype = cumavgcatnew or &&cov&i.etype = cumsumnew or &&cov&i.etype = cumsumcatnew %then %do;
	     		%if &&cov&i.etype_part2 = &timepoints %then %do;
					%let startindex = 0;
            		%let stopindex =  0;
			 	%end;
		 		%else %do;
		 			%let startindex = %eval(&timepoints - &&cov&i.etype_part2 - 1 ) ;
            		%let stopindex = %eval(&timepoints - &&cov&i.etype_part2 - 1 ) ;
		 		%end;
	%end;	
	%do timeindex = &startindex %to &stopindex ;
        %let timeindex_l1 = %eval(&timeindex - 1 );
        %if &&cov&i.etype = qdc or &&cov&i.etype = zqdc    %then %do;
            &&cov&i.._&timeindex._eofs = &&cov&i.._&timeindex._eof *&&cov&i.._&timeindex._eof ;  /* SQUARE */
			%if &&cov&i.etype = zqdc %then %do;
				/* IF COV&I NE 0, THEN ZCOV&I=1, IF COV&I=0, THEN ZCOV&I=0. */
				z&&cov&i.._&timeindex._eof = (&&cov&i.._&timeindex._eof ne 0);
				if z&&cov&i.._&timeindex._eof = 1 then l&&cov&i.._&timeindex._eof = log(&&cov&i.._&timeindex._eof);
			%end;	
        %end;

        %if &&cov&i.etype = cub  %then %do;
            &&cov&i.._&timeindex._eofs = &&cov&i.._&timeindex._eof *&&cov&i.._&timeindex._eof;             /* SQUARE */
            &&cov&i.._&timeindex._eofc = &&cov&i.._&timeindex._eof * &&cov&i.._&timeindex._eof * &&cov&i.._&timeindex._eof;     /* CUBE */
        %end;
            
        %if &&cov&i.etype = spl AND %numargs(&&cov&i.eknots) > 1   %then %do;
            %rcspline(&&cov&i.._&timeindex._eof ,&&cov&i.eknots);        /* SPLINE, SEE BELOW! */
        %end;
            
        %if &&cov&i.etype = cat AND %numargs(&&cov&i.eknots) > 1   %then %do;    
            %if &&cov&i.otype ne 5 %then %do;             /* MAKE CATEGORIES, SEE BELOW! */
                %makecat(&&cov&i.._&timeindex._eof , &&cov&i.eknots, &&cov&i.elev );
            %end;
			/*****/
            %if &&cov&i.otype = 5 %then %do;
                %do lev = 1 %to %eval(&&cov&i.lev - 1);
                     &&cov&i.._&timeindex._eof_&lev = (&&cov&i.._&timeindex._eof = &lev);
                %end;
           %end; 
          /****/ 
        %end;

/* IF TIMES WERE SKIPPED FOR THIS VARIABLE, WE NEED THE INTERACTION WITH TIME. */
 
  
      %if &&cov&i.etype = skpbin    %then %do;
           *%if &current = 1 %then  %maketi(&&cov&i.._&timeindex._eof , &timeindex, &timeindex_l1, &&cov&i.skip, &interval);
                  
      %end;

        %if &&cov&i.etype = skpqdc or &&cov&i.etype = skpzqdc   
               %then %do;
			 %if &&cov&i.etype = skpzqdc %then %do;
				/* IF COV&I NE 0, THEN ZCOV&I=1, IF COV&I=0, THEN ZCOV&I=0. */
				z&&cov&i.._&timeindex._eof = (&&cov&i.._&timeindex._eof ne 0);
				*if z&&cov&i.._&timeindex._eof = 1 then l&&cov&i.._&timeindex._eof = log(&&cov&i.._&timeindex._eof);
			%end;
            *%maketi(&&cov&i.._&timeindex._eof, &timeindex, &timeindex_l1, &&cov&i.skip, &interval);
            &&cov&i.._&timeindex._eofs = &&cov&i.._&timeindex._eof *&&cov&i.._&timeindex._eof ;;
            *%maketi(&&cov&i.._&timeindex._eofs , &timeindex, &timeindex_l1, &&cov&i.skip, &interval);
                     
       %end;


        %if &&cov&i.etype = skpcub   %then %do;
            *%maketi(&&cov&i.._&timeindex._eof , &timeindex, &timeindex_l1, &&cov&i.skip, &interval);
         
            
            &&cov&i.._&timeindex._eofs = &&cov&i.._&timeindex._eof * &&cov&i.._&timeindex._eof ;;
            &&cov&i.._&timeindex._eofc = &&cov&i.._&timeindex._eof*&&cov&i.._&timeindex._eof*&&cov&i.._&timeindex._eof ;;
            
            *%maketi(&&cov&i.._&timeindex._eofs, &timeindex, &timeindex_l1, &&cov&i.skip, &interval);
            *%maketi(&&cov&i.._&timeindex._eofc, &timeindex, &timeindex_l1, &&cov&i.skip, &interval);
           
       %end;

                        
        %if &&cov&i.etype = skpspl AND %numargs(&&cov&i.eknots) > 1   %then %do;
             *%maketi(&&cov&i.._&timeindex._eof,&timeindex,&timeindex_l1,&&cov&i.skip, &interval);
             %rcspline(&&cov&i.._&timeindex._eof ,&&cov&i.eknots);
             
            /***
            %do knot = 1 %to %eval(&&cov&i.lev - 2);
                %maketi(&&cov&i.._&timeindex._eof_spl&knot,&timeindex,&timeindex_l1,&&cov&i.skip, &interval);
            %end;
			 ***/
        %end;
            
        %if &&cov&i.etype = skpcat AND %numargs(&&cov&i.eknots) > 1   %then %do;

             %makecat(&&cov&i.._&timeindex._eof , &&cov&i.knots, &&cov&i.lev);
            /**
            %do lev = 1 %to %eval(&&cov&i.lev - 1);
                %maketi(&&cov&i.._&timeindex._eof_&lev ,&timeindex,&timeindex_l1,&&cov&i.skip, &interval);
              
            %end;
			 **/
        %end; /* SKPCAT */      
         

        %if ( &&cov&i.etype=cumavg |  &&cov&i.etype = cumavgcat | &&cov&i.etype = cumavgnew | &&cov&i.etype = cumavgcatnew | 
               &&cov&i.etype = cumsum | &&cov&i.etype = cumsumnew | &&cov&i.etype = cumsumcat | &&cov&i.etype = cumsumcatnew ) 
           %then %do;
            
		
		      %if &&cov&i.etype = cumsum or &&cov&i.etype = cumsumcat  %then %do;
				    %if &timeindex = 0 %then %do;
						&&cov&i.._cumsum_&timeindex._eof  = sum(of a&&cov&i{*});
					%end;
					%else %if &timeindex > 0 %then %do;
	                     mysum = 0 ;
						 do myi = &startindex to %eval(&timepoints - 1);
						     mysum = mysum + a&&cov&i [ myi ];
						  end;
						 &&cov&i.._cumsum_&timeindex._eof  = mysum  ;
						 DROP MYSUM MYI ;
					%end;
				%end ;
				%else %if &&cov&i.etype = cumsumnew or &&cov&i.etype = cumsumcatnew  %then %do;
               		 mysum = 0 ;
					 do myi = 0 to %eval(&startindex);
					     mysum = mysum + a&&cov&i [ myi ];
					  end;
					 &&cov&i.._cumsum_&timeindex._eof  = mysum  ;
					 DROP MYSUM MYI ;
			    %end;
                

               %if &&cov&i.etype = cumavg or &&cov&i.etype = cumavgcat %then %do;
				    %if &timeindex = 0 %then %do;
						&&cov&i.._cumavg_&timeindex._eof  = mean(of a&&cov&i{*});
					%end;
					%else %if &timeindex > 0 %then %do;
	                     mysum = 0 ;
						 do myi = &startindex to %eval(&timepoints - 1);
						     mysum = mysum + a&&cov&i [ myi ];
						  end;
						 &&cov&i.._cumavg_&timeindex._eof  = mysum / %eval(&timepoints  - &startindex );
						 DROP MYSUM MYI ;
					%end;
				%end ;
				%else %if &&cov&i.etype = cumavgnew or &&cov&i.etype = cumavgcatnew %then %do;
               		 mysum = 0 ;
					 do myi = 0 to %eval(&startindex);
					     mysum = mysum + a&&cov&i [ myi ];
					  end;
					 &&cov&i.._cumavg_&timeindex._eof  = mysum / %eval(&startindex + 1 );
					 DROP MYSUM MYI ;
			    %end;

             

				%if (&&cov&i.etype in cumavg cumavgnew )   AND  ( %numargs(&&cov&i.eknots) ge 3 ) /* AND &usespline=1 */  %then %do;
                    %rcspline(&&cov&i.._cumavg_&timeindex._eof , &&cov&i.eknots);
                   
                %end;
                %else %if (&&cov&i.etype in cumavgcat cumavgcatnew ) AND ( %numargs(&&cov&i.eknots) ge 3 )    %then %do ;                      
                     %makecat(&&cov&i.._cumavg_&timeindex._eof , &&cov&i.eknots, &&cov&i.elev);
                                                            
                 %end;
				 %else %if (&&cov&i.etype in cumsum cumsumnew )   AND ( %numargs(&&cov&i.eknots) ge 3 ) /* AND &usespline=1 */  %then %do;
                    %rcspline(&&cov&i.._cumsum_&timeindex._eof , &&cov&i.eknots);
                   
                %end;
				%else %if (&&cov&i.etype in cumsumcat cumsumcatnew ) AND ( %numargs(&&cov&i.eknots) ge 3 )    %then %do ;                      
                     %makecat(&&cov&i.._cumsum_&timeindex._eof , &&cov&i.eknots, &&cov&i.elev);
                                                            
                 %end;
                                                                                   
            
        %end;

        

       /* for user-defined functions of cov&i and any cov&j where j < i (or possibly lagged values of cov&h where h > i ) */
		/***
        %if &&cov&i.genmacro^=  %then %do;
            %&&cov&i.genmacro;
       %end;
	   ****/
    %end ; /* timeindexz */
	
%end ; /* useeof = 1  */




    
    %mend genpred;

/* THIS MACRO MAKES CATEGORIES. IF THE VARIABLE X IS ALREADY INTEGERS, BE CAREFUL ABOUT THE KNOTS. */
/* NOTE THAT IT LOOKS FOR X TO BE STRICTLY LESS THAN A KNOT TO DECIDE WHICH CATEGORY TO PUT IT IN. */
%MACRO MAKECAT(x,knots,levels);
    %local lev ;
    if &x = . then &x.gr = .;
    %do lev = 1 %to %eval(&levels - 1);
        else if  &x < %scan(&knots,&lev,%str( )) then &x.gr = &lev;
    %end;
    else &x.gr = &levels;
    %do lev = 1 %to %eval(&levels - 1);
        &x._&lev = (&x.gr = &lev);
        %end;
    drop &x.gr;
    %MEND MAKECAT;

/* TIME INTERACTION FOR WHEN THE DATA WASN'T COLLECTED IN THAT TIME PERIOD. */
/* NOTE: ASSUMED THAT THERE'S A CONSTANT INTERVAL BETWEEN MEASUREMENTS/QUESTIONNAIRES. */
/* ANOTHER ASSUMPTION: NOTHING IS SKIPPED MORE THAN 2 PERIODS IN A ROW. */



/* A new version of the maketi routine for arbitrary number of times to check for skipped times.
    if a starting time TIME1 is in the skipped list then at most n - 1 other times can be in the skipped list
    where n = the number of times in the skipped list. In the extreme case time1 is the value of the largest skipped time.
    Since the difference between time1 and time2 is always 1 all we need to check are the n-1 previous times prior to the 
    value of time1.
        
    To create the interaction term we only need to count the number of testimes that are also skip times
        0 1 2 3 4 ... TIME1 are checked working backwards in time to see if they are also skip times. Once we find one that
    is not a skip time there had to be a true measurement at that time and we will have found the number of consecutive times
    that were skipped. This will be _mycount_ . The interaction term is then x_ti = x * _mycount_ * interval .

    This provides the same when there were at most two consecutive times that could be skipped.
    ***/

%MACRO MAKETI(x,time1,time2,skipped,interval);

    %local nskip   ;

    %let nskip = %numargs(&skipped ) ;
     
    _testtime_ = &time1 ; 

    &x._ti = 0;
    
    * find how many periods are in the skipped list ;
    _mycount_ = 0 ;

    do _myind = 1 to &nskip ;
        if _testtime_ in (&skipped) then do ;
             _mycount_ = _mycount_ + 1 ;
             _testtime_ = _testtime_ - 1 ;
        end;
        else do ;
           leave ;
        end;  
    end;

    
   &x._ti = &x * _mycount_ * &interval ;
    
   drop _myind _testtime_ _mycount_  ;
  
    %MEND MAKETI;


/* CREATES VARIABLES NEEDED FOR RESTRICTED CUBIC SPLINE. */
%MACRO RCSPLINE(x,knots,norm=2);
 /* based on original rcspline macro with changes for use in gformula .
     knots are now included in a single list */
    /*MACRO RCSPLINE



   For a given variable named X and from 3-10 knot locations,
   generates SAS assignment statements to compute k-2 components
   of cubic spline function restricted to be linear before the
   first knot and after the last knot, where k is the number of
   knots given.  These component variables are named c1, c2, ...
   ck-2, where c is the first 7 letters of X.

   Usage:

   DATA; ....
   %RCSPLINE(x,knot1,knot2,...,norm=)   e.g. %RCSPLINE(x,-1.4,0,2,8)

        norm=0 : no normalization of constructed variables
        norm=1 : divide by cube of difference in last 2 knots
                 makes all variables unitless
        norm=2 : (default) divide by square of difference in outer knots
                 makes all variables in original units of x

   Reference:

   Devlin TF, Weeks BJ (1986): Spline functions for logistic regression
   modeling. Proc Eleventh Annual SAS Users Group International.
   Cary NC: SAS Institute, Inc., pp. 646-51.


   Author  : Frank E. Harrell Jr.
             Clinical Biostatistics, Duke University Medical Center
   Date    : 10 Apr 88
   Mod     : 22 Feb 91 - normalized as in S function rcspline.eval
             06 May 91 - added norm, with default= 22 Feb 91
             10 May 91 - fixed bug re precedence of <>

             28 Feb 08 - allowed for names of length 32 (roger logan )
                                                                      */
    %LOCAL a v7 b tk tk1 t k1 k2 iknot;
    
    %LET v7=&x; %IF %LENGTH(&v7)=32 %THEN %LET v7=%SUBSTR(&v7,1,31);

    %*Get no. knots, last knot, next to last knot;
    %let knotcount=0;
	%let knotcount = %numargs(&knots) ;
	%do iknot = 1 %to &knotcount ;
	     %let knot&iknot = %scan(&knots,&iknot,%str( ));
		 
	%end;
	/****
    %do %while (%scan(&knots,%eval(&knotcount+1),%str( )) ne %str());
        %let knotcount = %eval(&knotcount+1);
        %local knot&knotcount;
        %let knot&knotcount = %trim(%left(%qscan(&knots,&knotcount,%str( ))));
        %end;
    ***/
    %let a=&knotcount;
    %let k1=%eval(&knotcount-1);
    %let k2=%eval(&knotcount-2);

    %IF &a<3 %THEN %PUT ERROR: <3 KNOTS GIVEN.  NO SPLINE VARIABLES CREATED.;
    %ELSE %DO;
        %LET tk=&&knot&a;
        %LET tk1=&&knot&k1;
        DROP _kd_; _kd_=
            %IF &norm=0 %THEN 1;
        %ELSE %IF &norm=1 %THEN &tk - &tk1;
        %ELSE (&tk - &knot1)**.666666666666; ;
        %DO a=1 %TO &k2;
            %LET t=&&knot&a;
            &x._spl&a=max((&x-&t)/_kd_,0)**3+((&tk1-&t)*max((&x-&tk)/_kd_,0)**3
                -(&tk-&t)*max((&x-&tk1)/_kd_,0)**3)/(&tk-&tk1)%STR(;);
            %END;
        %END;
    %MEND rcspline;

/************************************************************************************************************************/


%macro simvar(start=1, stop=&ncov);  
    
    %*Generate random numbers for standard errors;
   %local i  lev;
   %if &outctype = conteofu4 %then U&outc = rand('uniform');;
   %if &outctype = cateofu %then %do;
       %do lev = 1 %to %eval(&outclev - 1);
                Uoutc_&lev=rand('uniform');
       %end;
   %end;
    %do i = &start %to &stop;
        %if &&usevisitp&i = 1 %then %do;            
             U&&cov&i.randomvisitp = rand('uniform');
        %end;
        %if &&cov&i.otype=1 or &&cov&i.otype=2  or  &&cov&i.otype=6  %then %do; /* include truncated normal here */
            U&&cov&i=rand('uniform');
        %end;
        %else %if &&cov&i.otype=3  or  &&cov&i.otype=7 %then %do;
            N&&cov&i= rand('normal');
        %end;
        %else %if &&cov&i.otype=4 %then %do;
            Uz&&cov&i=rand('uniform');
            N&&cov&i = rand('normal');
        %end;
        %else %if &&cov&i.otype=5 %then %do;
            %do lev = 1 %to %eval(&&cov&i.lev - 1);
                U&&cov&i.._&lev=rand('uniform');
            %end;
        %end;          
        %else %if &&cov&i.otype=-1 %then %do;
            Uclass&&cov&i =rand('uniform');
            N&&cov&i  = rand('normal');
       %end;
   %end;


    %* Generating time k covariates;

    
       
    %do i = &start %to &stop;        
        %*Generating Type 0 variables;
        %if &&cov&i.otype=0 %then %do;
            &&cov&i = &&cov&i  + &&cov&i.inc;
        %end;

        %*Carrying forward values if skipped year;
        /* since we use the same variable name for each time point, doing nothing is the same as 
           carrying forward a value, only need to look at non-skip times  */

        if   not( &time  in (&&cov&i.skip) ) then do;                                                             
            if &&cov&i.simstart > &time  then do; /* assign fixed value if simulation should start after time k*/                      
                &&cov&i = &&cov&i.else + 0;
                %if &intvisittype = 2 %then %do ;
					_&&&cov&i.._holder = &&cov&i ; * will hold the most recent simulated value. will always be in a non-skip time ;
				%end;
                %if &&usevisitp&i = 1 %then %do;
                    &&cov&i.randomvisitp = &&cov&i.visitpelse ;
                    if &&cov&i.randomvisitp = 1 then ts_last_&&cov&i = 0 ;
                    else ts_last_&&cov&i = ts_last_&&cov&i.._l1 + 1 ;
                %end;

                %if &&cov&i.otype=-1 %then %do; 
                    %&&cov&i.setblvar ;
                %end;
            end; 
            else do;          
                if %unquote(&&cov&i.wherenosim) then do;  
                    /* evaluate user defined macro  because of some condition other than its too early to start simulating*/
                    %if %bquote(&&cov&i.nosimelsemacro)^= %then  %&&cov&i.nosimelsemacro ;
                    
                end;
                else do; /*else simulate from the estimated model parameters*/   
                    %if &&usevisitp&i = 1 %then %do;                              
                        m&&cov&i.randomvisitp= 0 ;  
                        do j=1 to dim(abvisitp&&cov&i);
                            m&&cov&i.randomvisitp=sum(m&&cov&i.randomvisitp,abvisitp&&cov&i..(j)*as&&cov&i(j));
                        end;
                        p&&cov&i.randomvisitp.=1/(1+exp(-m&&cov&i.randomvisitp));
                        if U&&cov&i.randomvisitp <= p&&cov&i.randomvisitp    then &&cov&i.randomvisitp=1;
                        if U&&cov&i.randomvisitp >  p&&cov&i.randomvisitp >. then &&cov&i.randomvisitp=0;

                        if ts_last_&&cov&i.._l1  = &&cov&i.visitpmaxgap   then do;
                            &&cov&i.randomvisitp=1;                 
                        end;

                        if  &&cov&i.randomvisitp = 1 then do;
                            ts_last_&&cov&i = 0 ;
                        end;                                                        
                        else if &&cov&i.randomvisitp = 0  then do ;
                            ts_last_&&cov&i = ts_last_&&cov&i.._l1 + 1 ;  
							%if &intvisittype = 2 %then %do;
								&&cov&i =  _&&cov&i.._holder ; * uses the previous simulated value instead of a potential inervened on value ; 
							%end; 
                        end ;

                    %end;

                    %if &&usevisitp&i = 1 %then %do;                
                        if  &&cov&i.randomvisitp = 1 then do;
                    %end ;

                    %*Generating Type 1 variables;
                        %if &&cov&i.otype=1 %then %do;
                            m&&cov&i = 0 ;
                            do j=1 to dim(ab&&cov&i);
                                m&&cov&i=sum(m&&cov&i,ab&&cov&i(j)*as&&cov&i(j));
                            end;
                            p&&cov&i=1/(1+exp(-m&&cov&i));
                            if U&&cov&i <=p&&cov&i  then &&cov&i =1;
                            if U&&cov&i >p&&cov&i >. then &&cov&i=0;
                        %end;

                        %*Generating Type 2 variables;
                        %else %if &&cov&i.otype=2 %then %do;
                            m&&cov&i  = 0 ;
                            do j=1 to dim(ab&&cov&i..);                                
                                m&&cov&i =sum(m&&cov&i ,ab&&cov&i..(j)*as&&cov&i(j));
                            end; 
                            if &&cov&i =1 then &&cov&i =1;
                            if &&cov&i =0 then  do;
                                p&&cov&i =1/(1+exp(-m&&cov&i));
                                if U&&cov&i <=p&&cov&i   then &&cov&i =1;
                                if U&&cov&i >p&&cov&i >. then &&cov&i =0;
                            end;
                        %end;

                        %*Generating Type 3 variables;
                        %else %if &&cov&i.otype=3 %then %do;
                            m&&cov&i = 0 ;
                            do j=1 to dim(ab&&cov&i);
                                m&&cov&i =sum(m&&cov&i ,ab&&cov&i(j)*as&&cov&i (j));
                            end; 
                            &&cov&i = m&&cov&i  + N&&cov&i *se&&cov&i;
							%if &sim_trunc = 1 %then %do;
	                            if &&cov&i  < &&cov&i.._min then do;
	                                &&cov&i = &&cov&i.._min; &&cov&i.._limit = &&cov&i.._limit +1;
	                            end;

	                            if &&cov&i > &&cov&i.._max then do;
	                                &&cov&i  = &&cov&i.._max; &&cov&i.._limit = &&cov&i.._limit +1;
	                            end;
							%end; 
                        %end;

                        %*Generating Type 4 variables;
                        %else %if &&cov&i.otype=4 %then %do;
                            %if %bquote(&&cov&i.class)^= %then %do;
                                mz&&cov&i = 0 ;
                                if &&cov&i = 0 then do;
                                    do j=1 to dim(abz0&&cov&i);
                                        mz&&cov&i =sum(mz&&cov&i,abz0&&cov&i..(j)*as&&cov&i(j));
                                    end;
                                end;
                                if &&cov&i ^= 0 then do;
                                    mz&&cov&i = 0 ;
                                  do j=1 to dim(abz1&&cov&i);
                                      mz&&cov&i=sum(mz&&cov&i,abz1&&cov&i..(j)*as&&cov&i(j));
                                  end;
                                end;
                            %end;
                            %else %do;
                                mz&&cov&i = 0 ;
                                do j=1 to dim(abz&&cov&i);
                                    mz&&cov&i=sum(mz&&cov&i,abz&&cov&i..(j)*as&&cov&i(j));
                                end;
                            %end;
                            pz&&cov&i=1/(1+exp(-mz&&cov&i));
                            if Uz&&cov&i<=pz&&cov&i   then z&&cov&i=1;
                            if Uz&&cov&i>pz&&cov&i>. then z&&cov&i=0;

                            if z&&cov&i=0 then &&cov&i=0;
                            else do;
                                m&&cov&i = 0 ;
                                do j=1 to dim(ab&&cov&i..);
                                    m&&cov&i=sum(m&&cov&i,ab&&cov&i..(j)*as&&cov&i(j));
                                end;
                                &&cov&i=exp(m&&cov&i+N&&cov&i*se&&cov&i);

                                /*  move boundary checks to inside interval since min is now for non-zero values */
                                if &&cov&i < &&cov&i.._min then do;
                                    &&cov&i = &&cov&i.._min; &&cov&i.._limit = &&cov&i.._limit +1;
                                end;

                                if &&cov&i > &&cov&i.._max then do;
                                    &&cov&i = &&cov&i.._max; &&cov&i.._limit = &&cov&i.._limit +1;
                                end;     
                            end;


                        %end;

                        %*Generating Type 5 variables;
                        %else %if &&cov&i.otype=5 %then %do;
                            &&cov&i = . ;
                            %do lev = 1 %to %eval(&&cov&i.lev - 1);
                                if &&cov&i = . then do;
                                    m&&cov&i.._lev = 0 ;
                                    do j=1 to dim(ab&&cov&i.._&lev);
                                        m&&cov&i.._&lev=sum(m&&cov&i.._&lev,ab&&cov&i.._&lev(j)*as&&cov&i(j));
                                    end;
                                    p&&cov&i.._&lev=1/(1+exp(-m&&cov&i.._&lev));
                                    if U&&cov&i.._&lev<=p&&cov&i.._&lev  then &&cov&i=&lev;
                                end;
                           %end;

                           if &&cov&i = . then &&cov&i = &&cov&i.lev;

                        %end; /* type 5 */
                        %else %if &&cov&i.otype=6 %then %do; /* true truncated normal using estimates from proc qlim */
                            m&&cov&i = 0;
                            mubar = . ;
                            do j=1 to dim(ab&&cov&i);
                                m&&cov&i =sum(m&&cov&i ,ab&&cov&i(j)*as&&cov&i(j));
                            end; 

                            Fa =  probnorm((&&cov&i.._min - m&&cov&i) / se&&cov&i ) ;
                            Fb =  probnorm((&&cov&i.._max - m&&cov&i) / se&&cov&i );

                            mubar = (1-U&&cov&i ) *Fa +  U&&cov&i * Fb ;

                            &&cov&i= m&&cov&i +  se&&cov&i  * probit(mubar) ;

                       
                        %end;

                        %else %if &&cov&i.otype=7 %then %do; /* true tobit using estimates from proc qlim */
                            m&&cov&i = 0 ;
                            do j=1 to dim(ab&&cov&i);
                                m&&cov&i=sum(m&&cov&i ,ab&&cov&i(j)*as&&cov&i(j));
                            end; 
                            &&cov&i= m&&cov&i+ N&&cov&i*se&&cov&i;

                            if &&cov&i  < &&cov&i.._min then do;
                                &&cov&i = &&cov&i.._min; &&cov&i.._limit = &&cov&i.._limit +1;
                            end;

                            if &&cov&i > &&cov&i.._max then do;
                                &&cov&i = &&cov&i.._max; &&cov&i.._limit = &&cov&i.._limit +1;
                            end;            
                        %end;



                        %*Generating Type -1 override variables; 
                        %else %if &&cov&i.otype=-1 %then %do;
                            %&&cov&i.simusermacro ;
                        %end; /*type -1*/
                        %if &&usevisitp&i = 1 %then %do;
                            %if &enforcegap > 0 %then %do;
                                if ts_last_&&cov&i = &enforcegap AND &&cov&i.randomvisitp = 0 then do;
                                    mygood = 0 ; 
                                end;
                            %end;
							%if &intvisittype = 2  %then %do ;
				 				_&&&cov&i.._holder = &&cov&i ; * will hold the most recent simulated value. will always be in a non-skip time ;
							%end;
                    end ; /* this end is for the do loop that is conditional on the usevisitp macro variable, need an end here */
                        %end;

                end ;
            end;   

            s&&cov&i [ &time] = &&cov&i ;
			%if &intvisittype = 2  AND &&usevisitp&i = 0 %then %do ;
			      * when there is a visitprocess the value of &&cov&i at this point could be from a carried forward value when this timepoint there is no
			        visit. Only evaluate this definition when there is not a visit process. The code for the visit process situation is defined above. ;
				 _&&&cov&i.._holder = &&cov&i ; * will hold the most recent simulated value. will always be in a non-skip time ;
			%end;
            %if &&usevisitp&i = 1 %then  s&&cov&i.randomvisitp[&time] = &&cov&i.randomvisitp ;;
            %*Generating derived variables (categories, etc);

            %genpred(sim,lagtype=1);   
        end;
        else do ;
		   /* with a skip type variable there is no visit process. We do not need to condition on &&usevisitp&i */

		    %if &intvisittype=1 %then %do;
            	&&cov&i = &&cov&i.._l1 ; *original method to carry forward during skip time. This will use the intervened on value ;
			%end;
			%else %if &intvisittype = 2 %then %do;
				&&cov&i = _&&cov&i.._holder ; * this uses the previous simulated value need to keep track of this value even when changed due to intervention ;
			%end;
            s&&cov&i [&time] = &&cov&i ;
            /* with a skip type variable there is no rndom visit process , so we do not need to update the next line 
             %if &&usevisitp&i = 1 %then  s&&cov&i.randomvisitp[&time] = &&cov&i.randomvisitp ;; */
            %*Generating derived variables (categories, etc);

            %genpred(sim,lagtype=1);  
        end;
          
    %end; /* i */

%mend simvar;







 
%macro obscuminc(data=,time=,timepoints=,event=,compevent=compevent );    
 
        proc sql noprint ;
            create table summed  as
            select   &time  , mean(&outc) as proboutc ,                      
                      %if %bquote(&compevent)^= AND &compevent_cens = 0 %then mean(&compevent) as probcompevent  ;
                      %else 0 as probcompevent ;                                     
           from  &data(keep = &time  &outc &compevent )
        
           group by &time 
           order by &time 
           ;      
       quit;

       
        data cuminc;
        set summed;
        keep cuminc;
        by &time;
        retain cumsurv 1 cuminc 0;   
        inc = cumsurv * proboutc * (1.0 - probcompevent) ;         
        cuminc = cuminc + inc;
        surv = (1.0 - proboutc) * (1.0 - probcompevent) ;         
        cumsurv = cumsurv * surv;
        if _N_ = &timepoints then output; /* this is based on row number and not time value so there is no -1 here */
        run;

        proc datasets library=work nolist;
        delete    summed;
        run;
    
    %mend;



 

%macro simcuminc;  
          
    %local myi ;

    cumsurv = 1;
    cuminc = 0;
    cumincdead=0;
    array cumincr{0:%eval(&timepoints - 1)} ;
    array cumsurvr{0:%eval(&timepoints - 1)} ;   
    array cumcompeventr{0:%eval(&timepoints - 1)}  ;
    if calchazard = 0 then do;
         do time = 0 to %eval(&timepoints - 1) ;
		    /* when placing this code at the top before cumsurv has been updated has the same
		       effect as using time - 1 in the array for cumsurvr or uses the previous value of 
		       cumsurv */
         	%if &intno = 0 /**** AND %bquote(&censor) ^= ****/ %then %do;
				%do myi = 1 %to &ncov ;				   
					ncs&&cov&myi [ time] = s&&cov&myi [ time] * cumsurv ;
					%if &&usevisitp&myi = 1 %then ncs&&cov&myi.randomvisitp [ time ] = s&&cov&myi.randomvisitp [ time] * cumsurv  ;;
				%end;
			%end;

             inc = cumsurv*s&outc[time]*(1-scompevent[time]);  
             incdead = cumsurv * scompevent[time];
                 
             cuminc = cuminc + inc;
             cumincr[time] = cuminc ;

             cumincdead = cumincdead + incdead;      
             cumcompeventr[time] = cumincdead ;

             surv = (1-s&outc[time])*(1-scompevent[time]);
             cumsurv = cumsurv*surv;
             cumsurvr[time] = cumsurv ; 

		
 
         end;
     end; 
     drop time inc surv  ;    
%mend;



/* THESE LAST TWO MACROS ARE FOR CONTINUOUS OR DICHOTOMOUS END OF FOLLOW-UP ANALYSIS. */
%macro obstotcont(data=,time=,timepoints=,outc=,compevent=compevent );
   /* Note: this is a sort of crude average outcome at end of follow-up; this is not a 
    "fair comparison" for the simulated results. */
    
    proc sort data=&data(keep=&time &outc /* &compevent   */) out=limit;
        by &time;
    run;
    %if &outctype ne cateofu %then %do;
    data summed;
        set limit;
        by &time;
        keep &time n e;
        retain n e 0;   
        if &outc ne . then do;
            n = n+1;
            e = e+&outc;
        end;
   
       if last.&time then do; 
          output;
          n=0;
          e=0;
       end;
   
    run;

    
   
    data cuminc;
        set summed;
        keep cuminc;        
        if _N_=&timepoints then do;
           cuminc=(e/n);
           output;       
       end;
    run;

	 proc datasets library=work nolist;
        delete limit summed;
    run;
    %end;
	%else %do;
	   /* this code assumes there is no correction for ipw censoring variable */
		proc sql noprint  ;
		create table summed0 as 
		select &outc, count(*) as proportion  from limit where &time = %eval(&timepoints - 1)  group by &outc ;
		select sum(proportion) as nn into :total_num from summed0 ;
		create table proportions0 as select &outc label="&outc level" , round(proportion * (100/&total_num), 0.01) label="Proportion (%)" as proportion from summed0 ;
		quit  ;
        	
        proc datasets library=work nolist;
        delete limit summed0 ;
        quit ;
	%end;
   

    %mend;



%macro simtotcont;
    %local k ;
    data summed;
        set simulated;
        
        keep &time n e r c;
        indcumsurv=1;
   
        %do k = 3 %to %eval(2+&timepoints);

            retain r_&k n_&k e_&k c_&k 0;

          c_&k = c_&k + scensl&k; /*number of loss to followup: if scensl&k=1 then pchd&k=.*/ 

              if p&outc.&k ne . and pcompevent&k ne . then do; 
                                          /*note: below is conditional on not simulated as lost*/
           
                surv = (1 - (pcompevent&k)); 
                indcumsurv = indcumsurv * surv;

            r_&k = r_&k + pcompevent&k; /*expected number of competing cause deaths*/
                n_&k = n_&k + indcumsurv; /*expected number still in dataset: 
                         if scensl&k=1 then we don't add anything to n_&k*/
    
            e_&k = e_&k + (p&outc.&k * indcumsurv); /*total weighted sum of continuous outcome*/
                          /*(weighted by individual cumulative probability of survival)*/
                end;
       
            &time = &k;

            n = n_&k;
            e = e_&k;
        r = r_&k;
        c = c_&k;
       
       
            if _N_ = &nsimul then output;

            %end;

    run;

 
        
    data cuminc;
        set summed;
        keep cuminc _sample_;

        _sample_ = &bsample;

        if _N_ = &timepoints then do;
       cuminc=(e/n) ;  /*weighted average--the weights are already in e and n*/
       output;
       end;
   
    run;
   

    data simulated;
        merge simulated cuminc;
        by _sample_;
    run;

    %mend simtotcont;
/************************************************************************************************************/


/* ADDED BELOW: SUB-MACROS TO CLEAN UP THE RESULTS SECTION. */

%macro rescaleround;
   
%if &outctype=bineofu or &outctype=binsurv %then %do;
   
  /* FOR BINARY OUTCOME, MULTIPLY RISKS BY 100 TO GET % AND ROUND OFF TO TWO DECIMAL PLACES.*/
   obsp= round(&obsp*10000)/100;   
   call symput('obspm',obsp);

   &pd = round(&pd*10000)/100;

   &pd._mean = round(&pd._mean*10000)/100;
   &pd._std  = round(&pd._std*10000)/100;
   
   &pd._ulim95 = round(&pd._ulim95*10000)/100;
   &pd._llim95 = round(&pd._llim95*10000)/100;
   
   RD = round(RD*10000)/100;

   RD_mean = round(RD_mean*10000)/100;
   RD_std  = round(RD_std*10000)/100;

   
   RD_ulim95 = round(RD_ulim95*10000)/100;
   RD_llim95 = round(RD_llim95*10000)/100;
   

   NNT = round(NNT);

   NNT_mean = round(NNT_mean);
   NNT_std  = round(NNT_std);

   

   
      NNT_ulim95 = round(NNT_ulim95);
      NNT_llim95 = round(NNT_llim95);
   
   
%end;
%else %if &outctype=conteofu or &outctype=conteofu2 or &outctype = conteofu3 or &outctype=conteofu4  %then %do;
   
   obsp= round(&obsp*100)/100;
   call symput('obspm',obsp);
   s&outc = round(s&outc*100)/100;

   s&outc._mean = round(s&outc._mean*100)/100;
   s&outc._std  = round(s&outc._std*100)/100;
   
   s&outc._ulim95 = round(s&outc._ulim95*100)/100;
   s&outc._llim95 = round(s&outc._llim95*100)/100;
   

   RD = round(RD*100)/100;

   RD_mean = round(RD_mean*100)/100;
   RD_std  = round(RD_std*100)/100;
   
   if int^=0 then do;
      RD_ulim95 = round(RD_ulim95*100)/100;
      RD_llim95 = round(RD_llim95*100)/100;
      end;
   
%end;

/* in both cases do the following: */
   RR = round(RR*100)/100;

      RR_ulim95 = round(RR_ulim95*100)/100;
      RR_llim95 = round(RR_llim95*100)/100;
   
   intervened = round(intervened*10000)/100;
   averinterv = round(averinterv*10000)/100;
   
%mend rescaleround;

%macro labels;
  
   label RR_ulim95 = 'Upper limit 95% CI';
   label RR_llim95 = 'Lower limit 95% CI';
   label RD_ulim95 = 'Upper limit 95% CI';
   label RD_llim95 = 'Lower limit 95% CI';
   label intervened = '% Intervened On';
   label averinterv = 'Aver % Intervened On';
   label int       = 'Interv.';
   label int2      = 'Description';
   
%if &outctype=binsurv  %then %do;
 	label pD_ulim95 = 'Upper limit 95% CI';
   label pD_llim95 = 'Lower limit 95% CI';
   label pD        = 'Risk (%)';
   label pD_std    = 'Bootstrap Risk SE';
   label pD_mean   = 'Bootstrap Risk Mean';
   label rr        = 'Risk ratio';
   label RD        = 'Risk difference';
   label NNT       = '# Needed to Treat';
   label NNT_ulim95 = 'Upper limit 95% CI';
   label NNT_llim95 = 'Lower limit 95% CI';
   %end;

 
%else %if &outctype=bineofu %then %do;
 label s&outc._ulim95 = 'Upper limit 95% CI';
   label s&outc._llim95 = 'Lower limit 95% CI';
   label s&outc        = 'Proportion';
   label s&outc._std    = 'Bootstrap Proportion SE';
   label s&outc._mean   = 'Bootstrap Proportion Mean';
   label rr        = 'Ratio';
   label RD        = 'Difference';
   label NNT       = '# Needed to Treat';
   label NNT_ulim95 = 'Upper limit 95% CI';
   label NNT_llim95 = 'Lower limit 95% CI';
   %end;


   
%mend labels;



 
%macro construct_graphs (time = period , 
                         outcome= outcome ,
                         compevent = ,
                         outctype=binsurv ,
                        
                        
                        simsurv = ,
                        obssurv= ,
                        covmean=,
                         sixgraphs=0 ,
                         gfilename=fourgraphs.pdf ,
                         title1=,title2=,title3=,titledata=,tsize=0.75 ,
                         frombootstrap = 1) ;
 





/****
 
 proc greplay igout = GSEG nofs ;
 delete _all_ ;
 run;
 quit;
*****/
 %local nperiods mintime mydate mywork anydoubles  idouble word simlist i ivar vartmp ngraph mycount myextra mypair graphlist graph1
        graph2 graph3 graph4 graph5 graph6 ;


%let mydate = &sysdate9._%scan(&systime,1,':')_%scan(&systime,2,':') ;
%let mywork = %sysfunc(getoption(work)) ;       
       filename mygraphs "&mywork./graphs.pdf";
   
    proc sql noprint ;
    select count(&time) as mtime into :nperiods from &covmean ;
    select min(&time)   as mintime into :mintime from &covmean ;      
    quit;

    %let nperiods = %sysfunc(compress(&nperiods)) ;
    %let mintime = %sysfunc(compress(&mintime)) ;

 goptions reset=all noborder device=png gsfname=mygraphs  ;
 proc greplay tc=work.tempcat nofs;
 tdef newtemp des="my four panel template"
     1/llx=0   lly=52
       ulx=0   uly=90
       urx=48  ury=90
       lrx=48  lry=52
        

     2/llx=0   lly=10
       ulx=0   uly=48
       urx=48  ury=48
       lrx=48  lry=10
       

     3/llx=52 lly=52
       ulx=52 uly=90
       urx=100 ury=90
       lrx=100 lry=52
       

     4/llx=52 lly=10
       ulx=52 uly=48
       urx=100 ury=48
       lrx=100 lry=10
        
 
        ;
        tdef newtemp2 des="my five panel template"
     1/llx=0   lly=52
       ulx=0   uly=90
       urx=48  ury=90
       lrx=48  lry=52
        

     2/llx=0   lly=10
       ulx=0   uly=48
       urx=48  ury=48
       lrx=48  lry=10
       

     3/llx=52 lly=52
       ulx=52 uly=90
       urx=100 ury=90
       lrx=100 lry=52
       

     4/llx=52 lly=10
       ulx=52 uly=48
       urx=100 ury=48
       lrx=100 lry=10
      
      5/llx=0   lly = 0  
        ulx=0   uly=100
        urx=100 ury=100
        lrx=100 lry=0
        ;
tdef newtemp3 des="my six panel template"
     1/llx=10   lly=68
       ulx=10   uly=98
       urx=40  ury=98
       lrx=40  lry=68
        

     2/llx=50   lly=68 
       ulx=50   uly=98
       urx=80   ury=98
       lrx=80   lry=68
       

     3/llx=10 lly=34
       ulx=10 uly=64
       urx=40 ury=64
       lrx=40 lry=34
       

     4/llx=50 lly=34
       ulx=50 uly=64
       urx=80 ury=64
       lrx=80 lry=34
        
     5/llx=10 lly=0 
       ulx=10 uly=30
       urx=40 ury=30
       lrx=40 lry=0
       

     6/llx=50 lly=0
       ulx=50 uly=30
       urx=80 ury=30
       lrx=80 lry=0
        ;
  template newtemp;
 list template;
 quit;


 %if %bquote(&title1)= %then %do;

  %let period = ;
  %if &frombootstrap = 0 %then %let period= . ;
  	    %if %bquote(&censor) ^= %then %do;
    		title1 height= &tsize  'Left column: IP-weighted natrual course estimates (solid line), parameteric g-formula estimated natural 
                                      course estimates (dotted line) by follow-up period';
        	title2 height= &tsize   "Right column: differences between IP-weighted and parametric g-formula natural course estimates (solid lines) ";
		%end;
		%else %do;
        	title1 height= &tsize  'Left column: observed (solid line), natural course (dotted line) estimates by follow-up period.';
        	title2 height= &tsize   "Right column: differences between observed and natural course estimates (solid lines)&period  ";
		%end;
        %if &frombootstrap = 1 %then title3 height=&tsize  'and 95% pointwise confidence intervals (dotted lines).';;
%end;
%else %do;
        title1 height= &tsize &title1 ;
        title2 height= &tsize &title2 ;
        title3 height= &tsize &title3 ;
%end; 

proc gslide gout=work.gseg ;
run;
quit;

title ;



 proc contents data = &covmean out = _cont (keep = name)   noprint ;
 run ;


%let anydoubles = 0 ; /* default value since it is possible that _cont2 has no observations */

data _cont2 ;
set _cont  ( where = ( substr(name,1,1)='s' 
  and index(name,'_ub') = 0 and index(name,'_lb')= 0 and index(name,'_stddev')=0 and index(name,'_diff')=0 )) ;
  if index(name,'ss') = 1 then doubles = 1 ;
  else doubles = 0 ;

  run;



  proc means data = _cont2 noprint ;
  var doubles ;
  output out = _anydouble (keep = sum)  sum = sum ;
  run;

  data _null_ ;
  set _anydouble ;
  call symput('anydoubles',compress(sum)) ;
  run;
 
    %if &anydoubles > 0 %then %do;
       data _doubles_;
       set _cont2 ;
       if doubles = 1 ;
       run;

       proc sql noprint ;
       select name into :doubles separated by ' ' from _doubles_ ;
       quit;

       
       %do idouble = 1 %to &anydoubles ;
            %let word = %scan(&doubles,&idouble);
            

            data _cont2 ;
            set _cont2 ;
            test = "&word" ;
            ind = index(compress(test),compress(name)) ;
            if ind = 2 then delete ;
            drop ind test ; 
            run;

           
         
           
       %end;
    
    %end;

 title ;
    
   
 data _cont2 ;
 set _cont2 ;
 drop doubles ;
 name = substr(name,2);
 if compress(name) = compress("&time") then delete ;
 run;

 %let simlist =  ;
 proc sql noprint ;
 select name into :simlist separated by ' ' from _cont2 ;
 quit;

   
   
       
   



    %if &outctype = binsurv %then %do;

        data _datass ;
        set &simsurv ;
        run;

        data _dataos ;
        set &obssurv ;
        run;


   


        * natural course risk ;
        data _estncrsk ;
        set _datass ( keep = risk1 - risk&nperiods _sample_ int where = (_sample_ = 0 and int = 0));
        array risk{&nperiods} ;
        &time = 0 ;
        r = 0 ;
        output ;
        do &time = 1 to &nperiods;
           r = risk[&time];
           output;
        end;
        keep &time r ;
        run;


        data _estncrsk ;
        set _estncrsk ;
        rename r = ncrisk ;
        run;


        %if &frombootstrap = 1 %then %do;
            data _bsss ;
            set _datass ( keep = risk1 - risk&nperiods _sample_ int where = (_sample_ > 0 and int = 0) )  ;
            drop int ;
            run;

            data _ncss ;
            set _dataos (keep = obrisk&mintime - obrisk%eval(&mintime + &nperiods - 1) _sample_ where = (_sample_ > 0));
            run;

            data _bsres ;
            merge _bsss _ncss ;
            by _sample_ ;
            array obrisk{&nperiods} obrisk&mintime - obrisk%eval(&mintime+&nperiods-1) ;
            array risk{&nperiods} risk1 - risk&nperiods ;

            array rdiff{&nperiods} ;
            do &time = 1 to &nperiods ;
               rdiff[&time] = obrisk[&time] - risk[&time];
            end;
            run;

        

        
           proc univariate data = _bsres  noprint  ;
           var  rdiff1 - rdiff&nperiods  ;


           output out = _bsstd2  pctlpre =  rdiff1 - rdiff&nperiods
                pctlname = _pct025 _pct975 
                pctlpts = 2.5 97.5 ;
            run;
         
            data _bsstd2 ;
            set _bsstd2 ;
            %do i = 1 %to &nperiods ;
                &time = &i ;
                lb = rdiff&i._pct025 ;
                ub = rdiff&i._pct975 ;
                output ;
            %end;
            keep &time lb ub ;
            run;

        %end;

        data _estobsrisk ;
        set _dataos (keep = obrisk&mintime - obrisk%eval(&mintime+&nperiods-1) _sample_ where = (_sample_ = 0));
        array obrisk{&nperiods} obrisk&mintime - obrisk%eval(&mintime+&nperiods-1) ;
         &time = 0 ;
        obs = 0 ;
        output ;
        do &time = 1 to &nperiods ;
            obs = obrisk[&time] ;
            output ;
        end;
        keep &time obs ;
        run ;

        data _estobsrisk ;
        set _estobsrisk ;
        rename obs = obrisk ;
        run;




        data _both ;
        merge _estncrsk _estobsrisk  %if &frombootstrap = 1 %then _bsstd2 ; ;
        by &time ;
        riskdiff = obrisk - ncrisk  ;
        %if &frombootstrap = 1 %then %do;
            if &time = 0 then do ;
               lb = 0 ;
               ub = 0 ;
            end;
        %end;
       
        run ;


    %end ; /* end of binsurv */
    %else %do; /* all other outcome types, may not work for all types */
       

	    %if &outctype ne cateofu %then %do;
	        data _both1_ ;
	        merge &obssurv &simsurv ; * for continuous at eof type outcome which use the mean ;
	        by _sample_ ;
	        **mean_diff = &outcome - ncmean ;
			mean_diff = &outcome - s&outcome  ;
	        run;

	        data _both10_ ;
	        set _both1_ (where = (_sample_ = 0 and int = 0));
	        run;

	        %if &frombootstrap = 1 %then %do ;
	            proc univariate data = _both1_ (where = (_sample_ > 0 ))  noprint  ;
	            var  mean_diff  ;
	            output out = _meanstd  pctlpre =  meandiff  
	                pctlname = _pct025 _pct975 
	                pctlpts = 2.5 97.5 ;
	            run;

	            data _graph2 ;
	            merge _both10_   _meanstd ;
	            lb = meandiff_pct025 ;
	            ub = meandiff_pct975 ;
	            run;

	            data _null_ ;
	            set _graph2 ;
	            call symput('lb',trim(left(round(lb,0.001))));
	            call symput('ub',trim(left(round(ub,0.001))));
	            call symput('meandiff',trim(left(round(mean_diff,0.001))));
	            run;


	            proc gslide name="mdiff";         
	            note height=3
	            justify=right 
	            color="black" 
	            "Mean difference &meandiff (&lb , &ub )" ;        
	            run;
	            quit;
	 
	       %end;

 
	       data _null_ ;
	       set _both10_ ;
	       call symput('obsmean',trim(left(round(&outcome,0.001))));
	       call symput('simmean',trim(left(round(s&outcome,0.001))));
	       run;

	 /* IP-weighted natrual course estimates (solid line), parameteric g-formula estimated natural  */

		   %local titletmp ;

		    %if %bquote(&censor)= %then %let titletmp = Observed nc ; 
			%else %let titletmp = Observed IP-weighted nc ;

	       proc gslide name="mean";   
	       note height=3
	       justify=left 
	       color="black" 
	       "&titletmp &outcome : &obsmean"  ;
		   note height = 3 color="black"
	       justify=left 
	       "Parametric g-formula estimated nc &outcome : &simmean";
	       run;
	       quit;
	 

       
  	   %end ;
        





    %end;
  /* for competing risk graphs */


    %if %bquote(&compevent)^= AND &compevent_cens = 0  %then %do;
        data _estobssurv ;
        set _dataos  ;
        array obrisk{&nperiods} obrisk&mintime - obrisk%eval(&mintime+&nperiods-1)  ;
        array osurv{&nperiods} obsurv&mintime - obsurv%eval(&mintime+&nperiods-1) ;
        array death{&nperiods} estobsdeath1 - estobsdeath&nperiods  ;
        do i = 1 to &nperiods ;
           death[i] = 1 - obrisk[i] - osurv[i];
        end;
        keep _sample_  estobsdeath1 - estobsdeath&nperiods  ;
        run;


        data _estncdeath ;
        set _datass (keep = _sample_ compevent1 - compevent&nperiods int where = (   int = 0 ));
        drop  int ;
        run;



        %if &frombootstrap = 1 %then %do;
            data _forgraph3  ;
            merge _estobssurv   _estncdeath  ;
            by _sample_ ;
            array death{&nperiods} estobsdeath1 - estobsdeath&nperiods  ;
            array cdeath{&nperiods} compevent1 - compevent&nperiods ;
            array riskdiff{&nperiods} ;

            if _sample_ > 0 ;
            do &time = 1 to &nperiods ;     
               riskdiff[&time] = death[&time] - cdeath[&time] ;      
            end;
            keep _sample_ riskdiff1 - riskdiff&nperiods ;
            run;

        
 

       
            proc univariate data = _forgraph3  noprint  ;
            var  riskdiff1 - riskdiff&nperiods  ;
            output out = _sdgraph4  pctlpre =  riskdiff1 - riskdiff&nperiods
                 pctlname = _pct025 _pct975 
                 pctlpts = 2.5 97.5 ;
            run;
     


            data _sdgraph4 ;
            set _sdgraph4 ;
            %do i = 1 %to &nperiods ;
               rename riskdiff&i._pct025 = pctlb&i riskdiff&i._pct975 = pctub&i ; 
            %end;
            run;

        %end ;


        data _graph3  ;
        merge _estobssurv(where = (_sample_ = 0)) _estncdeath(where = (_sample_ = 0))  
             %if &frombootstrap = 1 %then _sdgraph4 ; ;
        %if &frombootstrap = 1 %then %do ;
            array  pctlb{&nperiods} ;
            array  pctub{&nperiods} ;
        %end ;
        array death{&nperiods} estobsdeath1 - estobsdeath&nperiods  ;
        array cdeath{&nperiods} compevent1 - compevent&nperiods ;
        &time = 0 ;
        obs = 0 ;
        sim = 0 ;
        riskdiff = 0 ;
        lb = 0 ;
        ub = 0 ;
        output ;
        do &time = 1 to &nperiods ;
           obs = death[&time] ;
           sim = cdeath[&time] ;
           riskdiff = obs - sim ;
           
           %if &frombootstrap = 1 %then %do ;
                lb = pctlb[&time] ;
                ub = pctub[&time] ; 
           %end ;
           output ;
        end;
        keep  &time obs sim riskdiff lb ub  ;
        run;

    
  %end;
 

    goptions reset = all  /* display */ hsize=8 in vsize = 8 in device=png  gsfname=mygraphs  ;
    symbol1 line = 1 width = 3 interpol = line  color = red ;
    symbol2 line = 2 width = 3 interpol = line color = blue ;
    symbol3 line = 2 width = 3 interpol = line color = blue ;

    axis1 order = 0 to &nperiods by 1 minor = none label=none value=(h=2) ;
    axis2 minor = none value=(h=2) label = (h = 3  justify=center angle=90 "Cummulative incidence %upcase(&outcome) " ) ;
    axis3 minor = none value=(h=2) label= none  ;
	/**
	%if %bquote(&censor ) ^= %then %do;
		axis3 order = -0.05 to 0.05 by 0.01 minor = none value=(h=2) label= none  ;
	%end;
    ***/
   
    %if &outctype = binsurv   %then %do;
 
         proc gplot data = _both  ;
         plot   obrisk * &time ncrisk * &time / overlay   vaxis = axis2 haxis = axis1 name="outc" noframe nolegend;  
         plot  riskdiff * &time  %if &frombootstrap = 1 %then ub * &time lb * &time ; / overlay  haxis = axis1 vaxis= axis3 name = "outcdiff" noframe nolegend  ;   
         run;
         quit;
 
    %end;
     %if %bquote(&compevent)^= AND &compevent_cens = 0  %then %do;
        axis2 minor = none value=(h=2) label = (h = 3 justify=center angle=90 "Cummulative incidence %upcase(&compevent) " ) ;
        proc gplot data = _graph3   ;
        plot   obs * &time sim * &time / overlay    vaxis = axis2 haxis = axis1 name='crisk' noframe nolegend;
        plot  riskdiff * &time  %if &frombootstrap = 1 %then ub * &time lb * &time ; / overlay  haxis = axis1 vaxis=axis3 name = 'criskd' noframe nolegend ;
        run;
        quit;
    %end;

    


 
     data _allgraphs ;
     set &covmean ;
     run;

     %let nvar = %numargs(&simlist) ;
     %let ngraph = 1 ;
     %do ivar = 1 %to &nvar ;
          %let vartmp = %scan(&simlist,&ivar) ;

          data _onegraph ;
          set _allgraphs ;
          keep &time &vartmp s&vartmp &vartmp._diff &vartmp._lbp &vartmp._ubp ;
          run;

          %let vartitle = Mean &vartmp ;
          %if %bquote(&titledata)^= %then %do;
            data _null_ ;
            set &titledata ;
            call symput('vartitle',trim(left(&vartmp)));
            run; 
          %end;
      
         
          axis2 minor = none value=(h=2) label = ( h = 3 justify=center angle=90 "&vartitle" ) ;
         

         
         proc gplot data = _onegraph  ;
         plot   &vartmp * &time s&vartmp  * &time / overlay  vaxis = axis2 haxis = axis1 name="p%eval(&ngraph)" noframe nolegend;
         plot  &vartmp._diff * &time %if &frombootstrap = 1 %then &vartmp._ubp * &time &vartmp._lbp * &time ; / overlay  haxis = axis1  vaxis = axis3 name = "p%eval(&ngraph + 1)" vref=0 nolegend noframe ;
         run;
         quit;


         

         %let ngraph = %eval(&ngraph + 2) ;
      %end;

      * create a page of graphs, 4 per page = 2 covariates per page ;

      %let mycount = %eval(&nvar / 2) ;
     
      %if %eval(&mycount * 2) = &nvar %then %let myextra = 0 ;
      %else %let myextra = 1 ;

      

       %if %bquote(&compevent)^= AND &compevent_cens = 0  %then %let ngraph = 5 ;
       %else %let ngraph = 3 ;
       
%let ngraph = 1 ;       

       %if &sixgraphs = 1 %then %do;
            %let mycount = %eval(&nvar / 3) ;
            %let myextra = %eval(&nvar - &mycount * 3 );

           
       %end;

       ods pdf body="&gfilename"  ;

       %if &outctype = binsurv %then %do;
          %if %bquote(&compevent)^= AND &compevent_cens = 0  %then %let graphlist =  1:outc  2:crisk  3:outcdiff  4:criskd  5:gslide  ;
           %else %let graphlist = 1:outc 3:outcdiff  5:gslide  ;
       %end ;
       %else %let graphlist =  1:mean 3:mdiff 5:gslide ;

       proc greplay igout = GSEG nofs tc=work.tempcat ;
       template= newtemp2 ;
       treplay &graphlist ;
       run;
       quit;

       %do mypair = 1 %to &mycount ;
 
            %if &sixgraphs = 0 %then %do ;
                %let graph1 = p&ngraph ;
                %let graph3 = p%eval(&ngraph + 1) ;
                %let graph2 = p%eval(&ngraph + 2) ;
                %let graph4 = p%eval(&ngraph + 3) ;
                %let template = newtemp ;
                %let graphlist =  1:&graph1 2:&graph2 3:&graph3  4:&graph4 ;
                 
            %end;
            %else %do ;
                %let graph1 = p&ngraph ;
                %let graph2 = p%eval(&ngraph + 1) ;
                %let graph3 = p%eval(&ngraph + 2) ;
                %let graph4 = p%eval(&ngraph + 3) ;
                %let graph5 = p%eval(&ngraph + 4) ;
                %let graph6 = p%eval(&ngraph + 5) ;

                %let template = newtemp3 ;
                %let graphlist =  1:&graph1 2:&graph2 3:&graph3  4:&graph4 5:&graph5 6:&graph6;
                
           %end;
             
           

             proc greplay igout = GSEG nofs tc= work.tempcat /*sashelp.templt */
             template=  &template /* l2r2s */ ;
             treplay &graphlist ; /* 1:&graph1 2:&graph2 3:&graph3  4:&graph4 ;*/
             run;
             quit;

             %if &sixgraphs = 0 %then %let ngraph = %eval(&ngraph + 4) ;
             %else %let ngraph = %eval(&ngraph + 6) ;
      %end;
 

      %if &myextra > 0 %then %do;
           %let graph1 = p&ngraph ;
           %let graph2 = p%eval(&ngraph + 1) ;
           %if &sixgraphs = 0 %then %do ;
             
                %let template = newtemp ;
                %let graphlist =   1:&graph1   3:&graph2 ;
               
            %end;
            %else %do ;
                %let template = newtemp3 ;
                %let graphlist = 1:&graph1 2:&graph2 ;
                %if &myextra = 2 %then %do;
                       %let graph3 = p%eval(&ngraph + 2) ;
                       %let graph4 = p%eval(&ngraph + 3) ;
                       %let &graphlist = &graphlist 3:&graph3 4:&graph4 ;
                %end;
                
             %end;   
                 
            proc greplay igout = GSEG nofs tc= work.tempcat /* sashelp.templt */
             template= &template /*l2r2s */ ;
             treplay  &graphlist ;
             run;
             quit;
      %end;

       

  ods pdf close  ;

  proc datasets library = work nolist ;
  delete _allgraphs _anydouble _both _bsres _bsss _bsstd2 _cont _cont2 _dataos _datass 
     _estobsrisk _estobssurv  _obssurv_ _onegraph 
   %if %bquote(&compevent)^= AND &compevent_cens = 0  %then _graph3 _estncdeath _estncrsk _forgraph3 _ncss _sdgraph4 ;
   %if &outctype ^= binsurv %then _both10_ _both1_ _graph2 %if &frombootstrap = 1 %then _meanstd ;;
   ;
  quit;
  
 %mend ;


%macro createhazard ;
     %local firstint secondint ;
      data _calchazard_ ;
      calchazard = 1 ;
      _sample_ = &bsample ;
      run;

      %let firstint = %sysfunc(compress(%scan(&intcomp,1)));
      %let secondint = %sysfunc(compress(%scan(&intcomp,2)));

      data hazard1 ;
      set simulated&firstint;
      int = &firstint;
      keep int _newtime  _censor ;
      run;

      data hazard2 ;
      set simulated&secondint;
      int = &secondint ;
      keep int _newtime _censor ;
      run;

      data both ;
      set hazard1 hazard2 ;
      run;

      data _calchazard_ ;
      calchazard = 0 ;
      _sample_ = &bsample ;
      run;


      
data both ;
set both ;
rename  _censor = event _newtime = newtime ;

if int = &firstint then int = 0;
else int = 1 ;
run;

 
 *ods select none ;
 proc phreg data = both   ;
 ods output ParameterEstimates=_inthr0_  ; 
 model newtime*event(0) =  int / rl  %if %bquote(&compevent)^= AND &compevent_cens = 0  %then eventcode=1 ; ;
 %if %bquote(&compevent)^= AND &compevent_cens = 0  %then hazardratio 'Subdistribution Hazards' int ;;
 run;
 *ods select all ;

%if &bsample = 0 %then %do;
    proc sql ;
    select HazardRatio into :sample_hazard from _inthr0_ ;
    quit;
    run;

    %put  hazard = &sample_hazard ;
%end;
/* need to include the case where the program is run in chunks and only the hazard ratio
   is needed for sample = 0. In this case we need to store the hazard ratio in a data set
   and then load the data when the pieces are put back together. */
%if &bootstrap_hazard = 1 OR &chunked = 1 %then %do; 

    data _inthr0_ ;
    set _inthr0_ (keep = HazardRatio) ;
    _sample_ = &bsample ;
    run;

    %if &bsample = &sample_start %then %do;
         data &hazardname ;
         set _inthr0_;
         run;
    %end;
    %else %do;
         data &hazardname ;
         set &hazardname  _inthr0_ ;
         run;
    %end;

    proc datasets library = work nolist ;
    delete both hazard1 hazard2 _inthr0_ ;
    quit;
 %end;


 
%mend ;
 
%macro addvarcheck(vartype=1 /* 0 for outcome type variables, 1 for modeled covariates */, outvar=) ;
   /* want the check covXaddvars for modeled variables to add in the correct predictors , used in dataprep or listpred(?) 
     THE VAR &i COMES FROM THE CALLING MACRO 

 */
       %if &vartype = 0 %then %do;
            %let addvarlisttmp = &&&outvar.addvars ;
       %end;
       %else %if &vartype = 1 %then %do;
            %let addvarlisttmp = &&cov&i.addvars ;
            %let cov&i.addvars = ;           
       %end;


       %let tmplist = ;
       %let naddvars = %numargs(&addvarlisttmp);
       %do ii = 1 %to &naddvars ;
            %let vartmp = %scan(&addvarlisttmp,&ii) ;
            %let incovlist = 0 ;
            %if &vartmp in (&covlist) %then %let incovlist = 1;
            %put  roger ii = &ii , vartmp = &vartmp , incovlist = &incovlist ;
            %if &incovlist = 0 %then %let tmplist = &tmplist &vartmp ;           
            %else %if &incovlist = 1 %then %do;
               %let varpos = -1 ;
               %let poscheck = 1 ;
               %do %while(&varpos = -1);
                  %let poscheck = %eval(&poscheck + 1);
                  %let covtmp = %scan(&covlist,&poscheck);
                  %if &vartmp = &covtmp %then %do;
                       %let varpos = %eval(&poscheck - 1) ;
                  %end;
                  %if &poscheck = &ncov and &varpos = -1 %then %let varpos = -2 ;
               %end;
               
               %if &varpos > 0 %then %do;
                   %put cov = &&cov&i ;
                   %put covlist = &covlist ;
                   %let mtypeholder = &&cov&varpos.mtype ;
                   %let cov&varpos.mtype = all ;
                   %put  vartmp = &vartmp , varpos = &varpos   , ptype = &&cov&varpos.ptype , &mtypeholder , &&cov&varpos.mtype    ;
                
                   %if &vartype = 0 %then %let predtypelist =  %listpred(main,,&varpos,&varpos,_l2,_l1) ;
                   %if &vartype = 1 %then %let predtypelist =  %listpred(contemp,,&varpos,&i-1) %listpred(main,&i,&varpos,&varpos,_l3,_l2,_l1)  ;  
                %put predtypelist = &predtypelist ;
                %let cov&varpos.mtype = &mtypeholder ;
                %let tmplist = &tmplist &predtypelist ;
               
               %end;
            %end;
       %end;
        %if &vartype = 0 %then %let &&outvar.addvars = &tmplist ;
        %if &vartype = 1 %then  %let cov&i.addvars =  &tmplist ;
%mend ;