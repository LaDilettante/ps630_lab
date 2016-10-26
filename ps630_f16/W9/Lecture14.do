/*Set-Up */
#delimit;
pause on;
clear all;
pause;
set mem 15m;
pause;
set more off;
pause;
cd "C:\data\QM2\Winter2011";
pause;
use saving.DTA;
save saving.DTA, replace;
pause;
save lec20.smcl, replace;
pause;


/*Wooldridge's notes on the data set: "Notes:  I remember entering this data set in the late 1980s, 
and I am pretty sure it came directly from an introductory econometrics text.  
But so far my search has been fruitless.  If anyone runs across this data set, 
I would appreciate knowing about it.*/
pause;

describe;
pause;

/*Wooldridge p. 287; Example 8.6*/
pause;

reg sav inc;
pause;

predict uhat1, r;
pause;

scatter uhat1 sav, xtitle("Dependent Variable=Saving") ytitle("Residuals");
pause;

/*This graph does not show the classic case of fanning-out along different
levels of the dependent varaible*/
pause;

rvfplot, xtitle("Fitted Values of Savings") ytitle("Residuals from Bivariate Regression");
pause;

/*The RVF plot does cause me to worry somewhat. We see very little variance along
levels of the fitted values until we hit about $1000 and then the variance increases, 
declines, and increases again.*/
pause;

/*Let's move to a better test of the heteroskedasticity, in order to diagnose the problem.*/
pause;

/*White Test*/
/*White, H. "A Heteroskedasticity-Consistent Covariance Matrix Estimator 
      and a Direct Test for Heteroskedasticity." Econometrica, 48, 1980,
      817-838.*/

pause;

/*W1. Estimate the original regression. Done and done!*/
pause;
/*W2. Square the residuals of the regession*/
pause;
generate uhat1_sq=uhat1^2;

/*W3. Generate the predicted values from the regression*/
pause;
predict yhat1;
pause;

/*W4. Square the predicted values.*/
generate yhat1_sq=yhat1^2;
pause;

/*W5. Estimate a regression of uhat1_sq on predicted values and the squared predicted values.*/
pause;
reg  uhat1_sq  yhat1 yhat1_sq;

/*W6. Now, we have two options. Calculate an F-distribution or a Chi2-distribution*/
pause;

/*W6a. F-test (R2/K)/[(1-R2)/(n-k-1)]*/
pause;

display (0.0185/2)/[(1-.0185)/(100-2-1)];
pause;

/*W6b. Calculate the p-value. I am too lazy to look in the back of the book.*/
display Ftail(2,97,.914);
pause;

/*W7. Using the F-test, we cannot reject the null-hypothesis of constant variance*/
pause;

/*How about the chi2 distribution*/
pause;

/*W7a. Calculate the Lagrange Multiplier = n*R-sq*/
pause;
display 100*(.0185);
pause;


/*W7b. Calculate the p-value from the chi-sq distribution using the LM and k+1*/
pause;
display chi2tail(2,1.85);
pause;

/*We cannot reject the null of constant variance*/

/*W8. Confirm using the white test*/
pause;
reg sav inc;
pause;
whitetst;
pause;

/*Breusch-Pagan test - You know this as hettest.*/
/*Breusch, T. and A. Pagan. "A Simple Test for Heteroskedasticity and Random
Coefficient Variation." Econometrica, 47, 1979, 1287-1294.*/
pause;

/*BPAGAN according to Chris Baum of Boston College*/
pause;

/*BP Baum 4. Summarize the squared residuals/Wooldridge 180*/
pause;
sum uhat1_sq;

/*BP Baum 5. Divide the squared residuals by the average squared residual. r() saves the result of the sum command*/
gen ravg=1/r(mean)*uhat1_sq;
pause;

/*BP Baum 6. Regress ravg on income*/
pause;
reg ravg inc;
pause;

/*BP Baum 7. Save a scaler called df which is equal to the degrees of freedom from the
above model in BP Baum 6*/
pause;
sca df=e(df_m);
pause;
sca list df;
pause;

/*BP Baum 8. Save a scalar called bpagan which is equal to half of the of the 
explained sum of squares from the model in step 6*/
pause;
sca bpagan=0.5*e(mss);
pause;


/*BP Baum 9. List the scaler Bagan*/
sca list bpagan;
pause;

/*BP Baum 10. Create a scalar called p based on the model degrees of freedom and the bpagan scalar*/ 
sca p=chiprob(df,bpagan);
pause;

/*BP Baum 11. List the p scaler, which is the p-statistic of the chi_sq distribution*/
sca list p;
pause;

/*Confirm*/
reg sav inc;
pause;
hettest;
pause;
/*The nice thing about the bpgan command is that it is possible to identify a particular variable.*/
bpagan inc;
pause;

/*The tests are contradictory, so let's move straight to robust*/ 
pause;


/*According to this test, we must reject the null-hypothesis of homoskedasticity*/

/*ROBUST STANDARD ERRORS - How are they caculated?*/
pause;
reg sav inc size educ age black;
pause;


/*Robust 1 - Save residuals*/
predict uhat_robust, resid;
pause;

/*Robust 2 - Square residuals*/
generate uhat_rob_sq= uhat_robust^2;
pause;

/*Robust 3 - To get the robust standard error for education, I regress education on all the other independent variables
in the model*/
pause;

reg educ inc size  age black;
pause;

/*Robust 4 - Save the residuals from that regression*/
pause;
predict rhat_educ, resid;
pause;

/*Robust 5 - Sqare rhat_income*/

generate rhat_edu_sq= rhat_educ^2;
pause;

/*Robust 6 - Generate the numerator of the robust equation*/
pause; 

generate numerator_rob= uhat_rob_sq*rhat_edu_sq;
pause;

/*Robust 7 - Sum up the numerator*/
pause;

 tabstat  numerator_rob, statistics( sum );
pause; 

/*Robust 8 - Calculate the formula on page 274 of Wooldridge*/
pause;

display  1.57e+10/(758.249262^2);
pause;

/*Robust 9 - Add a correction for Degrees of Freedom (Mackinnon and White, 1985). Why?  If the 
squared uhats for the same for all observations i (the strongerst homoskedasticity), we would get the usual OLS Standard Errors.*/
pause;

display 27307.105*(100/(100-5-1));
pause;


/*Robust 10 - Take the square root of the Robust Variance around Beta Hat to get the standard error*/
pause;

display sqrt(29050.112);
pause;


/*Compare to Robust*/
pause;
reg sav inc size educ age black, robust;
pause;
/*How do they compare?*/
pause;

/*The difference is Robust makes the correction for all of the variables.  Now, 
some of you asked why does Robust correct some and not others.  The answer, as you can see now, 
is the size of the rhat for each variable.  Varibles entirely explained by the other independents, 
will require only very small standard error corrections.*/
pause;

/*Notice that the robust correction alters the standard errors, but has no effect on the coefficients*/
pause;

/*Weighted Least Squares*/
pause;

reg sav inc;
pause;

hettest;
pause;


/*With the Breusch-Pagan test, we reject the null, so I decide to reweight regression by the variance
in income. Basically, I divide the variance of income out of the error-covariance matrix, 
so that only sigma squared remain on the diagnols.*/
pause;

generate h=sqrt(inc);
pause;
/*WLS2: Re-Weight DV, IVs, constant*/
pause;
generate hsav=sav/h;
pause;
generate hinc=inc/h;
pause;
generate weight=1/h;
pause;
/*Run weighted regression without the constant*/
pause;
reg hsav hinc  weight, nocon;
pause;
/*Using stata aweight subcommand
aweights, or analytic weights, are weights that are inversely
proportional to the variance of an observation; i.e., the variance of
 the j-th observation is assumed to be sigma^2/w_j, where w_j are the
 weights.  Typically, the observations represent averages and the
weights are the number of elements that gave rise to the average.
For most Stata commands, the recorded scale of aweights is
irrelevant; Stata internally rescales them to sum to N, the number of
observations in your data, when it uses them.*/
pause;
reg sav inc [aweight=1/inc];
pause;

/*Repeat using the full regression. Once again we assume income is 
the source of the variance*/
reg sav inc size educ age black;
pause;
whitetst;
pause;
hettest;
pause;

/*Generate weighte controls*/
pause;
generate hsize=size/h;
generate heduc=educ/h;
generate hage=age/h;
generate hblack=black/h;

/*Run WLS regression*/
pause;

reg  hsav hinc hsize heduc hage hblack weight, nocon;
pause;

/*Compare to STATA*/
pause;
reg sav inc size educ age black [aweight=1/inc];
pause;


/*What if the exact form of heteroskedasticity is not obvious*/
pause;

/*Feasible Generalizable Least Squares removes the variance from error-covariance
matrix, by estimating the contribution of the independent variables to the observed 
variance (the squared residuals) and re-weighting the model by the predicted values from 
this regression.*/
pause;

/*FGLS1: Run the regression of y on the x's and obtain the residual*/
reg sav inc size educ age black;
predict uhat_FGLS, resid;
pause;

/*FGLS2:Take the natural log of the squared residuals*/
generate luhat2_FGLS=ln(uhat_FGLS^2);
pause;

/*FGLS3: Regress the log of the squared residulas on the independent 
variables*/
pause;
reg  luhat2_FGLS inc size educ age black;
pause;

/*FGLS4:Obtain the fitted values*/
predict g_hat;
pause;

/*FGLS5: Exponentiate the fitted values*/
pause;
generate h_hat=exp(g_hat);
pause;

/*FGLS 6: Take the square root of h_hat*/
generate hhat_sqrt=sqrt(h_hat);
pause;

/*FGLS7: Estimate the new equation weighting by hhat_sqrt rather than income*/
pause;

generate savehat=sav/hhat_sqrt;
generate inchat=inc/hhat_sqrt;
generate sizhat=size/hhat_sqrt;
generate agehat=age/hhat_sqrt;
generate blackhat=black/hhat_sqrt;
generate educhat=educ/hhat_sqrt;
generate weight_FLGS=1/hhat_sqrt;
pause;

/*FGLS7:Run the weighted regression*/
pause;

reg savehat inchat sizhat educhat agehat blackhat weight_FLGS, nocons;
pause;

/*Compare the result to the weighting procedure*/
pause;

reg sav inc size educ age black [aweight=1/h_hat];
pause;

/*Notice that WLS and FGLS both bias the coefficients in the course
of making the standard error corrections.  We are trading unbiasedness for efficiency*/

reg sav inc size educ age black;
pause;






 







