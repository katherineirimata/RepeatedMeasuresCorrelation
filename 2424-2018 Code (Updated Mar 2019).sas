/*********************************************************************************************/
/** Code: 2424-2018 SAS Code		 														**/
/** Author: Katherine Irimata																**/
/** Description: SAS Global Forum 2018 Paper 2424 macro code 								**/
/**		Evaluate correlation using Hamlett et al (2004) mixed model approach				**/
/**		Macros for inference using the following approaches:								**/
/** 	1) Normal approximation (CI and p-value)		       								**/
/** 	2) Nonparametric Cluster Bootstrap (CI)  											**/
/** 	3) Nonparametric Cluster Bootstrap (p-value)										**/
/** Notes: Macro code updated Mar 2019 														**/
/** 		a) Resolved reference &ID. in bootsamp											**/
/**			b) Renamed "results" in macros to results_BSCI and results_BSpvalue				**/
/**			c) Updated %MMCorr_BSCI with missing rho12 values								**/
/**			d) Set seed for random shuffling in %MMCorr_BSpvalue							**/
/*********************************************************************************************/


/*********************************************************************************************/
/** Appendix I: Normal approximation (CI and p-value)		       							**/
/*********************************************************************************************/
/**MACRO CODE**/
%macro MMCorr_NormalApprox(data=,ID=,rep=,var1=,var2=);
/*Convert to long format*/
data &data._long (drop=&var1. &var2. i); *univariate format;
	set &data. (keep=&ID. &rep. &var1. &var2.);

	array var[2] &var1. &var2.;
	do i = 1 to 2;
		Vtype = i;
 		Response = var[i];
		output;
	end;
run;

/*Fit Hamlett (2004) mixed model*/
proc mixed data=&data._long method=ml covtest asycov asycorr;
	class &ID. vtype &rep.;
	model response = vtype / solution ddfm=kr;
	random vtype / type=un subject=&ID. v vcorr;
	repeated vtype / type=un subject=&rep.(&ID.);
	ods output VCorr=VCorr ConvergenceStatus=CS asycorr=asycorr asycov=asycov CovParms=CovParms;
run;

/*Use delta method to find SE(rho) and calculate CI and p-value based on normality*/
proc iml;
	use Covparms;
	read all var {Estimate} into Cov_estimate;
	close Covparms;

	use Asycov;
	read all var {CovP1 CovP2 CovP3 CovP4 CovP5 CovP6} into asycov;
	close Asycov;

	a = Cov_estimate[1];
	b = Cov_estimate[2];
	c = Cov_estimate[3];
	g = Cov_estimate[4];
	h = Cov_estimate[5];
	i = Cov_estimate[6];

	rhohat = (b+h)/sqrt((a+g)*(c+i));

	df_da = -0.5*((b+h)*(c+i))/sqrt(((a+g)*(c+i))**3);
	df_db = 1/sqrt((a+g)*(c+i));
	df_dc = -0.5*(b+h)*(a+g)/sqrt(((a+g)*(c+i))**3);
	df_dg = -0.5*(b+h)*(c+i)/sqrt(((a+g)*(c+i))**3);
	df_dh = 1/sqrt((a+g)*(c+i));
	df_di = -0.5*(b+h)*(a+g)/sqrt(((a+g)*(c+i))**3);

	create partialderiv var {df_da df_db df_dc df_dg df_dh df_di};
	append;
	close partialderiv;

	use partialderiv;
	read all var {df_da df_db df_dc df_dg df_dh df_di} into partialderiv;
	close partialderiv;

	Sigma = asycov;

	rho_var = partialderiv*Sigma*t(partialderiv);
	rho_SE = sqrt(rho_var);

	l95=rhohat-1.96*rho_SE;
	u95=rhohat+1.96*rho_SE;
	p = (1-cdf("normal",abs(rhohat),0,rho_SE))*2;

	NormalApproxOutput = j(1,5,.);
	NormalApproxOutput[,1]=rhohat;
	NormalApproxOutput[,2]=rho_SE;
	NormalApproxOutput[,3]=l95;
	NormalApproxOutput[,4]=u95;
	NormalApproxOutput[,5]=p;
	Outtitle={'Estimate'  'Std Error'  '95% CI Lower Bound'  '95% CI Upper Bound' 'Pvalue'};
	Varnames=t("&var1. and &var2.");

	PRINT NormalApproxOutput[colname=Outtitle rowname=Varnames];

quit;
%mend MMCorr_NormalApprox;


/*********************************************************************************************/
/** Appendix II: Nonparametric Cluster Bootstrap (CI)										**/
/*********************************************************************************************/
/**MACRO CODE**/
%macro MMCorr_BSCI(data=,ID=,rep=,var1=,var2=,numSamples=);
/*Select bootstrap samples*/
proc sort data=&data.(keep=&ID.) out=&data._ID nodupkey;
	by &ID.;
run;
proc surveyselect data=&data._ID NOPRINT seed=864132 out=&data._boot(drop=NumberHits)
	method=urs samprate=1 OUTHITS reps=&NumSamples.;
run;
data &data._boot;
	set &data._boot;
	if first.&ID. then RptNum = 1;
	else RptNum + 1;
	by Replicate &ID.;
	ID2 = &ID. || "_" || strip(RptNum);
run;
data &data._boot;
	set &data._boot;
	if first.Replicate then SubjNum = 1;
	else SubjNum + 1;
	by Replicate;
run;

data results_BSCI;
	input Replicate rho12 convstatus reason $200. rho1 rho2;
run;

/*Bootstrap: Sample Individuals*/
%do j=1 %to &NumSamples.;
	proc sql;
		create table bootsamp
		as select a.*,b.&rep.,b.&var1.,b.&var2. 
		from &data._boot a
		left join &data. b			
		on &data._boot.&ID. = &data..&ID.
		where Replicate = &j.;
	quit;
	proc sort data=bootsamp;
		by &ID. RptNum &rep.;
	run;

	data &data._long (drop=&var1. &var2. i); *univariate format;
		set bootsamp (keep=&ID. RptNum ID2 &rep. &var1. &var2.);

		if first.ID2 then Visit2 = 1;
		else Visit2 + 1;
		by ID2;

		array var[2] &var1. &var2.;
		do i = 1 to 2;
			Vtype = i;
			Response = var[i];
			output;
			end;
	run;

	/*Fit Hamlett (2004) mixed model for bootstrap sample*/
	proc mixed data=&data._long method=ml;
		class ID2 vtype visit2;
		model response = vtype / solution ddfm=kr;
		random vtype / type=un subject=ID2 v vcorr;
		repeated vtype / type=un subject=visit2(ID2);
		ods output VCorr=VCorr ConvergenceStatus=CS;
	run;

	/*Calculate correlation*/
	proc iml;
		use Vcorr(keep=COL:);
		read all var _ALL_ into Vcorr;
		close Vcorr;

		use CS;
		read all var {Status} into convstatus;
		read all var {Reason} into reason;
		close CS;

		R = Vcorr[1:2,1:2];
		v = cusum( 1 || (ncol(R):2) ); 
		rho12 = remove(vech(R), v); 

		R2 = Vcorr[1:2,(1+2):(2+2)];
		rho1 = R2[1,1];
		rho2 = R2[2,2];

		Replicate = j(2*(2-1)/2,1,&j.);

		create corrcalc var {Replicate rho12 convstatus reason rho1 rho2};
		append;
		close corrcalc;
	quit;

	data results_BSCI;
		set results_BSCI corrcalc;
	run;

	dm "output;clear;log;clear;odsresults;select all;clear;";
%end;

data results_BSCI_2;
	set results_BSCI;
	where convstatus = 0;
run;

data &data._long (drop=&var1. &var2. i); *univariate format;
	set &data. (keep=&ID. &rep. &var1. &var2.);

	array var[2] &var1. &var2.;
	do i = 1 to 2;
		Vtype = i;
 		Response = var[i];
		output;
	end;
run;

/*Fit Hamlett (2004) mixed model for original dataset*/
proc mixed data=&data._long method=ml covtest asycov asycorr;
	class &ID. vtype &rep.;
	model response = vtype / solution ddfm=kr;
	random vtype / type=un subject=&ID. v vcorr;
	repeated vtype / type=un subject=&rep.(&ID.);
	ods output CovParms=CovParms;
run;

proc iml;
	use results_BSCI_2;
	read all var {rho12} into rho;
	close results_BSCI_2;

	numsamp = nrow(rho);
	use Covparms;
	read all var {Estimate} into Cov_estimate;
	close Covparms;

	a = Cov_estimate[1];
	b = Cov_estimate[2];
	c = Cov_estimate[3];
	g = Cov_estimate[4];
	h = Cov_estimate[5];
	i = Cov_estimate[6];

	rhohat = (b+h)/sqrt((a+g)*(c+i));

	BS_mean = mean(rho);

	z = quantile("Normal",mean(rho < rhohat));
	a = 0;

	qlb = cdf("Normal",z+((z+quantile("Normal",0.025))/(1-a*(z+quantile("Normal",0.025)))));
	qub = cdf("Normal",z+((z+quantile("Normal",0.975))/(1-a*(z+quantile("Normal",0.975)))));

	call qntl(lb,rho,qlb);
	call qntl(ub,rho,qub);

	BootstrapOutput = j(1,4,.);
	BootstrapOutput[,1]=rhohat;
	BootstrapOutput[,2]=BS_mean;
	BootstrapOutput[,3]=lb;
	BootstrapOutput[,4]=ub;
	Outtitle={'Estimate'  'Bootstrap Estimate'  '95% CI Lower Bound'  '95% CI Upper Bound'};
	Varnames=t("&var1. and &var2.");

	PRINT BootstrapOutput[colname=Outtitle rowname=Varnames];

quit;

%mend MMCorr_BSCI;



/*********************************************************************************************/
/** Appendix III: Nonparametric Cluster Bootstrap (p-value)  								**/
/*********************************************************************************************/
/**MACRO CODE**/
%macro MMCorr_BSpvalue(data=,ID=,rep=,var1=,var2=,numSamples=);
/*Select bootstrap samples*/
proc sort data=&data.(keep=&ID.) out=&data._ID nodupkey;
	by &ID.;
run;
proc surveyselect data=&data._ID NOPRINT seed=864132 out=&data._boot(drop=NumberHits)
	method=urs samprate=1 OUTHITS reps=&NumSamples.;
run;
data &data._boot;
	set &data._boot;
	if first.&ID. then RptNum = 1;
	else RptNum + 1;
	by Replicate &ID.;
	ID2 = &ID. || "_" || strip(RptNum);
run;
data &data._boot;
	set &data._boot;
	if first.Replicate then SubjNum = 1;
	else SubjNum + 1;
	by Replicate;
run;

data results_BSpvalue;
	input Replicate rho convstatus reason $200.;
run;

/*Bootstrap: Sample Individuals*/
%do j=1 %to &NumSamples.;
	proc sql;
		create table bootsamp
		as select a.*,b.&rep.,b.&var1.,b.&var2. 
		from &data._boot a
		left join &data. b			
		on &data._boot.&ID. = &data..&ID.
		where Replicate = &j.;
	quit;

	proc sort data=bootsamp;
		by &ID. RptNum &rep.;
	run;

	/*shuffle Wi's within a subject, rematch with Ui's to form m_i new pairs of observations*/ 
	data w_order;
		set bootsamp (keep=ID2 SubjNum &rep. &var2.);
		order = ranuni(&j.);
	run;

	proc sort data=w_order;
		by ID2 order;
	run;

	/*create line number variable and merge on the line number*/
	data bootsamp;
		set bootsamp(drop=&var2.);
		linenum = _n_;
	run;		
	data w_order;
		set w_order(keep=&var2.);
		linenum = _n_;
	run;

	data bootsamp;
		merge bootsamp w_order;
		by linenum;
	run;

	data &data._long (drop=&var1. &var2. i); *univariate format;
		set bootsamp (keep=&ID. RptNum ID2 &rep. &var1. &var2.);

		if first.ID2 then Visit2 = 1;
		else Visit2 + 1;
		by ID2;

		array var[2] &var1. &var2.;
		do i = 1 to 2;
 			Vtype = i;
			Response = var[i];
			output;
		end;
	run;

	/*Fit Hamlett (2004) mixed model for bootstrap sample*/
	proc mixed data=&data._long method=ml;
		class ID2 vtype visit2;
		model response = vtype / solution ddfm=kr;
		random vtype / type=un subject=ID2 v vcorr;
		repeated vtype / type=un subject=visit2(ID2);
		ods output VCorr=VCorr ConvergenceStatus=CS V=V SolutionF=SolutionF;
	run;

	/*Calculate correlation for each bootstrap sample*/
	proc iml;
		use Vcorr(keep=COL:);
		read all var _ALL_ into Vcorr;
		close Vcorr;

		use V(keep=COL:);
		read all var _ALL_ into VCov;
 		close V;

		use CS;
		read all var {Status} into convstatus;
		read all var {Reason} into reason;
		close CS;

		use SolutionF;
		read all var {Estimate} into mu;
		close SolutionF;

		R = Vcorr[1:2,1:2];
		v = cusum( 1 || (ncol(R):2) );
		rho = remove(vech(R), v);

        Replicate = j(2*(2-1)/2,1,&j.);

		create corrcalc var {Replicate rho convstatus reason};
		append;
		close corrcalc;
	quit;

	data results_BSpvalue;
		set results_BSpvalue corrcalc;
	run;

	dm "output;clear;log;clear;odsresults;select all;clear;";
%end;

data results_BSpvalue_2;
	set results_BSpvalue;
	where convstatus = 0;
run;

data &data._long (drop=&var1. &var2. i); *univariate format;
	set &data. (keep=&ID. &rep. &var1. &var2.);

	array var[2] &var1. &var2.;
	do i = 1 to 2;
		Vtype = i;
 		Response = var[i];
		output;
	end;
run;

/*Fit Hamlett (2004) mixed model for original dataset*/
proc mixed data=&data._long method=ml covtest asycov asycorr;
	class &ID. vtype &rep.;
	model response = vtype / solution ddfm=kr;
	random vtype / type=un subject=&ID. v vcorr;
	repeated vtype / type=un subject=&rep.(&ID.);
	ods output CovParms=CovParms;
run;

/*Calculate p-value*/
proc iml;
	use results_BSpvalue_2;
	read all var {rho} into rho;
	close results_BSpvalue_2;

	BS_mean = mean(rho);

	numsamp = nrow(rho);

	use Covparms;
	read all var {Estimate} into Cov_estimate;
	close Covparms;

	a = Cov_estimate[1];
	b = Cov_estimate[2];
	c = Cov_estimate[3];
	g = Cov_estimate[4];
	h = Cov_estimate[5];
	i = Cov_estimate[6];

	rhohat = (b+h)/sqrt((a+g)*(c+i));

	if rhohat >0 then p = (sum(rho < 0)+sum(rho > 2*rhohat))/numsamp;
	else p = (sum(rho > 0)+sum(rho < 2*rhohat))/numsamp;

	BootstrapOutput = j(1,3,.);
	BootstrapOutput[,1]=rhohat;
	BootstrapOutput[,2]=BS_mean;
	BootstrapOutput[,3]=p;
	Outtitle={'Estimate'  'Bootstrap Estimate'  'Pvalue'};
	Varnames=t("&var1. and &var2.");

	PRINT BootstrapOutput[colname=Outtitle rowname=Varnames];
quit;

%mend MMCorr_BSpvalue;
