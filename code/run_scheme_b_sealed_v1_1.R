args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L) stop("usage: run_scheme_b_sealed_v1_1.R <data-dir> <sealed-dir>")
Sys.setlocale("LC_ALL","C"); Sys.setenv(TZ="Asia/Shanghai")
options(survey.lonely.psu="adjust",survey.adjust.domain.lonely=TRUE)
suppressPackageStartupMessages(library(survey)); suppressPackageStartupMessages(library(jsonlite))
set.seed(20260911); input_dir<-normalizePath(args[[1]],winslash="/",mustWork=TRUE); out_dir<-args[[2]]; dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
ck_levels<-c("00","10","01","11")

read_data<-function(country){
 d<-read.csv(file.path(input_dir,paste0("restricted_scheme_b_",country,".csv")),stringsAsFactors=FALSE,check.names=FALSE)
 d$CK_B<-sprintf("%02d",as.integer(d$CK_B)); d$CK<-factor(d$CK_B,levels=ck_levels)
 d$C<-factor(d$C_B,levels=c(0,1)); d$K<-factor(d$K_B,levels=c(0,1))
 d$age_group<-factor(d$age_group,levels=c("15-19","20-24","25-29","30-34","35-39","40-44","45-49")); d$pord_group<-factor(d$pord_group,levels=c("1","2-4","5+"))
 d$outcome_year<-factor(d$outcome_year); d$calendar_month<-factor(sprintf("%02d",as.integer(d$calendar_month)),levels=sprintf("%02d",1:12)); d$v024<-factor(d$v024); d$v025<-factor(d$v025); d$v106<-factor(d$v106); d
}
make_design<-function(d) svydesign(ids=~v021,strata=~v022,weights=~weight,data=d,nest=TRUE)
required<-c("Y","CK","weight","v021","v022","v024","v025","outcome_year","calendar_month","age_group","pord_group","v106")

support_gate<-function(d){
 use<-!is.na(d$Y)
 cells<-lapply(ck_levels,function(g){z<-d[use & d$CK_B==g,,drop=FALSE]; yy<-table(factor(z$Y,levels=c(0,1))); list(ck=g,outcome_known_n=nrow(z),independent_clusters=length(unique(z$v021)),smaller_outcome_category_n=min(as.integer(yy)),pass=nrow(z)>=100L&&length(unique(z$v021))>=50L&&min(as.integer(yy))>=30L)})
 names(cells)<-ck_levels; list(pass=all(vapply(cells,function(x)x$pass,logical(1))),cells=cells)
}
overlap_gate<-function(d,use){
 dims<-c("v024","outcome_year","calendar_month"); total<-sum(d$weight[use])
 items<-lapply(dims,function(v){lv<-unique(d[[v]][use]); ok<-lv[vapply(lv,function(a)all(ck_levels%in%unique(d$CK_B[use&d[[v]]==a])),logical(1))]; coverage<-sum(d$weight[use&d[[v]]%in%ok])/total; list(levels=length(lv),levels_with_all_four=length(ok),weighted_coverage=unname(coverage),pass=coverage>=.8)})
 names(items)<-dims; list(threshold=.8,pass=all(vapply(items,function(x)x$pass,logical(1))),dimensions=items)
}
binary_overlap<-function(d,use){
 dims<-c("v024","outcome_year","calendar_month"); out<-list(threshold=.8,C=list(),K=list())
 for(ex in c("C_B","K_B")) for(v in dims){lv<-unique(d[[v]][use]); ok<-lv[vapply(lv,function(a)all(c(0,1)%in%unique(d[[ex]][use&d[[v]]==a])),logical(1))]; cv<-sum(d$weight[use&d[[v]]%in%ok])/sum(d$weight[use]); out[[if(ex=="C_B")"C" else "K"]][[v]]<-list(weighted_coverage=unname(cv),pass=cv>=.8)}
 out$pass<-all(c(vapply(out$C,function(x)x$pass,logical(1)),vapply(out$K,function(x)x$pass,logical(1)))); out
}

num_list<-function(x){y<-as.list(as.numeric(x));names(y)<-names(x);y}
qsummary<-function(x){z<-quantile(as.numeric(x),c(0,.01,.25,.5,.75,.99,1),na.rm=TRUE,names=FALSE);names(z)<-c("min","p01","p25","p50","p75","p99","max");num_list(z)}

standardize<-function(fit,target,variable,levels_value,link="logit"){
 b<-coef(fit); vv<-vcov(fit); w<-target$weight; risks<-numeric(length(levels_value)); J<-matrix(NA_real_,length(levels_value),length(b),dimnames=list(levels_value,names(b)))
 for(i in seq_along(levels_value)){nd<-target; if(is.factor(nd[[variable]]))nd[[variable]]<-factor(levels_value[[i]],levels=levels(nd[[variable]])) else nd[[variable]]<-levels_value[[i]]; mm<-model.matrix(delete.response(terms(fit)),nd,contrasts.arg=fit$contrasts,xlev=fit$xlevels); mm<-mm[,names(b),drop=FALSE]; eta<-drop(mm%*%b); p<-if(link=="logit")plogis(eta) else eta; risks[[i]]<-sum(w*p)/sum(w); slope<-if(link=="logit")p*(1-p) else rep(1,length(p)); J[i,]<-colSums(mm*(w*slope))/sum(w)}
 names(risks)<-levels_value; list(risks=risks,J=J,cov=J%*%vv%*%t(J))
}
linear_est<-function(std,w,vv,df,level=.95){est<-sum(w*std$risks[names(w)]);g<-drop(w%*%std$J[names(w),,drop=FALSE]);se<-sqrt(drop(t(g)%*%vv%*%g));crit<-qt((1+level)/2,df);list(estimate=unname(est),standard_error=unname(se),conf_low=unname(est-crit*se),conf_high=unname(est+crit*se),p_value_unadjusted_supportive=unname(2*pt(-abs(est/se),df)))}

run_ng<-function(d,gate,overlap){
 use<-complete.cases(d[,required]); form<-Y~CK+outcome_year+calendar_month+v024+v025+age_group+pord_group+v106; mm<-model.matrix(form,d[use,]); rank0<-qr(mm)$rank
 pre<-list(support=gate,overlap=overlap,design_matrix=list(rows=sum(use),columns=ncol(mm),rank=rank0,full_rank=rank0==ncol(mm)))
 if(!gate$pass||!overlap$pass||rank0!=ncol(mm)) return(list(status="STOP_ITEM_SUPPORT_OR_OVERLAP",pre_model=pre,model_fitted=FALSE))
 des<-subset(make_design(d),use); target<-d[use,]; warns<-character(); fit<-withCallingHandlers(tryCatch(svyglm(form,design=des,family=quasibinomial(),influence=TRUE),error=function(e)e),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")})
 log_ok<-!inherits(fit,"error")&&isTRUE(fit$converged)&&fit$rank==length(coef(fit))&&all(is.finite(coef(fit)))&&all(is.finite(vcov(fit)))&&!any(grepl("0 or 1|did not converge|singular",warns,ignore.case=TRUE)); model_type<-"survey_weighted_logistic"
 if(!log_ok){warns<-c(warns,if(inherits(fit,"error"))conditionMessage(fit) else "logistic_failure_rule_triggered");fit<-svyglm(form,design=des,family=gaussian(),influence=TRUE);if(!isTRUE(fit$converged)||fit$rank!=length(coef(fit))||degf(des)<=0)stop("Scheme-B LPM fallback failed");model_type<-"survey_weighted_lpm_fallback"}
 std<-standardize(fit,target,"CK",ck_levels,if(model_type=="survey_weighted_logistic")"logit" else "identity"); vv<-vcov(fit);df<-degf(des);crit<-qt(.975,df)
 risks<-lapply(ck_levels,function(g){se<-sqrt(std$cov[g,g]);list(estimate=unname(std$risks[[g]]),standard_error=unname(se),conf_low=unname(std$risks[[g]]-crit*se),conf_high=unname(std$risks[[g]]+crit*se))});names(risks)<-ck_levels
 cw<-list(RD_10_00=c("00"=-1,"10"=1,"01"=0,"11"=0),RD_01_00=c("00"=-1,"10"=0,"01"=1,"11"=0),RD_11_00=c("00"=-1,"10"=0,"01"=0,"11"=1),IC_add=c("00"=1,"10"=-1,"01"=-1,"11"=1)); contrasts<-lapply(cw,function(w)linear_est(std,w,vv,df))
 infl<-attr(fit,"influence"); pinf<-if(is.null(infl))NULL else qsummary(sqrt(rowSums(rowsum(infl,target$v021,reorder=FALSE)^2)))
 co<-data.frame(tag="nigeria_scheme_B_v1_1",term=names(coef(fit)),estimate=unname(coef(fit)),standard_error=sqrt(diag(vv)),stringsAsFactors=FALSE); cv<-as.data.frame(as.table(vv),stringsAsFactors=FALSE);names(cv)<-c("row_term","column_term","covariance");cv$tag<-"nigeria_scheme_B_v1_1";cv<-cv[,c("tag","row_term","column_term","covariance")]
 esttab<-data.frame(tag="nigeria_scheme_B_v1_1",estimand=c(paste0("R_",ck_levels),names(contrasts)),estimate=c(std$risks,vapply(contrasts,function(x)x$estimate,numeric(1))),stringsAsFactors=FALSE)
 list(status="PASS",pre_model=pre,model_fitted=TRUE,model_type=model_type,n=nrow(target),risks=risks,contrasts=contrasts,diagnostics=list(converged=isTRUE(fit$converged),columns=ncol(model.matrix(fit)),rank=fit$rank,covariance_rank=qr(vv)$rank,coefficient_count=length(coef(fit)),design_df=df,warnings=unique(warns),prediction_range_pass=min(std$risks)>=0&&max(std$risks)<=1,residual_quantiles=qsummary(residuals(fit,type="pearson")),anonymous_psu_influence_quantiles=pinf),coefficients=co,covariance=cv,estimates=esttab)
}

run_bf_direct<-function(d){use<-!is.na(d$Y);des<-subset(make_design(d),use);fit<-svyglm(Y~0+CK,design=des,family=quasibinomial());b<-coef(fit);vv<-vcov(fit);r<-plogis(b);der<-r*(1-r);rc<-diag(der)%*%vv%*%diag(der);names(r)<-sub("^CK","",names(b));rownames(rc)<-colnames(rc)<-names(r);df<-degf(des);crit<-qt(.975,df);cells<-lapply(ck_levels,function(g){se<-sqrt(rc[g,g]);list(estimate=unname(r[[g]]),standard_error=unname(se),conf_low=unname(r[[g]]-crit*se),conf_high=unname(r[[g]]+crit*se))});names(cells)<-ck_levels;list(method="direct_survey_weighted_cells",p_values_reported=FALSE,cells=cells,design_df=df,rank=fit$rank,columns=ncol(model.matrix(fit)))}

run_bf_reduced<-function(d){req<-c("Y","C","K","weight","v021","v022","v024","v025","outcome_year","calendar_month","age_group","pord_group","v106");use<-complete.cases(d[,req]);form<-Y~C+K+outcome_year+calendar_month+v024+v025+age_group+pord_group+v106;mm<-model.matrix(form,d[use,]);ov<-binary_overlap(d,use);pre<-list(overlap=ov,design_matrix=list(rows=sum(use),columns=ncol(mm),rank=qr(mm)$rank,full_rank=qr(mm)$rank==ncol(mm)));if(!ov$pass||qr(mm)$rank!=ncol(mm))return(list(status="STOP_ITEM_REDUCED_SUPPORT",pre_model=pre,model_fitted=FALSE));des<-subset(make_design(d),use);fit<-svyglm(form,design=des,family=quasibinomial());vv<-vcov(fit);co<-data.frame(tag="burkina_scheme_B_reduced_C_plus_K",term=names(coef(fit)),estimate=unname(coef(fit)),standard_error=sqrt(diag(vv)),stringsAsFactors=FALSE);cv<-as.data.frame(as.table(vv),stringsAsFactors=FALSE);names(cv)<-c("row_term","column_term","covariance");cv$tag<-"burkina_scheme_B_reduced_C_plus_K";cv<-cv[,c("tag","row_term","column_term","covariance")];list(status="PASS",pre_model=pre,model_fitted=TRUE,diagnostics=list(converged=isTRUE(fit$converged),rank=fit$rank,columns=ncol(model.matrix(fit)),covariance_rank=qr(vv)$rank,design_df=degf(des)),coefficients=co,covariance=cv)}

ng<-read_data("nigeria");bf<-read_data("burkina_faso"); ng_use<-complete.cases(ng[,required]);bf_use<-complete.cases(bf[,required]);ng_gate<-support_gate(ng);bf_gate<-support_gate(bf);ng_overlap<-overlap_gate(ng,ng_use);bf_overlap<-overlap_gate(bf,bf_use)
ng_result<-run_ng(ng,ng_gate,ng_overlap);bf_direct<-run_bf_direct(bf);bf_reduced<-run_bf_reduced(bf)
coefs<-list();covs<-list();ests<-list();if(isTRUE(ng_result$model_fitted)){coefs[[1]]<-ng_result$coefficients;covs[[1]]<-ng_result$covariance;ests[[1]]<-ng_result$estimates;ng_result$coefficients<-NULL;ng_result$covariance<-NULL;ng_result$estimates<-NULL};if(isTRUE(bf_reduced$model_fitted)){coefs[[2]]<-bf_reduced$coefficients;covs[[2]]<-bf_reduced$covariance;bf_reduced$coefficients<-NULL;bf_reduced$covariance<-NULL}
manifest<-fromJSON(file.path(input_dir,"scheme_b_derivation_manifest.json"),simplifyVector=FALSE)
payload<-list(state="SCHEME_B_MODEL_OUTPUT_SEALED",amendment="v1.1",sample_flow=list(nigeria=list(eligible=nrow(ng),known=sum(!is.na(ng$Y)),complete=sum(ng_use)),burkina_faso=list(eligible=nrow(bf),known=sum(!is.na(bf$Y)),complete=sum(bf_use))),support_gates=list(nigeria=ng_gate,burkina_faso=bf_gate),overlap_gates=list(nigeria=ng_overlap,burkina_faso=bf_overlap),nigeria=ng_result,burkina_faso=list(adjusted_CxK_fitted=FALSE,direct=bf_direct,reduced_C_plus_K=bf_reduced),derivation_manifest=manifest,displacement_200="WITHDRAWN_WITHOUT_REPLACEMENT",survey_options=list(lonely_psu=getOption("survey.lonely.psu"),domain_lonely=isTRUE(getOption("survey.adjust.domain.lonely"))))
write_json(payload,file.path(out_dir,"scheme_b_sealed_results.json"),pretty=TRUE,auto_unbox=TRUE,digits=16,na="null")
if(length(coefs))write.csv(do.call(rbind,coefs),file.path(out_dir,"scheme_b_coefficients.csv"),row.names=FALSE,quote=TRUE,na="")
if(length(covs))write.csv(do.call(rbind,covs),file.path(out_dir,"scheme_b_covariances.csv"),row.names=FALSE,quote=TRUE,na="")
if(length(ests))write.csv(do.call(rbind,ests),file.path(out_dir,"scheme_b_estimates.csv"),row.names=FALSE,quote=TRUE,na="")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"R_sessionInfo.txt"),useBytes=TRUE)
