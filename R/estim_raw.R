#'@title Estimation of Four-Fold (2x2) Table Cell Frequencies (Raw Data) from Effect Size Measures
#'
#'@description Estimation of four-fold (2x2) table cell frequencies (raw data) from risk ratios (relative risks), risk differences and odds ratios.
#'
#'@param es Value of the effect size (summary statistic).
#'
#'@param lb Lower bound of the 95% confidence interval of the effect size.
#'
#'@param ub Upper bound of the 95% confidence interval of the effect size.
#'
#'@param m1 Total number of participants in the exposed/treated group.
#'
#'@param m2 Total number of participants in the unexposed/control group.
#'
#'@param e1 Total number of participants developing the outcome event (optional).
#'
#'@param dec Number of decimal places with which effect sizes are being presented.
#'
#'@param measure Character string indicating the type of effect size measure. Possible options are the risk ratio/relative risk ("rr"), the risk difference ("rd"), and the odds ratio ("or").
#'
#'@section Details and context
#'While raw data can be useful for doing meta-analysis, such data is often not provided by primary studies (with summary statistics being solely presented). Therefore, based on summary statistics (namely, risk ratios, risk differences and odds ratios), this function estimates the value of each cell in a 2x2 table according to the equations described by
#'C. Di Pietrantonj (Statistics in Medicine 2006;25:2299-2322). On a technical note, please note that raw data are estimated based on equations that give priority to computation of events number in the treatment group.
#'The following notation is used (see also table below):
#'\enumerate{
#'    \item **a:** Patients who are exposed/treated and develop the outcome event.
#'    \item **b:** Patients who are exposed/treated and do not develop the outcome event.
#'    \item **c:** Patients who are exposed/treated and do not develop the outcome event.
#'    \item **d:** Patients who are not exposed/treated and do not develop the outcome event.
#'    }
#'
#'|           | Event | No event | Total |
#'| --------- | ----- | -------- | ----- |
#'|  Exposed  |   a   |    b     |  m1   |
#'| --------- | ----- | -------- | ----- |
#'| Unexposed |   c   |    d     |  m2   |
#'| --------- | ----- | -------- | ----- |
#'|   Total   |  e1   |          |       |
#'
#'Estimating raw data from risk differences and odds ratio involves application of the quadratic formula. Therefore, for such effect sizes, two results sets ("solutions") may be presented, particularly, when there is no information on the total number of participants developing the outcome event (e1). In those cases, the researcher must choose the most adequate solution based on prior clinical knowledge.
#'
#'The precision with which effect sizes are provided conditions the uncertainty underlying raw data estimation. For example, effect sizes ranging from 1.45 to 1.54 can be reported as '1.5'. Therefore, when e1 is unknown, this function provides not only point estimates based on exact input values, but also a range of compatible values for a, b, c and d (considering possible roundings of effect sizes). When e1 is known, this function only provides the result(s) in which a+c=e1.
#'
#'@section Practical examples
#'
#'Herein, are discussed four practical examples, corresponding to four different scenarios
#'
#'## **Low effect size precision and no information on the number of participants developing the outcome of event:**
#'
#'Consider a primary study presenting the Risk Ratio as 0.6, with a 95% confidence interval of 0.4-0.9 (worked example 3.1. of C. Di Pietrantonj (Statistics in Medicine 2006;25:2299-2322)).
#'Let us assume that (i) We know that there are 352 participants in the treated/exposed group, and 376 participants in the control group;
#' (ii) we do not know the total number of participants developing the outcome of event.
#'Raw data can be estimated by:
#'
#'estim_rr <- estim_raw(es=0.6,lb=0.4,ub=0.9,m1=352,m2=376,dec=1,measure="rr")
#'estim_rr
#'
#'The results indicate that the value of a (corresponding to the number of treated patients who develop the outcome event) lies between
#'22 and 47, and that the value of c (corresponding to the number of control patients who develop the outcome event) lies between 38 and
#'87. This large range reflects the low estimate precision and lack of information on the number of participants developing the outcome event.
#'The real values are actually 37 and 65.
#'
#'
#'## **High effect size precision and no information on the number of participants developing the outcome of event:**
#'
#'Consider a primary study presenting the Odds Ratio as 0.6207, with a 95% confidence interval of 0.3382-1.1391 (worked example 4.2. of C. Di Pietrantonj
#'(Statistics in Medicine 2006;25:2299-2322)).
#'Let us assume that (i) We know that there are 355 participants in the treated/exposed group, and 366 participants in the control group;
#'(ii) we do not know the total number of participants developing the outcome of event.
#'Raw data can be estimated by:
#'
#'estim_or1 <- estim_raw(es=0.6207,lb=0.3382,ub=1.1391,m1=355,m2=366,dec=4,measure="or")
#'estim_or1
#'
#'The results indicate that the values of a and c are either respectively 18 and 29 (solution 1) or 327 and 348 (solution 2).
#'These two solutions are presented since the quadratic formula had been applied to estimate raw data.
#'The correct solution should be assumed based on prior clinical knowledge (i.e., whether the outcome event is a probable one or not).
#'The real values are actually a=18 and c=29.
#'
#'
#'## **Low effect size precision and information on the number of participants developing the outcome of event:**
#'
#'Consider a primary study presenting the Odds Ratio as 0.6, with a 95% confidence interval of 0.3-1.1
#'(worked example 4.2. of C. Di Pietrantonj (Statistics in Medicine 2006;25:2299-2322)).
#'Let us assume that (i) We know that there are 355 participants in the treated/exposed group, and 366 participants in the control group;
#'(ii) we know that, overall, 47 participants developed the outcome of event.
#'Raw data can be estimated by:
#'
#'estim_or2 <- estim_raw(es=0.6,lb=0.3,ub=1.1,m1=355,m2=366,e1=47,dec=1,measure="or")
#'estim_or2
#'
#'The results indicate that the values of a and c are either respectively 17 and 30 or 18 and 29.
#'The real values are actually a=18 and c=29.
#'
#'
#'## **High effect size precision and information on the number of participants developing the outcome of event:**
#'
#'Consider a primary study presenting the Risk Difference as -7.83%, with a 95% confidence interval of -13.87;-1.8%
#'(worked example 2.2. of C. Di Pietrantonj (Statistics in Medicine 2006;25:2299-2322)).
#'Let us assume that (i) We know that there are 373 participants in the treated/exposed group, and 357 participants in the control group;
#'(ii) we know that, overall, 163 participants developed the outcome of event.
#'Raw data can be estimated by:
#'
#'estim_rd <- estim_raw(es=-0.0783,lb=-0.1387,ub=-0.018,m1=373,m2=357,e1=163,dec=4,measure="rd")
#'estim_rd
#'
#'The results indicate that the values of a and c are respectively 69 and 94.
#'These actually correspond to the real values of a and c.
#'
#'@return
#'\enumerate{
#'\item Estimates from risk ratios: A dataframe. If there is information on the number of participants developing the outcome event (e1), this dataframe will list all sets of results in which a+c=e1. If no such information is provided, this dataframe will list a point estimate for each cell (calculated based on the exact input values), as well as a minimum and a maximum estimate.
#'\item Estimates from risk differences and odds ratios with known number of participants developing the outcome event (e1): A dataframe containing all sets of results in which a+c=e1.
#'\item Estimates from risk differences and odds ratios with unknown number of participants developing the outcome event (e1): A list consisting of two dataframes - solution1 (presenting the results of the first solution of the quadratic formula) and solution2 (presenting the results of the second solution of the quadratic formula). Each of these dataframes will list a point estimate for each cell (calculated based on the exact input values), as well as a minimum and a maximum estimate.
#'}
#'
#'@examples
#'
#'estim_rr <- estim_raw(es=0.6,lb=0.4,ub=0.9,m1=352,m2=376,dec=1,measure="rr")
#'estim_rr
#'
#'estim_or1 <- estim_raw(es=0.6207,lb=0.3382,ub=1.1391,m1=355,m2=366,dec=4,measure="or")
#'estim_or1
#'
#'estim_or2 <- estim_raw(es=0.6,lb=0.3,ub=1.1,m1=355,m2=366,e1=47,dec=1,measure="or")
#'estim_or2
#'
#'estim_rd <- estim_raw(es=-0.0783,lb=-0.1387,ub=-0.018,m1=373,m2=357,e1=163,dec=4,measure="rd")
#'estim_rd
#'
#'@section References
#'\enumerate{
#'\item Di Pietrantonj C. Four-fold table cell frequencies imputation in meta analysis. Statistics in Medicine 2006;25:2299-2322
#'}
#'
#'@importFrom dplyr distinct filter rename %>%
#'
#'@export
estim_raw <- function(es,lb,ub,m1,m2,e1,dec=1,measure=c("or","rd","rr")) {
  if(measure=="or"){
    if(missing(e1)){
      es_seq <- seq(from=es-(5*10^(-dec-1)),to=es+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
      lb_seq <- seq(from=lb-(5*10^(-dec-1)),to=lb+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
      ub_seq <- seq(from=ub-(5*10^(-dec-1)),to=ub+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
      framex <- expand.grid(es_seq,lb_seq,ub_seq)

      framex$Var4 <- ((log(framex[,3])-log(framex[,2]))/3.92)*((log(framex[,3])-log(framex[,2]))/3.92)
      framex$Var5 <- ((1-framex[,1])^2)+(framex[,1]*m2*framex[,4])
      framex$Var6 <- (framex[,1]*m1)*((2*(1-framex[,1]))-m2*framex[,4])
      framex$Var7 <- (framex[,1]*m1)*((framex[,1]*m1)+m2)

      framex$solution1_a <- round((-framex[,6]-sqrt((framex[,6]^2)-(4*framex[,5]*framex[,7])))/(2*framex[,5]),0)
      framex$solution1_b <- round((m1-framex[,8]),0)
      framex$solution1_c <- round((framex[,8]*m2)/((framex[,1]*m1)+(framex[,8]*(1-framex[,1]))),0)
      framex$solution1_d <- round((m2-framex[,10]),0)

      framex$solution2_a <- round((-framex[,6]+sqrt((framex[,6]^2)-(4*framex[,5]*framex[,7])))/(2*framex[,5]),0)
      framex$solution2_b <- round((m1-framex[,12]),0)
      framex$solution2_c <- round((framex[,12]*m2)/((framex[,1]*m1)+(framex[,12]*(1-framex[,1]))),0)
      framex$solution2_d <- round((m2-framex[,14]),0)

      subs <- framex[,8:15]
      subs1 <- distinct(subs)

      point_se <- ((log(ub)-log(lb))/3.92)*((log(ub)-log(lb))/3.92)
      point_alpha <- ((1-es)^2)+(es*m2*point_se)
      point_beta <- (es*m1)*((2*(1-es))-m2*point_se)
      point_gamma <- (es*m1)*((es*m1)+m2)

      pointsolution1_a <- round((-point_beta-sqrt((point_beta^2)-(4*point_alpha*point_gamma)))/(2*point_alpha),0)
      pointsolution1_b <- round((m1-pointsolution1_a),0)
      pointsolution1_c <- round((pointsolution1_a*m2)/((es*m1)+(pointsolution1_a*(1-es))),0)
      pointsolution1_d <- round((m2-pointsolution1_c),0)

      pointsolution2_a <- round((-point_beta+sqrt((point_beta^2)-(4*point_alpha*point_gamma)))/(2*point_alpha),0)
      pointsolution2_b <- round((m1-pointsolution2_a),0)
      pointsolution2_c <- round((pointsolution2_a*m2)/((es*m1)+(pointsolution2_a*(1-es))),0)
      pointsolution2_d <- round((m2-pointsolution2_c),0)

      subs3 <- data.frame("a"=c(pointsolution1_a,min(subs1[,1]),max(subs1[,1])),"b"=c(pointsolution1_b,min(subs1[,2]),max(subs1[,2])),"c"=c(pointsolution1_c,min(subs1[,3]),max(subs1[,3])),"d"=c(pointsolution1_d,min(subs1[,4]),max(subs1[,4])),row.names =c("Point estimates","Minimum values","Maximum values"))
      subs4 <- data.frame("a"=c(pointsolution2_a,min(subs1[,5]),max(subs1[,5])),"b"=c(pointsolution2_b,min(subs1[,6]),max(subs1[,6])),"c"=c(pointsolution2_c,min(subs1[,7]),max(subs1[,7])),"d"=c(pointsolution2_d,min(subs1[,8]),max(subs1[,8])),row.names =c("Point estimates","Minimum values","Maximum values"))
      list(solution1=subs3,solution2=subs4)

    }else{
      es_seq <- seq(from=es-(5*10^(-dec-1)),to=es+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
      lb_seq <- seq(from=lb-(5*10^(-dec-1)),to=lb+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
      ub_seq <- seq(from=ub-(5*10^(-dec-1)),to=ub+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
      framex <- expand.grid(es_seq,lb_seq,ub_seq)

      framex$Var4 <- ((log(framex[,3])-log(framex[,2]))/3.92)*((log(framex[,3])-log(framex[,2]))/3.92)
      framex$Var5 <- ((1-framex[,1])^2)+(framex[,1]*m2*framex[,4])
      framex$Var6 <- (framex[,1]*m1)*((2*(1-framex[,1]))-m2*framex[,4])
      framex$Var7 <- (framex[,1]*m1)*((framex[,1]*m1)+m2)

      framex$solution1_a <- round((-framex[,6]-sqrt((framex[,6]^2)-(4*framex[,5]*framex[,7])))/(2*framex[,5]),0)
      framex$solution1_b <- round((m1-framex[,8]),0)
      framex$solution1_c <- round((framex[,8]*m2)/((framex[,1]*m1)+(framex[,8]*(1-framex[,1]))),0)
      framex$solution1_d <- round((m2-framex[,10]),0)

      framex$solution2_a <- round((-framex[,6]+sqrt((framex[,6]^2)-(4*framex[,5]*framex[,7])))/(2*framex[,5]),0)
      framex$solution2_b <- round((m1-framex[,12]),0)
      framex$solution2_c <- round((framex[,12]*m2)/((framex[,1]*m1)+(framex[,12]*(1-framex[,1]))),0)
      framex$solution2_d <- round((m2-framex[,14]),0)

      subs <- framex[,8:11]
      subs1 <- distinct(subs)
      subs1 <- subs1 %>% rename(a = solution1_a, b=solution1_b, c=solution1_c, d=solution1_d)
      subs2 <- framex[,12:15]
      subs3 <- distinct(subs2)
      subs3 <- subs3 %>% rename(a = solution2_a, b=solution2_b, c=solution2_c, d=solution2_d)
      subs4 <- subs1%>%filter(a+c==e1)
      subs5 <- subs3%>%filter(a+c==e1)
      subs6 <- rbind(subs4,subs5)
      prefixo <- "Solution"
      nro <- nrow(subs6)
      suffixo <- seq(from=1,to=nro)
      if (nrow(subs6)>0) (data.frame(subs6, row.names = paste(prefixo,suffixo))) else (message("No raw data satisfying the provided estimates. Please recheck your input values."))

    }
  }else{
    if(measure=="rd"){
      if(missing(e1)){
        es_seq <- seq(from=es-(5*10^(-dec-1)),to=es+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
        lb_seq <- seq(from=lb-(5*10^(-dec-1)),to=lb+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
        ub_seq <- seq(from=ub-(5*10^(-dec-1)),to=ub+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
        framex <- expand.grid(es_seq,lb_seq,ub_seq)

        framex$Var4 <- ((framex[,3]-framex[,2])/3.92)*((framex[,3]-framex[,2])/3.92)
        framex$Var5 <- m1+m2
        framex$Var6 <- -m1*(m1*(1+2*framex[,1])+m2)
        framex$Var7 <- ((m1)^3)*((framex[,1]*(1+framex[,1]))+(framex[,4]*m2))

        framex$solution1_a <- round((-framex[,6]-sqrt((framex[,6]^2)-(4*framex[,5]*framex[,7])))/(2*framex[,5]),0)
        framex$solution1_b <- round((m1-framex[,8]),0)
        framex$solution1_c <- round(((framex[,8]*(m2/m1))-(framex[,1]*m2)),0)
        framex$solution1_d <- round((m2-framex[,10]),0)

        framex$solution2_a <- round((-framex[,6]+sqrt((framex[,6]^2)-(4*framex[,5]*framex[,7])))/(2*framex[,5]),0)
        framex$solution2_b <- round((m1-framex[,12]),0)
        framex$solution2_c <- round(((framex[,12]*(m2/m1))-(framex[,1]*m2)),0)
        framex$solution2_d <- round((m2-framex[,14]),0)

        subs <- framex[,8:15]
        subs1 <- distinct(subs)

        point_se <- ((ub-lb)/3.92)*((ub-lb)/3.92)
        point_alpha <- m1+m2
        point_beta <- -m1*(m1*(1+2*es)+m2)
        point_gamma <- ((m1)^3)*((es*(1+es))+(point_se*m2))

        pointsolution1_a <- round((-point_beta-sqrt((point_beta^2)-(4*point_alpha*point_gamma)))/(2*point_alpha),0)
        pointsolution1_b <- round((m1-pointsolution1_a),0)
        pointsolution1_c <- round(((pointsolution1_a*(m2/m1))-(es*m2)),0)
        pointsolution1_d <- round((m2-pointsolution1_c),0)

        pointsolution2_a <- round((-point_beta+sqrt((point_beta^2)-(4*point_alpha*point_gamma)))/(2*point_alpha),0)
        pointsolution2_b <- round((m1-pointsolution2_a),0)
        pointsolution2_c <- round(((pointsolution2_a*(m2/m1))-(es*m2)),0)
        pointsolution2_d <- round((m2-pointsolution2_c),0)

        subs3 <- data.frame("a"=c(pointsolution1_a,min(subs1[,1]),max(subs1[,1])),"b"=c(pointsolution1_b,min(subs1[,2]),max(subs1[,2])),"c"=c(pointsolution1_c,min(subs1[,3]),max(subs1[,3])),"d"=c(pointsolution1_d,min(subs1[,4]),max(subs1[,4])),row.names =c("Point estimates","Minimum values","Maximum values"))
        subs4 <- data.frame("a"=c(pointsolution2_a,min(subs1[,5]),max(subs1[,5])),"b"=c(pointsolution2_b,min(subs1[,6]),max(subs1[,6])),"c"=c(pointsolution2_c,min(subs1[,7]),max(subs1[,7])),"d"=c(pointsolution2_d,min(subs1[,8]),max(subs1[,8])),row.names =c("Point estimates","Minimum values","Maximum values"))
        list(solution1=subs3,solution2=subs4)


      }else{
        es_seq <- seq(from=es-(5*10^(-dec-1)),to=es+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
        lb_seq <- seq(from=lb-(5*10^(-dec-1)),to=lb+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
        ub_seq <- seq(from=ub-(5*10^(-dec-1)),to=ub+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
        framex <- expand.grid(es_seq,lb_seq,ub_seq)

        framex$Var4 <- ((framex[,3]-framex[,2])/3.92)*((framex[,3]-framex[,2])/3.92)
        framex$Var5 <- m1+m2
        framex$Var6 <- -m1*(m1*(1+2*framex[,1])+m2)
        framex$Var7 <- ((m1)^3)*((framex[,1]*(1+framex[,1]))+(framex[,4]*m2))

        framex$solution1_a <- round((-framex[,6]-sqrt((framex[,6]^2)-(4*framex[,5]*framex[,7])))/(2*framex[,5]),0)
        framex$solution1_b <- round((m1-framex[,8]),0)
        framex$solution1_c <- round(((framex[,8]*(m2/m1))-(framex[,1]*m2)),0)
        framex$solution1_d <- round((m2-framex[,10]),0)

        framex$solution2_a <- round((-framex[,6]+sqrt((framex[,6]^2)-(4*framex[,5]*framex[,7])))/(2*framex[,5]),0)
        framex$solution2_b <- round((m1-framex[,12]),0)
        framex$solution2_c <- round(((framex[,12]*(m2/m1))-(framex[,1]*m2)),0)
        framex$solution2_d <- round((m2-framex[,14]),0)

        subs <- framex[,8:11]
        subs1 <- distinct(subs)
        subs1 <- subs1 %>% rename(a = solution1_a, b=solution1_b, c=solution1_c, d=solution1_d)
        subs2 <- framex[,12:15]
        subs3 <- distinct(subs2)
        subs3 <- subs3 %>% rename(a = solution2_a, b=solution2_b, c=solution2_c, d=solution2_d)
        subs4 <- subs1%>%filter(a+c==e1)
        subs5 <- subs3%>%filter(a+c==e1)
        subs6 <- rbind(subs4,subs5)
        prefixp <- "Solution"
        nrp <- nrow(subs6)
        suffixp <- seq(from=1,to=nrp)
        if (nrow(subs6)>0) (data.frame(subs6, row.names = paste(prefixp,suffixp))) else (message("No raw data satisfying the provided estimates. Please recheck your input values."))
      }
    }else{
      if(measure=="rr"){
        if(missing(e1)){
          es_seq <- seq(from=es-(5*10^(-dec-1)),to=es+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
          lb_seq <- seq(from=lb-(5*10^(-dec-1)),to=lb+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
          ub_seq <- seq(from=ub-(5*10^(-dec-1)),to=ub+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
          framex <- expand.grid(es_seq,lb_seq,ub_seq)

          framex$Var4 <- ((log(framex[,3])-log(framex[,2]))/3.92)*((log(framex[,3])-log(framex[,2]))/3.92)

          framex$a <- round(((m2+(framex[,1]*m1))/(m2*(framex[,4]+((m2+m1)/(m1*m2))))),0)
          framex$b <- round((m1-framex[,5]),0)
          framex$c <- round(((framex[,5]*m2)/(framex[,1]*m1)),0)
          framex$d <- round((m2-framex[,7]),0)

          subs <- framex[,5:8]
          subs1 <- distinct(subs)

          point_se <- ((log(ub)-log(lb))/3.92)*((log(ub)-log(lb))/3.92)

          pointsolution1_a <- round(((m2+(es*m1))/(m2*(point_se+((m2+m1)/(m1*m2))))),0)
          pointsolution1_b <- round((m1-pointsolution1_a),0)
          pointsolution1_c <- round(((pointsolution1_a*m2)/(es*m1)),0)
          pointsolution1_d <- round((m2-pointsolution1_c),0)

          subs3 <- data.frame("a"=c(pointsolution1_a,min(subs1[,1]),max(subs1[,1])),"b"=c(pointsolution1_b,min(subs1[,2]),max(subs1[,2])),"c"=c(pointsolution1_c,min(subs1[,3]),max(subs1[,3])),"d"=c(pointsolution1_d,min(subs1[,4]),max(subs1[,4])),row.names =c("Point estimate","Minimum value","Maximum value"))
          subs3


        }else{
          es_seq <- seq(from=es-(5*10^(-dec-1)),to=es+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
          lb_seq <- seq(from=lb-(5*10^(-dec-1)),to=lb+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
          ub_seq <- seq(from=ub-(5*10^(-dec-1)),to=ub+(4*10^(-dec-1)),by=(1*10^(-dec-1)))
          framex <- expand.grid(es_seq,lb_seq,ub_seq)

          framex$Var4 <- ((log(framex[,3])-log(framex[,2]))/3.92)*((log(framex[,3])-log(framex[,2]))/3.92)

          framex$a <- round(((m2+(framex[,1]*m1))/(m2*(framex[,4]+((m2+m1)/(m1*m2))))),0)
          framex$b <- round((m1-framex[,5]),0)
          framex$c <- round(((framex[,5]*m2)/(framex[,1]*m1)),0)
          framex$d <- round((m2-framex[,7]),0)

          subs <- framex[,5:8]
          subs1 <- distinct(subs)
          subs2 <- subs1%>%filter(a+c==e1)
          prefix <- "Solution"
          nr <- nrow(subs2)
          suffix <- seq(from=1,to=nr)
          if (nrow(subs2)>0) (data.frame(subs2, row.names = paste(prefix,suffix))) else (message("No raw data satisfying the provided estimates. Please recheck your input values."))
        }
      }else{
        warning("Please provide a valid effect size measure")
      }
    }
  }
}
