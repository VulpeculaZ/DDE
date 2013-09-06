times = seq(0,20,0.5)
FHN.fn = make.fhn()
Treated <- Puromycin[Puromycin$state == "treated", ]
weighted.MM <- function(resp, conc, Vm, K)
{
    ## Purpose: exactly as white book p. 451 -- RHS for nls()
    ##  Weighted version of Michaelis-Menten model
    ## ----------------------------------------------------------
    ## Arguments: 'y', 'x' and the two parameters (see book)
    ## ----------------------------------------------------------
    ## Author: Martin Maechler, Date: 23 Mar 2001

    pred <- (Vm * conc)/(K + conc)
    (resp - pred) / sqrt(pred)
}

Pur.wt <- nls( ~ weighted.MM(rate, conc, Vm, K), data = Treated,
              start = list(Vm = 200, K = 0.1))
summary(Pur.wt)

res = nls(~ProfileSSE(pars, allpars, times, data, coefs, lik, proc,
in.meth, control.in, active),
data = list(allpars=pars, times=times, data=data, coefs=ncoefs,
lik=lik, proc=proc,
in.meth=in.meth,control.in=control.in,active=active),
start = list(pars=pars[active]),trace=control.out$trace,control=control.out)
