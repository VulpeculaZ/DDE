library(CollocInfer)
SEIRtimes = SEIRtimes
SEIRdata = SEIRdata
data = cbind(matrix(NA,length(SEIRdata),2),SEIRdata)
logdata = log(data)
SEIRvarnames = SEIRvarnames
SEIRparnames = SEIRparnames
SEIRpars = SEIRpars
SEIRfn = make.SEIR()
