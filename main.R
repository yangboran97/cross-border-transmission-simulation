#Main proceeding code
setwd("XXX")
q_e=1#proportion of infectious outliers of particular 7day(0.309)14days(0.105)21days(0.045)28days(0.023)
t_r=1#Tourist reduction
source("Functional code.R")
source("Vaccination code.R")
source("Parameter code.R")
this_C <- C/scale_50 #R0 settings
ptm <- proc.time()  
for (i in seq(0, 100, by =1)){
  j <- i/100
  list_health_workers[[paste0(i)]] <- my_sim(this_C, j, "health_workers", num_perday, v_e_type, this_v_e)
}
time<-0:365
options(scipen=200)

