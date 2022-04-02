#Basic paramter settings
library(tidyverse) 
library(deSolve) # ode solver
library(gridExtra)
library(RColorBrewer)
library(wesanderson)
library(egg)
library(foreach)
library(doParallel)
library(ggplot2)
library(gplots)
d_E <- 1/4.4 
d_I <- 1/5 
C_o <- read.csv("XXX.csv",header=T)# Upload country-specific contact matrix
C <- as.matrix(C_o[1:11,2:12])
age_demo   <- read.csv("XXX.csv",header=T)# Upload country-specific demographics
pop_total  <- age_demo[12]
age_demo   <- age_demo[1:11]

N_i        <- pop_total*age_demo
rate_I <-0.00016 #proportion of infections among all oversea input
oversea_input<-14797000
Long_term<-0.02
Dynamic<-0.27
oversea_input<-oversea_input*Long_term #border strategy


num_groups <- length(age_demo) # num  groups

IFR        <- c(9.807703e-04, 3.277686e-03, 1.095386e-02, 3.660727e-02, 1.223397e-01, 4.088531e-01, 1.366367e+00,
                4.566330e+00, 1.526045e+01,0.6,9.807703e-04) 
IFR        <- IFR/100 # as decimal


u_var      <- c(0.4, 0.38, 0.79, 0.86, 0.8, 0.82, 0.88, 0.74, 0.74,0.04,0)


C[1:9,11]<-C[1:9,11]*q_e
C[1:9,11]<-C[1:9,11]*t_r

num_perday <- 1
list_all<-list_health_workers <- list_normal <-  vector(mode = "list")

scale_10 <- scale_u_for_R0(u_var, C, 1.0)
scale_11 <- scale_u_for_R0(u_var, C, 1.1)
scale_115 <- scale_u_for_R0(u_var, C, 1.15)
scale_12 <- scale_u_for_R0(u_var, C, 1.2)
scale_15 <- scale_u_for_R0(u_var, C, 1.5)
scale_16 <- scale_u_for_R0(u_var, C, 1.6)
scale_17 <- scale_u_for_R0(u_var, C, 1.7)
scale_18 <- scale_u_for_R0(u_var, C, 1.8)
scale_19 <- scale_u_for_R0(u_var, C, 1.9)
scale_25 <- scale_u_for_R0(u_var, C, 2.5)
scale_26 <- scale_u_for_R0(u_var, C, 2.6)
scale_20 <- scale_u_for_R0(u_var, C, 2.0)
scale_25 <- scale_u_for_R0(u_var, C, 2.5)
scale_30 <- scale_u_for_R0(u_var, C, 3.0)
scale_35 <- scale_u_for_R0(u_var, C, 3.5)
scale_40 <- scale_u_for_R0(u_var, C, 4.0)
scale_45 <- scale_u_for_R0(u_var, C, 4.5)
scale_50 <- scale_u_for_R0(u_var, C, 5.0)

C_10 <- C/scale_10
C_115 <- C/scale_115
C_12  <- C/scale_12
C_15  <- C/scale_15
C_16  <- C/scale_16
C_17  <- C/scale_17
C_18  <- C/scale_18
C_19  <- C/scale_19
C_20  <- C/scale_20
C_25  <- C/scale_25
C_26  <- C/scale_26
C_20  <- C/scale_20
C_30  <- C/scale_30
C_35  <- C/scale_35
C_45  <- C/scale_45
C_40  <- C/scale_40
C_50  <- C/scale_50

