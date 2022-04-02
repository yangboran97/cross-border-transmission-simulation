#Basic simulation frame
framework = function(C, percent_vax, strategy, num_perday, v_e_type, v_e = v_e_constant,
                  u = u_var, sero = sero_none, sp = 1, se = 0,  syn_sero_compartments = NA){
  E_0 <- Sv_0 <- R_0 <- D_0 <- rep(0, num_groups)
 
  I_0 <- rep(0, num_groups) 
  I_0[11]<-1
  S_0 <- N_i - I_0 - E_0 - R_0
  S_0[11]<-oversea_input
  E_0[11] <- S_0[11]*rate_I
  S_0[11] <- 0
  S_0<-S_0*(1-v_ac)
  if (strategy == "all"){ 
    groups <- 1:10 
  } else if (strategy == "health_workers"){ # Vaccination was assumed prioritizing health workers
    groups <- 10
  } else if (strategy == "normal") { 
    groups <- 1:9
  } 
  people_to_vax <- sum(S_0[groups])
  vax_proportion <- rep(0, num_groups)
  vax_proportion[groups] <-S_0[groups]/people_to_vax
  vax_supply <- percent_vax*pop_total*0.72533
  if (num_perday == 1){ 
    nvax <- vax_supply 
    vax_distribution <- nvax*vax_proportion[1:11]
    S <- S_0
    E <- E_0
    R <- R_0
    vax_eligible <- S
    if (any(vax_distribution > vax_eligible)){ 
      if (!all(vax_distribution[groups] > vax_eligible[groups])){ 
        temp <- vax_distribution
        temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
        leftover_vax <- sum(vax_distribution - temp)
        full_groups <- 1:11
        full_groups <- full_groups[vax_distribution > vax_eligible]
        leftover_groups <- groups[!groups %in% full_groups]
        people_to_vax <- sum(vax_eligible[leftover_groups])
        vax_proportion <- rep(0, num_groups)
        vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax
        vax_leftover_dist <- leftover_vax*vax_proportion
        vax_distribution <- temp + vax_leftover_dist
      }
      if (any(vax_distribution > vax_eligible)){ 
        temp <- vax_distribution
        temp[vax_distribution > vax_eligible] <- vax_eligible[vax_distribution > vax_eligible]
        leftover_vax <- sum(vax_distribution - temp)
        leftover_groups <- 1:11
        leftover_groups <- leftover_groups[!leftover_groups %in% groups]
        people_to_vax <- sum(vax_eligible[leftover_groups])
        vax_proportion <- rep(0, num_groups)
        vax_proportion[leftover_groups] <- vax_eligible[leftover_groups]/people_to_vax 
        vax_leftover_dist <- leftover_vax*vax_proportion
        vax_distribution <- temp + vax_leftover_dist
      }
    }
    
    alpha <- as.matrix(vax_distribution)/(as.matrix(vax_eligible))
    alpha[alpha == Inf] <- 0
    alpha[is.nan(alpha)] <- 0
    if (v_e_type == "aorn") {
      S_0  <- S - as.matrix(S*alpha*v_e) 
      Sv_0 <- as.matrix(S*alpha*v_e)
    }
    vax_supply <- 0 
  }
  compartments_initial <- c(S_0,Sv_0,E_0,I_0,R_0,D_0,vax_supply)
  if (length(syn_sero_compartments) > 1){
    compartments_initial <- c(as.numeric(syn_sero_compartments), vax_supply)
  }            
  parameters = list(u=u, C=C, d_E=d_E, d_I=d_I, v_e=v_e, v_e_type = v_e_type, num_groups = num_groups, 
                    N_i = N_i, num_perday=num_perday, vax_proportion=vax_proportion, groups=groups, pop_total=pop_total)
  t <- seq(from=0, to=365, by=1) 
  event_times <- seq(from=0, to=floor(vax_supply/(num_perday*pop_total)), by=1)
  df <- as.data.frame(deSolve::lsoda(y=compartments_initial, times=t, func=calculate_derivatives, parms=parameters, events=list(func=move_vaccinated_event, time=event_times)))
  names(df) <- c("time", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9","S10","S11",
                 "Sv1", "Sv2", "Sv3", "Sv4", "Sv5", "Sv6", "Sv7", "Sv8", "Sv9", "Sv10","Sv11",
                 "E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9","E10","E11",
                 "I1", "I2", "I3", "I4", "I5", "I6", "I7", "I8", "I9","I10","I11",
                 "R1", "R2", "R3", "R4", "R5", "R6", "R7", "R8", "R9","R10","R11",
                 "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9","D10","D11")
  
  return(df)
}
