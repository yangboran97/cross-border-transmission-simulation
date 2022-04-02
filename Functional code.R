when_strat_switch = function(list_df, groups){
  # find the vaccine supply where you switch from strategy to vaccinating everyone else
  # input: list_df is a list with df for a strategy for vax supply 0-50
  #        groups: vector index of the strategy
  # output: vax supply where switch occurs
  
  other_groups <- 1:11
  other_groups <- other_groups[!other_groups %in% groups]
  x_switch <- -1
  
  for (i in 1:length(list_df)){
    temp <- list_df[[i]]
    if (any(temp[dim(temp)[1], other_groups + 10] > 0)){
      x_switch <- i-1 
      break
    }
  }
  return(x_switch)
}



get_v_e = function(p, y0, hinge_age){
  # INPUT: p is the final v_e % (as a decimal) 
  #       i.e. p = 1 -> perfect vaccine at all ages
  #       y0 is the baseline v_e 
  #       hinge_age is the age group that v_e begins decreasing after
  # ASSUMPTIONS: perfect v_e up to hinge age then linear decline to p
  # OUTPUT: vector v_e 
  
  x <- c(0, 10, 20, 30, 40, 50, 60, 70, 80,60,60)
  y <- y0 + ((p-y0)/(80-hinge_age))*(x - hinge_age)
  index <- match(c(hinge_age), x)
  y[1:index] <- y0
  
  y
}


scale_u_for_R0 = function(u, C, wanted_R0){
  scalehigh <- 10
  scalelow <- 100
  
  R0_high <- compute_R0(u/scalehigh, C)
  R0_low <- compute_R0(u/scalelow, C)
  
  if (R0_high < wanted_R0 || R0_low > wanted_R0){
    print("Error: wanted_R0 is not in the original range of R0 values")
    return(-1)
  }
  
  running = TRUE
  while(running){
    scale_next = (scalehigh + scalelow)/2
    temp <- compute_R0(u/scale_next, C)
    if (temp > wanted_R0){
      scalehigh <- scale_next
    } else {
      scalelow <- scale_next
    }
    
    if (abs(temp - wanted_R0) < 0.0001){
      return(scale_next)
    }
  }
}



compute_R0 = function(u, C){
  gamma <- 1/5 # recovery period (I -> R), ref: Davies
  
  # Davies NGM
  Du <- diag(u, 11)
  Dy <- diag(1/gamma, 11)
  
  NGM <- Du %*% C %*% Dy
  R0  <- abs(eigen(NGM)$values[1])
}


calculate_derivatives = function(t, x, parameters){
  # x is a vector of length (# model compartment types)*(# age groups)
  # S, E, I, R etc. are vectors of length num_groups
  num_compartment <- 6
  num_groups <- (length(x)-1)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  E    <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  I    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  R    <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  D    <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  vax_supply <- x[6*num_groups+1]
  
  S[S < .Machine$double.eps] <- 0
  Sv[Sv < .Machine$double.eps] <- 0
  
  E[E < .Machine$double.eps] <- 0
  
  I[I < .Machine$double.eps] <- 0
  
  R[R < .Machine$double.eps] <- 0
  
  D[D < .Machine$double.eps] <- 0
  
  u <- parameters$u#È¡³öu
  C <- parameters$C
  d_E <- parameters$d_E
  d_I <- parameters$d_I
  v_e <- parameters$v_e
  v_e_type <- parameters$v_e_type
  N_i[11]<-oversea_input
  lambda <- C%*%(I/(N_i-D))*u
  
  if (v_e_type == "aorn") {
    # all-or-nothing
    dSv <- rep(0, num_groups)
    
  }
  
  dS <- -(S*lambda)
  
  
  dE  <- S*lambda - d_E*E
  dI  <- E*d_E - I*d_I
  dR  <- I*d_I*(1-IFR)
  dD  <- I*d_I*IFR
  
  dvax_supply <- 0
  
  out <- c(dS,dSv,dE,dI,dR,dD,dvax_supply)
  list(out)
}


move_vaccinated_event <- function(t, x, parms){
  v_e <- parms$v_e
  v_e_type <- parms$v_e_type
  num_perday <- parms$num_perday
  vax_proportion <- parms$vax_proportion
  groups <- parms$groups
  pop_total <- parms$pop_total
  
  # move those who are vaccinated in a given day
  num_compartment <- 6
  num_groups <- (length(x)-1)/num_compartment
  S    <- as.matrix(x[1:num_groups])
  Sv   <- as.matrix(x[(1*num_groups+1):(2*num_groups)])
  E    <- as.matrix(x[(2*num_groups+1):(3*num_groups)])
  I    <- as.matrix(x[(3*num_groups+1):(4*num_groups)])
  R    <- as.matrix(x[(4*num_groups+1):(5*num_groups)])
  D    <- as.matrix(x[(5*num_groups+1):(6*num_groups)])
  vax_supply <- x[6*num_groups+1]
  
  if (vax_supply >= num_perday*pop_total){
    nvax <- num_perday*pop_total
    vax_supply <- vax_supply - num_perday*pop_total
  } else {
    nvax <- vax_supply
    vax_supply <- 0
  }
  
  vax_distribution <- nvax*vax_proportion
  vax_eligible <- S
  if (any(vax_distribution > vax_eligible)){ 
    # make sure everyone in the specificed age groups are vaccinated
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
    #vax_distribution[vax_distribution > (S+E+R)] <- S[vax_distribution > (S+E+R)] + E[vax_distribution > (S+E+R)] + R[vax_distribution > (S+E+R)]
    #distribute vaccines to all other age groups if there's doses left over after strategy specified
    #age groups
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
    # don't go over the number eligible
    if (any(vax_distribution > vax_eligible)){
      vax_distribution[vax_distribution > vax_eligible] = vax_eligible
    }
  }
  
  alpha <- vax_distribution/vax_eligible
  alpha[alpha == Inf] <- 0 # no people in S,E,R
  alpha[is.nan(alpha)] <- 0 # no vax left and no people in S,E,R
  if(any(alpha > 1)){print("WARNING: alpha > 1 in move_vaccinated")
    alpha[alpha>1] <- 0} # more vaccines avaliable than eligible people
  
  if (v_e_type == "aorn") {
    # all-or-nothing
    dS  <- -S*alpha 
    dSv <- (S*alpha*v_e)
  }
  
  # update compartments
  S    <- S + dS
  Sv   <- Sv + dSv
  
  
  # output updated compartments
  out <- c(S,Sv,E,I,R,D,vax_supply)
}




get_reduction_df = function(outcome){
  total <- rep(NA, 202)
  
  if (outcome == "cases"){
    total <- compile_total_cases(total)
    baseline <- compute_total_cases_new(list_normal$`0`)
  } else if (outcome == "deaths"){
    total <- compile_total_deaths(total)
    baseline <- compute_total_deaths_new(list_normal$`0`)
  } 
  
  num_strategies <- 2
  vax_avail <- c(rep(seq(0, 100, by = 1), num_strategies))
  num_per_list <- 101
  strat <- c(rep("health_workers", num_per_list), rep("normal", num_per_list))
  variable <-  c(rep("constant", num_per_list))
  
  
  
  baseline <- c(rep(baseline, num_per_list))
  reduction <- (1-(total/baseline))*100
  
  df <- data.frame(vax_avail, strat, reduction, variable)
}

get_rowindex <- function(df, rName) {
  which(match(colnames(df), rName) == 1)
}


compile_total_cases = function(total){
  count <- 1
  for (i in list_health_workers){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  for (i in list_normal){
    total[count] <- compute_total_cases_new(i)
    count <- count + 1
  }
  total
}



compile_total_deaths = function(total){
  count <- 1
  for (i in list_health_workers){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  for (i in list_normal){
    total[count] <- compute_total_deaths_new(i)
    count <- count + 1
  }
  
  total
}

compute_total_cases_new = function(df){
  infections <- rep(0,num_groups) # number of age groups
  
  R_index <- 46 # col number for R 
  D_index <- 57
  
  for (i in 1:num_groups) {
    # infections = total # recovered - initial recovered (seropositive)
    infections[i] <- df[dim(df)[1], R_index] - df[1, R_index] + df[dim(df)[1], D_index]
    R_index  <- R_index + 1
    D_index  <- D_index + 1
  }
  tot_infections <- sum(infections)/pop_total * 100
}

compute_total_deaths_new = function(df){
  deaths <- rep(0,num_groups)
  D_index <- 57
  
  for (i in 1:num_groups) {
    deaths[i] <- df[dim(df)[1], D_index]
    D_index <- D_index + 1
  }
  tot_deaths <- sum(deaths)/pop_total * 100
}

gather_compartments_overtime <- function(df, compartment, strat){
  final_time <- as.numeric(dim(df)[1])-1
  
  if (compartment == "I"){
    start <- "I1"
    end <- "I10"
  } else if (compartment == "R"){
    start <- "R1"
    end <- "D10"
  } else if (compartment == "D"){
    start <- "D1"
    end <- "D10"
  }
  
  
  total_df <- data.frame(time = 0:final_time,
                         percent = ((apply(df[, get_rowindex(df, start):get_rowindex(df, end)], 
                                           1, sum))/pop_total)*100,
                         strat = rep(strat, final_time +1))
}