get_rt = function(beta,
                  phi, # av contact matrix
                  pop,
                  pop_ratio_mat,
                  sigma,
                  tau,
                  rho_A,
                  rho_Y,
                  gamma_A,
                  gamma_Y,
                  omega_A,
                  omega_Y,
                  omega_P,
                  S_prop) {
  
  if( is.na(beta) | sum(is.na(S_prop)>0) )
    return (NA)
  n_age = nrow(phi)
  n_c = 5  # [E,PA,PY,IA,IY]
  
  # zero matrix useful recycle below
  Zm = matrix(0, n_age, n_age)
  
  # Proportion of susceptibles
  S_prop_mat = matrix(S_prop,nrow = 5,ncol = 5,byrow=FALSE)
  
  
  # tiled av contact matrix
  big_C = matrix(0, n_c * n_age, n_c * n_age)
  for (i in 1:n_c)
    for (j in 1:n_c) {
      start_row = n_c * (i - 1) + 1
      end_row = n_c * i
      start_col = n_c * (j - 1) + 1
      end_col = n_c * j
      big_C[start_row:end_row, start_col:end_col] = phi 
    }
  
  # Transmission and Transition Matrices
  
  block11 = beta * Zm                                                                     * S_prop_mat  
  block12 = beta * pop_ratio_mat * matrix( omega_P * omega_A, n_age, n_age, byrow = TRUE) * S_prop_mat 
  block13 = beta * pop_ratio_mat * matrix( omega_P * omega_Y, n_age, n_age, byrow = TRUE) * S_prop_mat 
  block14 = beta * pop_ratio_mat * matrix(           omega_A, n_age, n_age              ) * S_prop_mat 
  block15 = beta * pop_ratio_mat * matrix(           omega_Y, n_age, n_age              ) * S_prop_mat 
  row1 = cbind(block11, block12, block13, block14, block15)
  row2 = matrix(0, 4 * n_age, 5 * n_age)
  T_mat = big_C * rbind(row1, row2)
  
  I = diag(1, n_age); O = Zm
  row1 = cbind(          -sigma *I,        O,         O,          O,          O)  # E
  row2 = cbind( (1 -tau)* sigma *I, -rho_A*I,         O,          O,          O)  # P_A
  row3 = cbind(     tau * sigma *I,        O, -rho_Y *I,          O,          O)  # P_Y
  row4 = cbind(                  O,  rho_A*I,         O, -gamma_A*I,          O)  # I_A 
  row5 = cbind(                  O,        O,  rho_Y *I,          O, -gamma_Y*I)  # I_Y
  #                              
  Sigma_mat = rbind(row1, row2, row3, row4, row5)
  
  
  # R0 computation
  K_L = - T_mat %*% solve(Sigma_mat)
  r0 = max(as.numeric(eigen(K_L)$values))
  
  r0
}


