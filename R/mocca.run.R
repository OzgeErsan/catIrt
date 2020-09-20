
######################################
# Create a wrapper function for the entire MOCCA process 
# max test length for MOCCA is 40 here 
######################################

mocca.run <- function(cutoff, method, delta, at.select, 
                  theta1, theta2, 
                  mocca.n.max, 
                  params1, params2, 
                  catIrt.object,
                  est.theta1.no.class){
  
  
  #*******************************
  # Transitioning from phase 1 to phase 2
  #*******************************
  
  # N and J are specified by the original function. 
  N <- dim(catIrt.object$full_resp)[1] # N is sample size (person)
  J <- dim(catIrt.object$full_resp)[2] # J is item bank size
 
  # Run mocca.transit function
  mocca.transit <- mocca.transit(theta2=theta2, params2=params2, catIrt.object=catIrt.object)
  
  # Get phase2 graded item responses and related info
  cat_indiv2_phase1 <- mocca.transit$cat_indiv
  
  # Get it.vec that indicates which item will be blocked during phase2
  it.vec <- mocca.transit$it.vec
  
  # Estimate theta2 for initial theta values for phase2
  init.theta2.vec <- vector(mode = "numeric", length = N)
  
  # Use catIrt built-in function to estimate theta for dimension2 after phase1 &
  # Add phase1 wleEst's sem output to cat_indiv2_phase1. This will be used during termCI after phase1
  wleEst.list <- vector(mode="list", length = N)
  for(i in 1:N){
    wleEst.list[[i]] <- wleEst(cat_indiv2_phase1[[i]]$cat_resp2, params=params2[cat_indiv2_phase1[[i]]$cat_it,], range = c(-4.5, 4.5), mod="grm")
    init.theta2.vec[i] <- cat_indiv2_phase1[[i]]$cat_theta2 <- wleEst.list[[i]]$theta
    cat_indiv2_phase1[[i]]$cat_sem <- wleEst.list[[i]]$sem
    cat_indiv2_phase1[[i]]$cat_info <- wleEst.list[[i]]$info
  }
  
  # Add dimension 2 item parameters used after phase1
  for( i in 1:N){
    cat_indiv2_phase1[[i]]$cat_params <- params2[cat_indiv2_phase1[[i]]$cat_it,]
  }
  
  # In phase2 cat_terminate, n.max will be a vector. For instance, entire MOCCA test length will be 40.
  mocca.test.len.vec = vector(mode = "numeric", length = N)
  mocca.test.len.vec <- rep(mocca.n.max, N)
  
  # Get test length from phase1 and substract it from mocca.test.len.vec. That is n.max for phase2.
  test.len.phase2 <- mocca.test.len.vec - catIrt.object$cat_length
  
  #*******************************
  #      Specify phase2 catIrt function inputs earlier
  #*******************************
  
  # Specify start and middle
  cat_start2 <- list(n.start=1, init.theta=init.theta2.vec, select="UW-FI", at=at.select, n.select=1, delta=NULL, score="WLE", range=c(-4.5, 4.5), step.size=NULL, leave.after.MLE= FALSE)
  cat_middle2 <- list(select="UW-FI", at=at.select, it.range=NULL, n.select=1, delta=NULL, score="WLE", range=c(-4.5, 4.5), expos="none")
  
  # Specify how to terminate the test
  class.term <- NULL
  if(method=="CI"){
    class.term <- list(method=method, bounds=cutoff, categ = c(0,2), alpha=0.1, beta =0.1, conf.level=0.90, indeterminate=TRUE)
  } else {
    class.term <- list(method=method, bounds=cutoff, categ = c(0,2), delta=delta, alpha=0.1, beta=0.1, indeterminate=TRUE)
  }
  
  cat_terminate2<- list(term=c("fixed", "class"), score="WLE", n.min=1, n.max=test.len.phase2, c.term =class.term)
  
  
  #*******************************
  #      Can we classify some of the simulees after phase1?
  #*******************************
  
  # If a simulee's cat_theta (after phase1) is equal or larger than a est.theta1.no.class value (for example cat_theta >= 1.2), do not classify the simulee.
  for(i in 1:N){
    if(catIrt.object$cat_theta[i] >= est.theta1.no.class){
      cat_indiv2_phase1[[i]]$cat_categ <- NA
    } else {
      # termGLR/termSPRT/termCI is a built-in functionof catIrt. However, brm or grm version of logLik function was automatically selected. E.g. termGLR is used during catIrt simulations, termGLR2 was duplicated with manual change in logLik.grm
      if(method=="SPRT"){
        cat_indiv2_phase1[[i]]$cat_categ <- termSPRT(cat_par = cat_indiv2_phase1[[i]]$cat_params, 
                                                     cat_theta = init.theta2.vec[i], 
                                                     cat_resp = cat_indiv2_phase1[[i]]$cat_resp2, 
                                                     catMiddle = cat_middle2, 
                                                     catTerm = cat_terminate2)
      }else if(method=="GLR"){
        cat_indiv2_phase1[[i]]$cat_categ <- termGLR(cat_par = cat_indiv2_phase1[[i]]$cat_params, 
                                                    cat_theta = init.theta2.vec[i], 
                                                    cat_resp = cat_indiv2_phase1[[i]]$cat_resp2, 
                                                    catMiddle = cat_middle2, 
                                                    catTerm = cat_terminate2)
      } else if(method=="CI"){
        cat_indiv2_phase1[[i]]$cat_categ <- termCI(cat_par = cat_indiv2_phase1[[i]]$cat_params,
                                                   cat_theta=init.theta2.vec[i], 
                                                   cat_resp = cat_indiv2_phase1[[i]]$cat_resp2,
                                                   cat_sem=cat_indiv2_phase1[[i]]$cat_sem, 
                                                   catTerm=cat_terminate2)
      }
    }
  }
  
  
  # Add true category after phase1
  for( i in 1:N){
    cat_indiv2_phase1[[i]]$true_categ <- class.term$categ[sum( cat_indiv2_phase1[[i]]$true_theta2 > class.term$bounds ) + 1]
  }

  
  #*******************************
  #       Start Phase 2 
  #*******************************
  # Specify person.vec. This is for deciding which simulee will take phase2 (1), which will not (0).
  # Who will not take phase2:
  # a) If est.theta1 is bigger than est.theta1.no.class (these simulees are not classified at all)
  # b) If a simulee already reached to mocca test length during phase1
  # c) If a simulee already classified on dimension2 after phase1
  
  person.vec <- vector(mode = "numeric", length = N)
  for(i in 1:N){
    if(catIrt.object$cat_theta[i] >= est.theta1.no.class | catIrt.object$cat_length[i] == mocca.n.max | cat_indiv2_phase1[[i]]$cat_categ %in% class.term$categ ){
      person.vec[i] <- 0
    } else {
      person.vec[i] <- 1
    }
  }
  
  # Run catIrt for dimension 2
  dim2.object <- catIrt2(theta = theta2, params=params2, resp=NULL, it=it.vec, person.vec=person.vec, mod="grm", catStart = cat_start2, catMiddle = cat_middle2, catTerm=cat_terminate2, progress=T, ddist=NULL)
  
  
  #*******************************
  #       Combine initial (after phase1) and final classification (after phase2)
  #*******************************
  
  # If a categorization is estimated after phase1 for a simulee or if test length is reached to mocca test length in phase1, use estimated category after phase1.
  for (i in 1:N){
    if(catIrt.object$cat_theta[i] >= est.theta1.no.class | catIrt.object$cat_length[i] == mocca.n.max | cat_indiv2_phase1[[i]]$cat_categ %in% class.term$categ ){
      # Use dimension 2 after phase1 outputs as mocca outputs
      dim2.object$mocca_theta2[i] <- cat_indiv2_phase1[[i]]$cat_theta
      dim2.object$mocca_categ[i] <- cat_indiv2_phase1[[i]]$cat_categ
      dim2.object$mocca_length[i] <- cat_indiv2_phase1[[i]]$cat_length <- catIrt.object$cat_length[i]
    } else {
      # Use dimension 2 after phase2 outputs as mocca outputs  
      dim2.object$mocca_theta2[i] <- dim2.object$cat_theta[i]
      dim2.object$mocca_categ[i] <- dim2.object$cat_categ[i]
      # mocca test length as phase1 + phase2 test length
      dim2.object$mocca_length[i] <- catIrt.object$cat_length[i] + dim2.object$cat_length[i]
    }
  }
  
  # Add cat_indiv2_phase1 object to dim2 object
  dim2.object$cat_indiv2_phase1 <- cat_indiv2_phase1
  
  # Add dimension 1 outputs to a final mocca object
  mocca <- list(dim1 = catIrt.object, dim2 = dim2.object)
  
  # Return
  return(mocca)

}

