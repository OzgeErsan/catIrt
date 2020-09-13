
######################################
# Create a wrapper function for the entire MOCCA process 
# max test length for MOCCA is 40 here 
######################################

mocca <- function(cutoff, method, delta, at.select, 
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
 
 # ****************************
 # Code below will be deleted, we included it.vec inside mocca.transit
 # ******************************************
  
  # Get polytomous item responses on phase1 items
#  resp <- matrix(0, nrow = N, ncol = J)
#  it.vec <- matrix(0, nrow = N, ncol = J)
  
#  # Simulate type of incorrect response.Instead of 0/1, they are now 1,2,3
#  for(i in 1:N){
#    for(j in 1:J){
#      if(catIrt.object$full_resp[i,j]==1){
#        resp[i, j]<-2
#      } else {
#        params2.vec <- params2[j, ]
#        p1 <- exp(params2.vec[1]*(theta2[i]-params2.vec[2]))/(1+exp(params2.vec[1]*(theta2[i]-params2.vec[2])))
#        p0 <- 1-p1
#        p2 <- exp(params2.vec[1]*(theta2[i]-params2.vec[3]))/(1+exp(params2.vec[1]*(theta2[i]-params2.vec[3])))
#        p2.star <- p2/(p0+p2)
#        u <- runif(n=1,min=0,max=1)
#        if(p2.star > u){resp[i, j] <- 3} else {resp[i, j] <- 1}
#      }
#    }
#    # Fill it.vec to block the items already administered in phase1 during phase2
#    for(j in 1:length(catIrt.object$cat_indiv[[i]]$cat_it)){
#      it.vec[i, catIrt.object$cat_indiv[[i]]$cat_it[j]] <- 1
#    } 
#  }
  
  # Run mocca.transit function
  mocca.transit <- mocca.transit(theta2=theta2, params2=params2, catIrt.object=catIrt.object)
  
  # Get phase2 graded item responses and related info
  cat_indiv2_phase1 <- mocca.transit$cat_indiv2
  
  # Get it.vec that indicates which item will be blocked during phase2
  it.vec <- mocca.transit$it.vec
  
  # Add true theta2 to cat_indiv2_phase1
  cat_indiv2_phase1$true_theta2 <- theta2
  
  # Add phase1 test length to cat_indiv2_phase1
  cat_indiv2_phase1$cat_length <- catIrt.object$cat_length
  
  # Estimate theta2 for initial theta values for phase2
  init.theta2.vec <- vector(mode = "numeric", length = N)
  
  # Use catIrt built-in function to estimate theta for dimension2 after phase1 &
  # Add phase1 wleEst's sem output to cat_indiv2_phase1. This will be used during termCI after phase1
  wleEst.list <- vector(mode="list", length = N)
  for(i in 1:N){
    wleEst.list[[i]] <- wleEst(cat_indiv2_phase1[[i]]$cat_resp2, params=params2[cat_indiv2_phase1[[i]]$cat_it,], range = c(-4.5, 4.5), mod="grm")
    init.theta2.vec[i]<-cat_indiv2_phase1[[i]]$init.theta2 <- wleEst.list[[i]]$theta
    cat_indiv2_phase1[[i]]$cat_sem <-  wleEst.list[[i]]$sem
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
  
  # Create a vector to fill with classification decision with phase1 item responses
  phase1.class.vec <- vector(mode = "numeric", length = N)
  
  # If a simulee's cat_theta (after phase1) is equal or larger than a est.theta1.no.class value (for example cat_theta >= 1.2), do not classify the simulee.
  for(i in 1:N){
    if(cat_theta[i] >= est.theta1.no.class){
      phase1.class.vec[i] <- NA
    } else {
      # termGLR/termSPRT/termCI is a built-in functionof catIrt. However, brm or grm version of logLik function was automatically selected. E.g. termGLR is used during catIrt simulations, termGLR2 was duplicated with manual change in logLik.grm
      if(method=="SPRT"){
        for(i in 1:N){
          phase1.class.vec[i] <- cat_indiv2_phase1[[i]]$cat_categ <- termSPRT(cat_par = params2[cat_indiv2_phase1[[i]]$cat_it,], 
                                                                              cat_theta = init.theta2.vec[i], 
                                                                              cat_resp = cat_indiv2_phase1[[i]]$cat_resp2, 
                                                                              catMiddle = cat_middle2, 
                                                                              catTerm = cat_terminate2)}
      }else if(method=="GLR"){
        for(i in 1:N){
          phase1.class.vec[i] <- cat_indiv2_phase1[[i]]$cat_categ <- termGLR(cat_par = params2[cat_indiv2_phase1[[i]]$cat_it,], 
                                                                             cat_theta = init.theta2.vec[i], 
                                                                             cat_resp = cat_indiv2_phase1[[i]]$cat_resp2, 
                                                                             catMiddle = cat_middle2, 
                                                                             catTerm = cat_terminate2)}
      } else if(method=="CI"){
        for(i in 1:N){
          phase1.class.vec[i] <- cat_indiv2_phase1[[i]]$cat_categ <- termCI(cat_par=params2[cat_indiv2_phase1[[i]]$cat_it,],
                                                                            cat_theta=init.theta2.vec[i], 
                                                                            cat_resp = cat_indiv2_phase1[[i]]$cat_resp2,
                                                                            cat_sem=cat_indiv2_phase1[[i]]$cat_sem, 
                                                                            catTerm=cat_terminate2)}
      }
    }
  }
  
  
  # Finally attached current classification decisions to phase1 (catIrt.object) output
  cat_indiv2_phase1$cat_categ <- phase1.class.vec
  
  #*******************************
  #       Start Phase 2 
  #*******************************
  # Rest of the code includes temporary solutions. They will be changed.
  #********************************************
  
  # Run phase2
  #phase2 <- catIrt(theta = theta2, params=params2, resp=NULL, it=it.vec, mod="grm", catStart = cat_start2, catMiddle = cat_middle2, catTerm=cat_terminate2, progress=T, ddist=NULL)
  
  ## If max test length of mocca is also 40:
  # respecify cat_terminate2, since when n.max=40, current code gives error
  # It is okey for now as temporary solution, since results for these simulees will be replaced by phase1 results
  cat_terminate2_2 <- cat_terminate2
  for(i in 1:N){
    if(cat_terminate2_2$n.max[i]==0){
      cat_terminate2_2$n.max[i]=1
    }
  }
  
  phase2 <- catIrt(theta = theta2, params=params2, resp=NULL, it=it.vec, mod="grm", catStart = cat_start2, catMiddle = cat_middle2, catTerm=cat_terminate2_2, progress=T, ddist=NULL)
  
  
  #*******************************
  #       Combine initial (after phase1) and final classification (after phase2)
  #*******************************
  
  # If a categorization is estimated after phase1 for a simulee or if test length is reached to 40 in phase1, use that category.
  for (i in 1:N){
    if(cat_indiv2_phase1$cat_categ[i] != "ID" | cat_indiv2_phase1$cat_length[i] == mocca.n.max){
      #Use phase1 categ estimate and test length as mocca length
      phase2$mocca_categ[i] <- cat_indiv2_phase1$cat_categ[i]
      phase2$mocca_length[i] <- cat_indiv2_phase1$cat_length[i]
    } else {
      #Use phase2 categ estimate and compute mocca test length as phase1 + phase2 test length
      phase2$mocca_categ[i] <- phase2$cat_categ[i]
      phase2$mocca_length[i] <- catIrt.object$cat_length[i] + phase2$cat_length[i]
    }
  }
  
  # Add cat_indiv2_phase1 object to phase2 object
  phase2$cat_indiv2_phase1 <- cat_indiv2_phase1
  
  # Return
  return(phase2)
}

