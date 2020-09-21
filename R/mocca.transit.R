#############################
# Generate item responses for phase1 for "all" response options
#############################

# n.options can be 3 or 5. Currently this function works only for 3-alternative items.
mocca.transit <- function(theta2, params2, catIrt.object){
  
  # N and J are specified by the original function. 
  N <- dim(catIrt.object$full_resp)[1] # N is sample size (person)
  J <- dim(catIrt.object$full_resp)[2] # J is item bank size
  
  # Create a blank list to fill with selected information. cat_indiv2 is for dimension 2.
  cat_indiv2 <- vector(mode="list", length=N)
  cat_indiv <- catIrt.object$cat_indiv
  
  # Create a blank vector that will indicate which items were administered in phase1. So block them for phase2.
  it.vec <- matrix(0, nrow = N, ncol = J)
  
  # Run over each person
  for(i in 1:N){
    # cat_indiv[[i]][c(2,5)] passes item IDs, and examinee dimension 1 responses (0 or 1) to cat_indiv2
    cat_indiv2[[i]] <- cat_indiv[[i]][c(2,5)]   
    # cat_resp2 is for dimension 2 responses (e.g. 0,1,2)
    cat_indiv2[[i]]$cat_resp2 <- vector("numeric", length=length(cat_indiv[[i]]$cat_resp))
    # specify the class of cat_resp2 as grm (graded response model)
    class(cat_indiv2[[i]]$cat_resp2) <- "grm"
    
    # Insert simulated responses for dimension2 in "cat_resp2" (run over each response)
    for(n in 1:length(cat_indiv2[[i]]$cat_resp2)){
      # If cat_indiv[[i]]$cat_resp=1, cat_indiv2[[i]]$cat_resp2=2
      # If cat_indiv[[i]]$cat_resp=0, for cat_indiv2[[i]]$cat_resp2 generate responses (1 or 3). 
      # (Instead of {0,1,2}, {1,2,3} were used, since catIrt uses 1,2,3... as grm responses).
      if(cat_indiv[[i]]$cat_resp[n]==1){
        cat_indiv2[[i]]$cat_resp2[n]<-2} else {
          # params2 is item parameters of an item for phase2
          params2.vec <- params2[cat_indiv2[[i]]$cat_it[n], ]
          # Compute category probabilities (p0, p1, p2)
          # params2.vec first element is a-parameters; second is b1-parameters; third is b2-parameters
          p1.temp <- exp(params2.vec[1]*(theta2[i]-params2.vec[2]))/(1+exp(params2.vec[1]*(theta2[i]-params2.vec[2])))
          p2 <- exp(params2.vec[1]*(theta2[i]-params2.vec[3]))/(1+exp(params2.vec[1]*(theta2[i]-params2.vec[3])))
          p1 <- p1.temp-p2
          p0 <- 1-p1.temp
          # Compute relative probability of category 2
          p2.star <- p2/(p0+p2)
          u <- runif(n=1,min=0,max=1)
          if(p2.star > u){cat_indiv2[[i]]$cat_resp2[n] <- 3} else {cat_indiv2[[i]]$cat_resp2[n] <- 1}
        }
    }
    
    # Indicate which items were administered already by "1".
    for(j in 1:length(catIrt.object$cat_indiv[[i]]$cat_it)){
      it.vec[i, catIrt.object$cat_indiv[[i]]$cat_it[j]] <- 1
    }
  }
  return(list(cat_indiv=cat_indiv2, it.vec=it.vec))
}
  
  



