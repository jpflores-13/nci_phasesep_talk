## https://github.com/ManjariKiran/DiffLTS/blob/main/MH_innerg.R
## Grabs pixels greater than a specific manhattan distance
mh_index <- function(buffer, loop, inner){  #buffer: pixels sourrounding the matrix, loop: loop, inner: MH distance from central pixel
  m=(buffer*2)+1 #Set the matrix size
  center <- buffer+1 #Set the senter pixel
  M <- matrix(data=NA,nrow=m,ncol=m) #Create a matrix with m column and row filled with NA
  for(j in 1:m){                     #For loop to fill numbers as MH distance
    l=buffer+j
    for(i in 1:m){
      k=(m+1)-j
      if((i <= (buffer+1)) && (j <= (buffer+1))){
        M[i,j]<-k-i
      }
      if((i <= (buffer+1)) && (j > (buffer+1))){
        M[i,j]<-j-i
      }
      if((i > (buffer+1)) && (j <= (buffer+1))){
        M[i,j]<-i-j
      }
      if((i > (buffer+1)) && (j > (buffer+1))){
        M[i,j]<-i-k
      }
    }
  }
  inner_val <- which(M >= inner) #Filter the cells with values greater than inner 
  new <- loop[inner_val]  #Extract the signals for all the pixel outside inner
  return(new)   #return the vector with signal values
}
