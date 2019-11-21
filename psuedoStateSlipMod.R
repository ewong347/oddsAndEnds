

slip <- F
norm <- T

s <- "SAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAE"
char <- "S"
position <- 1

while (char!="E") {
  char <- substr(s, position,1)
  i <- samp(1 in 100)
  
  #Slip event
  if (norm&i==1) {
    slip <- T
    norm <- F
    #slip()
    
  }else {
    
    if (slip&i<50) {
      slip <- T
      norm <- F
      #slip()
    } else ()
    
  }
  
  #To Be completed...
  
  
}