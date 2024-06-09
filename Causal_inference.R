rm(list=ls())

library(pcalg)
library(igraph)
db <- read.csv("~/results.csv")  

############ load variables
k <- as.numeric(db$Reg_BS_BM_coef)
tsa <- as.numeric(db$SumTSAfull)
hd <- as.numeric(db$HumanPopDensity)
pro <- as.numeric(db$IUCN_Ia)
coral <- as.numeric(db$Buffer10Km)

########## transform continuous to binary variables using median

bk <- k
bk[k > median(k,na.rm=T)] <- 1
bk[k <= median(k,na.rm=T)] <- 0

bhd <- hd
bhd[hd > median(hd,na.rm=T)] <- 1
bhd[hd <= median(hd,na.rm=T)] <- 0

btsa <- tsa
btsa[tsa > median(tsa,na.rm=T)] <- 1
btsa[tsa <= median(tsa,na.rm=T)] <- 0

bcoral <- coral
bcoral[coral >median(coral,na.rm=T)] <- 1
bcoral[coral <=median(coral,na.rm=T)] <- 0

bpro <- pro

dat0 <- as.matrix(data.frame(bhd,bk,bpro,bcoral,btsa))
dat <- na.omit(dat0)
head(dat)

### Causal discovery with IC algorithm
IC <- function(V, alpha=0.05){
  
  l <- ncol(V)
  co <- combn(l, 2)
  E <- v_struc <- list()
  
  # Step 1: finding the skeleton
  for(i in 1:ncol(co)){
    vars <- co[,i]
    a <- vars[1]
    b <- vars[2]
    c <- seq(l)[-c(a,b)]
    
    if(length(c)==1){ 
      S <- c
    } else{
      S <- do.call("c", lapply(seq(length(c)), function(i) combn(c, i, FUN = list)))
    }
    
    t1 <- gSquareBin(a,b,c(),V) 
    t2 <- lapply(S, function(x) gSquareBin(a,b,c(x),V) )
    
    E[[i]] <- ifelse( all(unlist(t2) < alpha) & t1 < alpha, 1, 0)
    
    if(t1 > alpha){
      t3 <- lapply(c, function(x) gSquareBin(a,b,c(x),V) )
      if(any(unlist(t3) < alpha)) v_struc[[i]] <- paste0(a, "-->", c[which(unlist(t3) < alpha)], "<--", b)
    } 
    
  }
  
  Edges <- co[, which( unlist(E)==1 )]
  Sk <- lapply(seq(ncol(Edges)), function(y) paste0(Edges[1,y], "--", Edges[2,y]  ))
  Sk <- unlist(Sk)
  
  # Step 2: finding colliders
  
  v_struc <- unlist(v_struc)
  if( length(v_struc) > 0 ){
    check_neighbor_V <- lapply(seq(length(v_struc)), function(z){ 
      vec <- as.numeric(strsplit(v_struc[z], "-->|<--" )[[1]])
      a <- vec[1]
      c <- vec[2]
      b <- vec[3]
      neigh1 <- lapply(seq(ncol(Edges)), function(u) ifelse( all(c(a,c) %in% Edges[,u]),1,0))
      neigh2 <- lapply(seq(ncol(Edges)), function(u) ifelse( all(c(b,c) %in% Edges[,u]),1,0))
      N <- ifelse( sum(unlist(neigh1))==1 & sum(unlist(neigh2))==1, 1, 0)
    }) 
    V_struc <- v_struc[which(check_neighbor_V==1)]
  } else{
    V_struc <- "none"
  } 
  
  
  # Step 3: allocating directions to remaining undirected edges
  
  if(  !(length(V_struc)==0) &  all(!(V_struc == "none")) ){
    v0 <- lapply(seq(length(V_struc)), function(g){
      n <- as.numeric(strsplit(V_struc[g], "-->|<--" )[[1]]) 
      m <- matrix( c(n[1], n[2], n[3], n[2]), ncol=2)  })
    mat_E <- unique(do.call("cbind", v0), MARGIN=2)
    rem_E <- lapply(seq(ncol(Edges)), function(u) 
      ifelse( any(unlist(lapply(seq(ncol(mat_E)), function(k) all(mat_E[,k] %in% Edges[,u]) ))) ,0,1))
    remaining_edges <- matrix(Edges[,which(unlist(rem_E)==1)], nrow=2)
    
    n <- 0
    new_edge <- rep(NA, ncol(remaining_edges) )
    while( n < ncol(remaining_edges) ){
      
      for(z in seq(ncol(remaining_edges)) ){
        #z <- 1  
        if( is.na(new_edge[[z]]) ){  
          vec <- remaining_edges[,z]
          b <- vec[1]
          c <- vec[2]
          
          # Rule 1. Orient b − c into b → c whenever there is a → b such that a and c are non-adjacents
          neigh_b <- lapply(seq(ncol(mat_E)), function(u){ 
            x <- match(b , mat_E[,u])
            adj <- lapply(seq(ncol(Edges)), function(m) ifelse( all(c(mat_E[1,u],c) %in% Edges[,m]),1,0))
            ifelse(is.na(x) | any(unlist(adj)==1), 0, 1)  
          })
          neigh_c <- lapply(seq(ncol(mat_E)), function(u){ 
            x <- match(c , mat_E[,u])
            adj <- lapply(seq(ncol(Edges)), function(m) ifelse( all(c(mat_E[1,u],b) %in% Edges[,m]),1,0))
            ifelse(is.na(x) | any(unlist(adj)==1), 0, 1)  
          })
          
          if( any(unlist(neigh_b)==1) ){
            new_edge[z] <- paste0(b, "-->", c)
            mat_E <- cbind(mat_E, c(b,c))
          } 
          
          if( any(unlist(neigh_c)==1) ){  
            new_edge[z] <- paste0(c, "-->", b)
            mat_E <- cbind(mat_E, c(c,b))
          } 
          
        }    
      }
      
      n <- length(na.omit(new_edge))
      
    }  
    
    ## Directed edges
    dir_E <- unique(mat_E, MARGIN=2)
    dir_edges <- c(paste0(t(dir_E)[,1], "-->", t(dir_E)[,2]))
  }else{
    dir_E <- NA
    dir_edges <- 0
  }    
  
  # plot graph
    g <- igraph::graph.data.frame(t(dir_E), directed=T  ) 
    V(g)$name <- colnames(V)[as.numeric(V(g)$name)]
  
  
  
  return(list(Skeleton=Sk, V_Struc=V_struc, Dir_Edges=dir_edges, G=g))
  
}

dg <- IC(dat)

dev.off()
plot(dg$G, layout=layout.fruchterman.reingold,
     vertex.size = 50,      
     vertex.label = V(dg$G)$name,
     vertex.label.cex = 1,
     vertex.label.color = "black",
     vertex.color="white" )

# Testing genuine effect between X=MPA and Y=structure (dep low indp high)

# i: X=Protection and Y=Structure are dependent
gSquareBin(3,2,c(1,4,5),dat)  

# iia: Z=Human and X=Protection are dependent
gSquareBin(3,1,c(5,4,2),dat) 

# iia.1: W=TSA and Z=Human are independent
gSquareBin(5,1,c(),dat) 
# iia.1: X=Protection and W=TSA are dependent
gSquareBin(3,5,c(),dat)  

# iia.2: Z=Human and W=Coral are independent
gSquareBin(1,4,c(),dat) 
# iia.2: X=Protection and W=Coral are dependent
gSquareBin(3,4,c(),dat) 

# iiia: Z=Human and Y=Structure are dependent
gSquareBin(1,2,c(4,5),dat) 
# iiia: Z=Human and Y=Structure are independent
gSquareBin(1,2,c(4,5,3),dat) 

# iib: Z=TSA and X=Protection are dependent
gSquareBin(5,3,c(4,1,2),dat) 

# iib.1: Z=TSA and W=Human are independent
gSquareBin(5,1,c(),dat) 
# iib.1: X=Protection and W=Human are dependent
gSquareBin(3,1,c(),dat) 

# iib.2: Z=TSA and W=Coral are independent
gSquareBin(5,4,c(),dat) 
# iib.2: X=Protection and W=Coral are dependent
gSquareBin(3,4,c(),dat) 

# iiib: Z=TSA and Y=Structure are dependent
gSquareBin(5,2,c(4,1),dat) 
# iiib: Z=TSA and Y=Structure are independent
gSquareBin(5,2,c(4,1,3),dat) 

# iic: Z=Coral and X=Protection are dependent
gSquareBin(4,3,c(5,1,2),dat) 

# iic.1: Z=Coral and W=Human are independent
gSquareBin(4,1,c(),dat) 
# iic.1: X=Protection and W=Human are dependent
gSquareBin(3,1,c(),dat) 

# iic.2: Z=Coral and W=TSA are independent
gSquareBin(4,5,c(),dat) 
# iic.2:X=Protection and W=Coral are dependent
gSquareBin(4,3,c(),dat) 

# iiic: Z=Coral and Y=Structure are dependent
gSquareBin(4,2,c(5,1),dat) 
# iiic: Z=Coral and Y=Structure are independent
gSquareBin(4,2,c(5,1,3),dat) 


# Direct effect of MPA on structure
d1 <- mean(bk[bpro==1],na.rm = T)
d2 <- mean(bk[bpro==0],na.rm = T)
if (is.na(d1)) {d1 <- 0}
if (is.na(d2)) {d2 <- 0}
d <- d1 - d2
print(d) 




