library(gurobi)
library(data.table)
library(dplyr)

multiway_partitions <- function(I, K) {
  sizes <- data.table(i=1:length(I),size=I)
  combos <- data.table(crossing(sizes,k=1:K))
  
  model = list()
  A <- NULL;  b <- NULL;
  sense <- NULL;
  current.row=1;
  
  #Each item can only be in one partition
  rows <- combos[,i]
  columns <- combos[,(i-1)*max(k)+k]
  entries <- rep(1,times=combos[,.N])
  new.for.b <- rep(1,times=length(I))
  new.for.sense <- rep("=",times=length(I))
  
  b <- c(b, c(new.for.b));
  sense <- c(sense,new.for.sense)
  new.for.A = data.frame(cbind(rows, columns, entries));
  A <- rbindlist(list(A,new.for.A), use.names=FALSE, fill=FALSE, idcol=NULL)
  new.for.A <- c(); new.for.b <- c(); new.for.sense <- c();
  rows <- c(); columns <- c(); entries <- c();
  current.row <- max(A$rows)+1
  
  #Define size of each partition
  rows <- c(current.row-1 + combos[,k],
            current.row-1 + 1:K)
  columns <- c(combos[,(i-1)*max(k)+k],
               combos[,.N] + 1:K)
  entries <- c(combos[,size],
               rep(-1,times=K))
  new.for.b <- rep(0,times=K)
  new.for.sense <- rep("=",times=K)
  
  b <- c(b, c(new.for.b));
  sense <- c(sense,new.for.sense)
  new.for.A = data.frame(cbind(rows, columns, entries));
  A <- rbindlist(list(A,new.for.A), use.names=FALSE, fill=FALSE, idcol=NULL)
  new.for.A <- c(); new.for.b <- c(); new.for.sense <- c();
  rows <- c(); columns <- c(); entries <- c();
  current.row <- max(A$rows)+1
  
  #Break symmetry and define min/max sized partitions
  rows <- c(current.row-1 + 1:(K-1),
            current.row-1 + 1:(K-1))
  columns <- c(combos[,.N] + 1:(K-1),
               combos[,.N] + 2:K)
  entries <- c(rep(-1,times=K-1),
               rep(1,times=K-1))
  new.for.b <- rep(0,times=K-1)
  new.for.sense <- rep(">=",times=K-1)
  b <- c(b, c(new.for.b));
  sense <- c(sense,new.for.sense)
  new.for.A = data.frame(cbind(rows, columns, entries));
  A <- rbindlist(list(A,new.for.A), use.names=FALSE, fill=FALSE, idcol=NULL)
  new.for.A <- c(); new.for.b <- c(); new.for.sense <- c();
  rows <- c(); columns <- c(); entries <- c();
  current.row <- max(A$rows)+1
  
  model$vtype <- c(
    rep('B',times=combos[,.N]),
    rep('C',times=K))
  
  model$lb <- c(
    rep(0,times=combos[,.N]),
    rep((sizes[,sum(size)]/K)*0,times=K))
  model$ub <- c(
    rep(1,times=combos[,.N]),
    rep((sizes[,sum(size)]/K)*5,times=K))
  
  model$obj <- c(
    rep(0,times=combos[,.N]),
    -1,
    rep(0,times=K-2),
    1
  )
  
  model$A          <- sparseMatrix(i=A$rows,j=A$columns,x=A$entries)
  model$modelsense <- 'min'
  model$rhs        <- b
  model$sense      <- sense
  rm(A); rm(b); rm(sense); 
  
  gc();
  solution <- gurobi(model,params=list(
    TimeLimit=2700,
    Heuristics=0.5,
    Threads=6,
    PreDepRow=1,
    Presolve=2,
    PreDual=2,
    PrePasses=20,
    PreSparsify=1,
    Cuts=2,
    CutPasses=500,
    CutAggPasses=1
  ))
  if(length(solution$x) == length(model$obj)){
    return(which(solution$x[1:combos[,.N]] ==1)%%K + 1)
  }
}
