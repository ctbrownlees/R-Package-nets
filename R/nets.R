
.packageName <- "nets"

.First.lib <- function(lib, pkg)
{
     library.dynam("nets", pkg, lib)

}

nets <- function(){
     a <- c(1,2,3)
     results <- .C("adding", as.double(a), as.integer(length(a)), ab = double(length(a)))
     print( results )
}
