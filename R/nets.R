
nets <- function(y,lambda="auto",q="auto"){

	# checks
	is.matrix(y)

	# problem dimensions
	T <- nrow(y)
	n <- ncol(y)

	# parsing options
	if( !is.character(lambda) ){ 
		lambda.selection <- "fixed"
	}
	else{
		lambda.selection <- lambda
	}
	


}

