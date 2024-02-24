df.to.dt = function(idf, rn = NULL){
	require(data.table)
	if(any(class(idf) == 'data.table')){ return(idf) }

	idt = as.data.table(idf, keep.rownames=T)
	if(!is.null(rn)){
		setnames(idt, 'rn', rn)
	}
	return(idt)
}
