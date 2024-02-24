dt.to.df = function(xx){
	# convert data.table to dataframe
	if(all(class(xx) == 'data.frame')) { return(xx) }
	# xx = dsn
	cn = xx[, 1, with=F]; cn
	cn = as.character(unlist(cn))
	ret = data.frame(row.names=cn)
	ret = cbind(ret, xx[, 2:ncol(xx)])
	return(ret)
}
