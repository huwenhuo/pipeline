fsave = function(iobject, iname = deparse(subsititute(iobject)), ipath=NULL){
	if(is.null(ipath)){
		ipath = getwd()
	}
	gsub("-", '_', iname) -> iname; 
	gsub("\\.", '_', iname) -> iname
	save(iobject, file= paste0(ipath, '/', iname, '.RData'))
}
