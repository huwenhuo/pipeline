finder = function(loc='', pattern='*', ttype='f'){
	loc = 'res/gsea'
	pattern=c('*h.all*', 'lps')
	if(loc == ''){loc = getwd()}
	cmd = paste0('find ', loc, ' -type ', ttype, ' -iname "', pattern[1], '"'); 
	if(length(pattern) > 1){
		for(i in 2:length(pattern)){
			cmd = paste0(cmd, ' | grep -i ', pattern[i])
		}
	}
	#cat(cmd, '\n')
	system(cmd, intern=T)
}
