fetch.file.fun = function(fname, project = 'micropapillary', 
			  remote.dir = '/juno/work/solitlab/huw/solit/study/hiseq/', 
			  lcl.dir = NULL, hostname = 'xbio.mskcc.org'){

	if(is.null(lcl.dir)){
		lcl.dir = dirname(fname)
		if(!dir.exists(lcl.dir)){
			stop('local directory not existed!')
		}
	}

	fname = basename(fname)
	
	remote.fname = paste0(remote.dir, '/', project, '/', fname)
	cmd = paste0('ssh -p 2222 -A huw@', hostname, ' scp selene:', remote.fname, ' .')
	message(cmd)
	system(cmd)

	cmd = paste0('scp -P 2222 huw@', hostname, ':~/', fname, ' ', lcl.dir, '/')
	message(cmd)
	system(cmd)

	cmd = paste0('ssh -p 2222 -A huw@', hostname, ' rm ', fname)
	message(cmd)
	system(cmd)
}
