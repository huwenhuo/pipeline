#bamcount.jobname='a'
#peakfile.cmp.sorted='aa'
#bamfile='cc'
#peak.coverage.file = 'd'
#param.list = list(jobname = bamcount.jobname, mem=10, cpu=2, We='10:10', cmd=paste0(cfg$bedtools2, '/coverageBed'), a=peakfile.cmp.sorted, b=bamfile, paste0(' -sorted -counts > ', peak.coverage.file));
#param.list
bsub.invoke = function(param.list, cfg, ifrun=F){
	param.list
	names(param.list)
	noname = param.list[names(param.list)=='']
	noname = paste(noname, collapse=' ')
	noname
	param.list = param.list[names(param.list) != '']
	param.list
	## bsub parameter
	bsub.p = c('cpu', 'mem', 'postdone', 'jobname', 'We', 'statsDir', 'cwd')
	bsub.param = param.list[bsub.p]
	names(bsub.param)[is.na(names(bsub.param))] = setdiff(bsub.p, names(bsub.param))
	## cmd
	cmd = param.list$cmd
	## all left 
	param.list = param.list[setdiff(names(param.list), c('cmd', bsub.p))]
	param.list
	## make the head of bsub
	if(is.null(bsub.param$We)){ bsub.param$We = '10:10' }
	if(is.null(bsub.param$cwd)){ bsub.param$cwd = getwd() }
	if(is.null(bsub.param$statsDir)){  
		if(dir.exists(cfg$statsDir)){
			bsub.param$statsDir = cfg$statsDir
		}else{
			bsub.param$statsDir = '.'
		}
	}
	bsub.hd = paste0(cfg$BSUB, " -J ", bsub.param$jobname, " -e ", bsub.param$statsDir, '/', bsub.param$jobname, ".err -o ", bsub.param$statsDir, '/', bsub.param$jobname, ".std ")
	bsub.hd = paste0(bsub.hd, " -cwd ", bsub.param$cwd, ' -We ', bsub.param$We, ' -R "rusage[mem=', bsub.param$mem, ']" -R "rusage[iounits=0]" -n ', bsub.param$cpu, ' ')
	bsub.hd
	bsub.param$postdone
	if(!is.null(bsub.param$postdone)){
		postdone = unlist(strsplit(bsub.param$postdone, split=" "))
		postdone = paste(postdone, collapse='|')
		runningjobs = system('bjobs -w', intern=T)
		runningjobs = runningjobs[2:length(runningjobs)]
		nr = length(runningjobs)
		as.data.table(matrix(unlist(lapply(runningjobs, strsplit, " +")), byrow=T, nrow=nr)) -> runningjobs
		runningjobs = runningjobs$V7
		intersect(runningjobs, postdone) -> ip
		if(length(ip) > 0){
			for( i in 1:length(ip)){
				bsub.hd = paste0(bsub.hd, ' -w "post_done(', ip[i], ')" ')
			}
		}
	}
	## make the right side of the bsub
	bar = c('-', '--')
	bar = bar[as.numeric(nchar(names(param.list)) > 1) + 1]
	param = paste0(bar, names(param.list), ' ',  param.list)
	param = paste(param, collapse=' ')
	param = paste(param, noname)
	param
	## put them together
	cmd.full = paste0(bsub.hd, ' "', cmd, ' ', param, ' "')
	cmd.full
	if(ifrun == T){
		system(cmd.full)
	}
	cmd.full
}
