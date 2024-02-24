options(width = 200)
library(data.table)

pid = '6048'
system(paste0('find /ifs/projects/CMO/ -name "*_sample_patient.txt" > ~/pid.list'))
pid.all = scan('~/pid.list', character())
head(pid.all)
pid.sel = pid.all[grep(pid, pid.all)]
pid.sel.dt = lapply(pid.sel, fread)
pid.sel.dt[[1]]
rbindlist(pid.sel.dt) -> pid.sel.samples
pid.sel.samples[grep('JuB3_P10_', Collab_ID), ]

colid = c(
	  "SuB27_T", 
	  "SuB27_Org_P0", 
	  "SuB27_Org_P7", 
	  "JuB3_slides", 
	  "JuB3_P1_Pellet", 
	  "JuB3_P10_Pellet", 
	  "SuB19_T", 
	  "SuB19_Org_P0", 
	  "SuB19_Org_P9") 
colid = paste0('s_', colid)
colid

fwrite(pid.sel.samples[Collab_ID %in% colid, ], file='~/pid.xls', sep="\t", quote=F)
system('rsync -avur ~/pid.xls mski1925:~')


