gg.stat = function(gg){
	gg = gg + stat_summary(fun.data = "mean_cl_boot", colour = "red", size = .2) 
	gg
}
