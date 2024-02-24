generate_z_score = function(x) {
	rm = rowMeans(x, na.rm=na.rm)
	sd = apply(x, 1, sd, na.rm = na.rm)
	sweep(x, 1, rm) -> x
	x <- sweep(x, 1, sd, "/")
}
