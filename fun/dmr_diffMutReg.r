# ## for DMP
# wes.maf
# sampleID # with Tumor_Sample_Barcode and sample.class which is either 'Primary', or 'Tumor'
# geneset # to plot
# plot.sample # order of samples to plot, could be direct from the oncoplot of all the mutations
# pid # to extract patient name from Tumor_Sample_Barcode 
# sid # to extract patient name from Tumor_Sample_Barcode 

dmr = function(maf.dt, sample.id, geneset, plot.sample, pid, sid, ifonly = F, fname = 'dmr_wes.pdf'){
	# ifonly discordance genes are plot

	pid = "bla-[0-9]+"
	sid = "[T|M][0-9]?$"

	plot.id = c(190, 188, 186, 208, 205, 204, 187, 200, 210, 202, 185, 194,203, 193, 195, 201, 211, 206, 191, 189, 192)
	plot.id = paste0(c('DS-bla-'), rep(plot.id, each=2), c('-T1', '-T2'))
	ind = c(4,20)
	mm = c('DS-bla-188-M', 'DS-bla-202-M1')
	plot.id = c(plot.id, mm)[order(c(seq_along(plot.id), ind-0.5))]
	plot.id
	plot.sample = paste0('bla-', plot.id)
	plot.sample

	source('~/pipeline/fun/param.r') # for other parameters needed here: color.table, simple.table, color.value et al

	sampleID = data.table(`Sample ID` = unique(maf.dt$Tumor_Sample_Barcode))
	sampleID[, `Sample Class` := 'Primary']
	sampleID[grep("T2", `Sample ID`), `Sample Class` := 'Tumor']
	sampleID[grep("M", `Sample ID`), `Sample Class` := 'Meta']
	sampleID[, patient := stringr::str_extract(`Sample ID`, pid)]
	sampleID[, sample := stringr::str_extract(`Sample ID`, sid)]
	sampleID[, Tumor_Sample_Barcode := `Sample ID`]
	sampleID[, sample.class := `Sample Class`]
	sample.id = sampleID
	sample.id

	gg.sel = c('TP53', 'KMT2D', 'PIK3CA','ATM','EP300','ARID1A','FGFR3','KMT2C','STAG2','RHOA','ERBB3','NFE2L2','CDKN1A','CREBBP','CUL1','ERBB2','KMD6A','NF1','TSC1','ACTB')

	x.dt = merge(maf.dt, sample.id[, .(Tumor_Sample_Barcode, sample.class)], all.x=T, by = 'Tumor_Sample_Barcode')
	x.dt = x.dt[Hugo_Symbol %in% gg.sel,  ][Variant_Classification %in% simple.table$Variant.Classification, ]
	x.dt[, sample.name := str_extract(Tumor_Sample_Barcode, sid), by=1:nrow(x.dt)]
	x.dt[, patient.name := str_extract(Tumor_Sample_Barcode, pid), by=1:nrow(x.dt)]
	x.dt[oncogenic == '', oncogenic := 'Unknown']
	x.dt = merge(x.dt, oncogenic.table, by = 'oncogenic', all.x=T)
	x.dt = merge(x.dt, simple.table, by.x = 'Variant_Classification', by.y = 'Variant.Classification', all.x=T)
	x.dt[, code := paste0(simple.value, ' ', oncogenic.value)] # for plot color
	x.dt[, Hugo_Symbol := factor(Hugo_Symbol, levels=rev(gg.sel))]
	x.dt[, patient.name := factor(patient.name, levels=plot.sample)]
	x.dt[, tag := paste0(patient.name, ':', Hugo_Symbol)]; # id for each paired square to plot
	x.dt[, genes := Hugo_Symbol] # an alias easier for typing
	x.dt

	# filter out the exact same mutations in maf
	x.dt[, tag2 := paste0(Tumor_Sample_Barcode, '_', Chromosome, '_', Start_Position, '_', Hugo_Symbol, '_', HGVSp_Short)] # 
	x.dt = x.dt[!duplicated(tag2), ]
	x.dt[, tag2 := NULL]
	x.dt
	
	x.dt[grep('KMT2D', Hugo_Symbol), ][, .(Tumor_Sample_Barcode, patient.name, Hugo_Symbol, HGVSp_Short, sample.class)]

	## make a table for cooridnations and colors(code1/2) of rect to plot
	dt.xy = data.table(gene = character(), patient.name=character(),
			   xleft1=integer(), xright1=integer(),
			   xleft2=integer(), xright2=integer(),
			   ybottom1=integer(), ytop1=integer(),
			   ybottom2=integer(), ytop2=integer(),
			   tag2= character(), code1 = character(), code2 = character(), concordance=integer())
	dt.xy

	ii
	for(ii in unique(x.dt$tag)){
		dt.1 = x.dt[tag == ii, ]; dt.1
		# tag2 for variants to plot in the rect
		dt.1[, tag2 := paste0(Chromosome, ':', Start_Position, ':', End_Position)]; 
		dt.1[, sample.class := factor(sample.class, levels=c('Primary', 'Tumor', 'Meta'))]
		# data table for number of rows in the big rect
		sq.dt = dcast(tag2 ~ sample.class, data=dt.1, value.var='code', drop=F, fun.aggregate = function(x){paste(unique(x), collapse=':')}); 
		sq.dt[is.na(Primary), Primary := adjustcolor('grey80')]
		sq.dt[is.na(Tumor),   Tumor   := adjustcolor('grey80')]; sq.dt
		sq.dt
		nn 	= nrow(sq.dt)
		istep 	= 0.8 / nn; istep
		yy.seq	= seq(from = 0, to = .8, length.out = nn + 1);  yy.seq
		yy.seq	= yy.seq[1:(length(yy.seq)-1)]; yy.seq
		xleft1 	= rep(as.integer(dt.1$patient.name[1]), nn); xleft1
		ybottom1= unique(as.integer(dt.1$genes)) + yy.seq; ybottom1
		xright1 = xleft1 + .35; xright1
		ytop1 	= ybottom1 + istep; ytop1
		xleft2 	= rep(as.integer(dt.1$patient.name[1]), nn)+.45; xleft2
		ybottom2= ybottom1; ybottom2
		xright2 = xleft2 + .35; xright2
		ytop2 	= ytop1;
		tag2 	= sq.dt$tag2
		code1	= sq.dt[, Primary]; code1
		code2	= sq.dt[, Tumor]; code2
		gene 	= dt.1[1, genes]
		patient.name	= dt.1[1, patient.name]
		dt.xy.1  =data.table(
				     xleft1 = xleft1, ybottom1 = ybottom1, xright1 = xright1, ytop1 = ytop1, 
				     xleft2 = xleft2, ybottom2 = ybottom2, xright2 = xright2, ytop2 = ytop2,
				     code1 = code1, code2 = code2, tag2 = tag2, patient.name = patient.name, 
				     gene = gene, concordance = '') 
		dt.xy.1[, concordance := ifelse(code1 == code2, 'yes', 'no')]
		dt.xy.1
		dt.xy = rbind(dt.xy, dt.xy.1)
	}
	dt.xy

	dt.xy[, color1 :=  color.table[dt.xy$code1, color.value]]
	dt.xy[, color2 :=  color.table[dt.xy$code2, color.value]]
	dt.xy[order(patient.name), ]

	n.genes = length(gg.sel); n.genes
	n.patients = length(levels(x.dt$patient.name)); 
	
	fname = 'dmg.pdf'; fname
	pdf(file=fname, width=7, height=6)
	plot('', type='n', xlim=c(0, n.patients + 4), ylim=c(0, n.genes + 2), xlab='', ylab='', xaxt = 'n', yaxt = 'n', bty = 'n', pch = '')
	# all the rect background color
	xleft = rep(1:max(as.integer(x.dt$patient.name)), times=n.genes);
	ybottom = rep(1:max(as.integer(x.dt$genes)), each=n.patients);
	xright = xleft +0.8;
	ytop = ybottom + .8;
	# gene names
	text(.5, 0.5+n.genes:1, label=gg.sel, xpd=T, adj=c(1, 1))
	text(1:n.patients, n.genes + 1, label=levels(x.dt$patient.name), srt=90, adj=c(0, 1), xpd=T)
	rect(xleft = xleft, ybottom = ybottom, xright = xright, ytop = ytop, col = 'grey70', border='NA')
	rect(xleft = dt.xy$xleft1, ybottom = dt.xy$ybottom1, xright = dt.xy$xright1, ytop = dt.xy$ytop1, col = dt.xy$color1, border='grey30')
	rect(xleft = dt.xy$xleft2, ybottom = dt.xy$ybottom2, xright = dt.xy$xright2, ytop = dt.xy$ytop2, col = dt.xy$color2, border='grey30')
	dev.off()

	scp(fname)

	tmp.x = copy(tmb.dt)
	tmp.x[, tag := paste0(patient, '_', sample.class)]
	tmp.x = tmp.x[order(tag), ]
	tmp.x[, perMb.m := mean(perMb), by=tag]
	tmp.x = tmp.x[!duplicated(tag), ]
	tmp.x [, patient := factor(patient, levels=plot.sample)]
	tmp.x = tmp.x[order(patient), ]
	tmp.x[, tsb := factor(tsb, levels=tsb)]

	gg = ggplot(tmp.x, aes(x=tsb, y=perMb.m)) + geom_bar(stat='identity') +
		theme(axis.text.x=element_text(angle=90))
	fname = 'tmb_wes.pdf'
	ggsave(gg, file=fname)
	scp(fname)

}

dmr.dis = function(maf.dt, sample.id, geneset, plot.sample, pid, sid, ifonly = T, fname = 'dmr_wes.pdf'){
	# ifonly discordance genes are plot

	source('~/pipeline/fun/param.r') # for other parameters needed here: color.table, simple.table, color.value et al

	x.dt = merge(maf.dt, sample.id[, .(Tumor_Sample_Barcode, sample.class)], all.x=T, by = 'Tumor_Sample_Barcode')
	x.dt[, sample.name := str_extract(Tumor_Sample_Barcode, sid), by=1:nrow(x.dt)]
	x.dt[, patient.name := str_extract(Tumor_Sample_Barcode, pid), by=1:nrow(x.dt)]
	x.dt[is.na(patient.name), .(Tumor_Sample_Barcode, patient.name)]
	#x.dt = x.dt[Hugo_Symbol %in% geneset,  ][Variant_Classification %in% simple.table$Variant.Classification, ]
	x.dt = x.dt[Variant_Classification %in% simple.table$Variant.Classification, ]
	x.dt[oncogenic == '', oncogenic := 'Unknown']
	x.dt = merge(x.dt, oncogenic.table, by = 'oncogenic', all.x=T)
	x.dt = merge(x.dt, simple.table, by.x = 'Variant_Classification', by.y = 'Variant.Classification', all.x=T)
	x.dt[, code := paste0(simple.value, ' ', oncogenic.value)] # for plot color
	x.dt[, Hugo_Symbol := factor(Hugo_Symbol)]
	#x.dt[, patient.name := factor(patient.name, levels=plot.sample)]
	x.dt[, tag := paste0(patient.name, ':', Hugo_Symbol)]; # id for each paired square to plot
	x.dt[, genes := Hugo_Symbol] # an alias easier for typing
	x.dt$code

	#x.dt.1 = x.dt[oncogenic.value != 'Unknown', ]
	x.dt.1 = x.dt[oncogenic %in% onco.sel, ]
	x.dt.1[, tag2 := paste0(patient.name, 'XX', Hugo_Symbol, 'XX',  Chromosome, 'XX', Start_Position, 'XX', End_Position, 'XX', code )]
	sq.dt = dcast(tag2 ~ sample.class, data=x.dt.1, value.var='tag2', drop=F, fun.aggregate = function(x){paste(unique(x), collapse=':')}); 
	sq.dt = sq.dt[Primary != Tumor, ]
	sq.dt[, patient := unlist(lapply(strsplit(tag2, 'XX'), '[[', 1))]
	sq.dt[, symbol := unlist(lapply(strsplit(tag2, 'XX'), '[[', 2))]
	sq.dt[, tag3 := paste0(patient, '_', symbol)]
	sq.dt

	plot.name = factor(unique(sq.dt$patient))
	plot.name
	geneset = factor(unique(sq.dt$symbol))
	geneset

	x.dt[, code := paste0(simple.value, ' ', oncogenic.value)] # for plot color
	x.dt[, tag2 := paste0(patient.name, 'XX', Hugo_Symbol, 'XX',  Chromosome, 'XX', Start_Position, 'XX', End_Position, 'XX', code )]
	x.dt[, tag3 := paste0(patient.name, '_', Hugo_Symbol)]
	x.dt.2 = x.dt[tag3 %in% unique(sq.dt$tag3), ]
	x.dt.2[, n.patient := factor(patient.name, levels=plot.name)]
	x.dt.2[, n.gene := factor(Hugo_Symbol, levels=geneset)]
	x.dt.2[, sample.class := factor(sample.class, levels=c('Primary', 'Tumor'))]
	x.dt.2[, tag3 := paste0(patient.name, '_', Hugo_Symbol)]
	x.dt.2[, tag3]
	x.dt.2[, nn := length(tag2), by = 'tag3']
	x.dt.2[, ni := 1:length(tag2), by = 'tag3']
	x.dt.2 = x.dt.2[order(tag3), ]
	x.dt.2[, .(n.gene, n.patient, nn, ni)]

	x.dt.2[, x.left := as.numeric(n.patient) + 0.4 * (as.numeric(sample.class) - 1)]
	x.dt.2[, x.right := x.left + .4]
	x.dt.2[, y.bot := as.numeric(n.gene) + .8*(ni-1)/nn]
	x.dt.2[, y.top := as.numeric(n.gene) + .8*ni/nn]
	x.dt.2[Hugo_Symbol == 'FGFR3', .(n.gene, n.patient, nn, ni, x.left, x.right, y.bot, y.top, color, HGVSp_Short)]
	x.dt.2[, color :=  color.table[code, color.value]]
	fname = 'discordance.pdf'
	pdf(file=fname, width=8, height=14)
	plot('', type='n', xlim=c(0, length(plot.name) + 4), ylim=c(0, length(geneset) + 2), xlab='', ylab='', xaxt = 'n', yaxt = 'n', bty = 'n', pch = '')
	text(.8, 0.5+1:length(geneset), label=geneset, xpd=T, adj=c(1, 1))
	text(as.numeric(plot.name)+.4, length(geneset) + 1, label=plot.name, srt=90, adj=c(0, .5), xpd=T)
	xleft = rep(1:max(as.integer(plot.name)), times=length(plot.name));
	ybottom = rep(1:max(as.integer(geneset)), each=length(plot.name));
	xright = xleft + 0.8;
	ytop = ybottom + 0.8;
	rect(xleft = xleft, ybottom = ybottom, xright = xright, ytop = ytop, col = 'grey80', border=NA)
	rect(xleft = x.dt.2$x.left, ybottom = x.dt.2$y.bot, xright = x.dt.2$x.right, ytop = x.dt.2$y.top, col = x.dt.2$color, border='white')
	#rect(xleft = x.dt.2$xleft2, ybottom = x.dt.2$ybottom2, xright = x.dt.2$xright2, ytop = x.dt.2$ytop2, col = x.dt.2$color2, border=NA)
	dev.off()
	scp(fname)

	color.table

	sq.dt[, p.t := paste0(Primary, ':', Tumor), by=tag]
	sq.dt[, p.t := paste0(unique(unlist(strsplit(p.t, ':'))), collapse=':'), by=tag]
	sq.dt = separate_rows(sq.dt, p.t, sep=':')
	sq.dt[, primary2 := p.t[grep(p.t, Primary)], by=1:nrow(sq.dt)]
	sq.dt[, tumor2 := p.t[grep(p.t, Tumor)], by=1:nrow(sq.dt)]
	sq.dt

	sq.dt[, patient.name := sub(':.*', '', tag)]
	sq.dt[, gene := sub('.*:', '', tag)]
	sq.dt[, mut.nn := .N, by = 'patient.name']
	sq.dt = sq.dt[order(mut.nn, decreasing=T), ]
	sq.dt[, patient.name := factor(patient.name, levels=unique(patient.name))]
	sq.dt[, gene.nn := .N, by = 'gene']
	sq.dt = sq.dt[order(gene.nn, decreasing=T), ]
	sq.dt[, gene := factor(gene, levels=unique(gene))]
	sq.dt

	sq.dt.1
	sq.dt.2[, xleft1 := as.integer(patient.name)]
	sq.dt.2[, xright1 := xleft1 + 0.35]
	sq.dt.2[, xleft2 := 0.5 + as.integer(patient.name)]
	sq.dt.2[, xright1 := xleft2 + 0.35]
	sq.dt.2[, ybottom1 := {
		yy.seq	= seq(from = 0, to = .8, length.out = .N + 1);  yy.seq
		yy.seq	= yy.seq[1:(length(yy.seq)-1)]; yy.seq
		yy.seq + as.integer(gene)
		yy.seq } , by = 'tag']
	sq.dt.2

	## make a table for cooridnations and colors(code1/2) of rect to plot
	dt.xy = data.table(gene = character(), patient.name=character(),
			   xleft1=integer(), xright1=integer(),
			   xleft2=integer(), xright2=integer(),
			   ybottom1=integer(), ytop1=integer(),
			   ybottom2=integer(), ytop2=integer(),
			   tag2= character(), code1 = character(), code2 = character(), concordance=integer())
	dt.xy
	for(ii in unique(x.dt$tag)){
		dt.1 = x.dt[tag == ii, ]; dt.1
		# tag2 for variants to plot in the rect
		dt.1[, tag2 := paste0(Chromosome, ':', Start_Position, ':', End_Position)]; 
		dt.1[, sample.class := factor(sample.class, levels=c('Primary', 'Tumor'))]
		# data table for number of rows in the big rect
		sq.dt = dcast(tag2 ~ sample.class, data=dt.1, value.var='code', drop=F, fun.aggregate = function(x){paste(unique(x), collapse=':')}); 
		sq.dt[is.na(Primary), Primary := adjustcolor('grey80')]
		sq.dt[is.na(Tumor),   Tumor   := adjustcolor('grey80')]; sq.dt
		nn 	= nrow(sq.dt)
		istep 	= 0.8 / nn; istep
		yy.seq	= seq(from = 0, to = .8, length.out = nn + 1);  yy.seq
		yy.seq	= yy.seq[1:(length(yy.seq)-1)]; yy.seq
		xleft1 	= rep(as.integer(dt.1$patient.name[1]), nn); xleft1
		ybottom1= unique(as.integer(dt.1$genes)) + yy.seq; ybottom1
		xright1 = xleft1 + .35; xright1
		ytop1 	= ybottom1 + istep; ytop1
		xleft2 	= rep(as.integer(dt.1$patient.name[1]), nn)+.45; xleft2
		ybottom2= ybottom1; ybottom2
		xright2 = xleft2 + .35; xright2
		ytop2 	= ytop1;
		tag2 	= sq.dt$tag2
		code1	= sq.dt[, Primary]; code1
		code2	= sq.dt[, Tumor]; code2
		gene 	= dt.1[1, genes]
		patient.name	= dt.1[1, patient.name]
		dt.xy.1  =data.table(
				     xleft1 = xleft1, ybottom1 = ybottom1, xright1 = xright1, ytop1 = ytop1, 
				     xleft2 = xleft2, ybottom2 = ybottom2, xright2 = xright2, ytop2 = ytop2,
				     code1 = code1, code2 = code2, tag2 = tag2, patient.name = patient.name, gene = gene) 
		dt.xy.1[, concordance := ifelse(code1 == code2, 'yes', 'no')]
		dt.xy.1
		dt.xy = rbind(dt.xy, dt.xy.1)
	}
	dt.xy[, color1 :=  color.table[dt.xy$code1, color.value]]
	dt.xy[, color2 :=  color.table[dt.xy$code2, color.value]]
	dt.xy[order(patient.name), ]

	n.genes = length(geneset); n.genes
	n.patients = length(levels(x.dt$patient.name)); 
	pdf(file=fname, width=12, height=8)
	plot('', type='n', xlim=c(0, n.patients + 4), ylim=c(0, n.genes + 2), xlab='', ylab='', xaxt = 'n', yaxt = 'n', bty = 'n', pch = '')
	# all the rect background color
	xleft = rep(1:max(as.integer(x.dt$patient.name)), times=n.genes);
	ybottom = rep(1:max(as.integer(x.dt$genes)), each=n.patients);
	xright = xleft +0.8;
	ytop = ybottom + .8;
	# gene names
	text(.5, 0.5+n.genes:1, label=gs, xpd=T, adj=c(1, 1))
	text(1:n.patients, n.genes + 1, label=levels(x.dt$patient.name), srt=90, adj=c(0, 1), xpd=T)
	rect(xleft = xleft, ybottom = ybottom, xright = xright, ytop = ytop, col = 'grey80', border=NA)
	rect(xleft = dt.xy$xleft1, ybottom = dt.xy$ybottom1, xright = dt.xy$xright1, ytop = dt.xy$ytop1, col = dt.xy$color1, border=NA)
	rect(xleft = dt.xy$xleft2, ybottom = dt.xy$ybottom2, xright = dt.xy$xright2, ytop = dt.xy$ytop2, col = dt.xy$color2, border=NA)
	dev.off()
	scp(fname)

}
