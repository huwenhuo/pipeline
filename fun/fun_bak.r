## maftools/R/oncomatrix.R

createOncoMatrix = function(m, g = NULL, chatty = TRUE, add_missing = FALSE, cbio = FALSE, verbose = F){

  if(is.null(g)){
    stop("Please provde atleast two genes!")
  }

  subMaf = subsetMaf(maf = m, genes = g, includeSyn = FALSE, mafObj = FALSE)

  if(nrow(subMaf) == 0){
    if(add_missing){
      numericMatrix = matrix(data = 0, nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(numericMatrix) = g
      colnames(numericMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      oncoMatrix = matrix(data = "", nrow = length(g), ncol = length(levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])))
      rownames(oncoMatrix) = g
      colnames(oncoMatrix) = levels(getSampleSummary(x = m)[,Tumor_Sample_Barcode])

      vc = c("")
      names(vc) = 0

      return(list(oncoMatrix = oncoMatrix, numericMatrix = numericMatrix, vc = vc))
    }else{
      return(NULL)
    }
  }

  if(add_missing){
    subMaf[, Hugo_Symbol := factor(x = Hugo_Symbol, levels = g)]
  }


  cnv_events = c(c("Amp", "Del"), as.character(subMaf[Variant_Type == "CNV"][, .N, Variant_Classification][, Variant_Classification]))
  cnv_events = unique(cnv_events)

  if(cbio){
    vc = c("Nonstop_Mutation", "Frame_Shift_Del", "Missense_Mutation",
           "Nonsense_Mutation", "Splice_Site", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins")
    vc.cbio = c("Truncating", "Truncating", "Missense", "Truncating", "Truncating", "Truncating",
                "In-frame", "In-frame")
    names(vc.cbio) = vc
    subMaf[,Variant_Classification_temp := vc.cbio[as.character(subMaf$Variant_Classification)]]
    subMaf$Variant_Classification_temp = ifelse(test = is.na(subMaf$Variant_Classification_temp), yes = as.character(subMaf$Variant_Classification), no = subMaf$Variant_Classification_temp)
    subMaf[,Variant_Classification := as.factor(as.character(Variant_Classification_temp))]
    subMaf[,Variant_Classification_temp := NULL]
  }

  oncomat = data.table::dcast(data = subMaf[,.(Hugo_Symbol, Variant_Classification, Tumor_Sample_Barcode)], formula = Hugo_Symbol ~ Tumor_Sample_Barcode,
                              fun.aggregate = function(x, cnv = cnv_events){
                                #x = unique(as.character(x)) #>=2 distinct variant classification = Multi_Hit
                                x = as.character(x) # >= 2 same/distinct variant classification = Multi_Hit See #347
                                xad = x[x %in% cnv]
                                xvc = x[!x %in% cnv]

                                if(length(xvc)>0){
                                  xvc = ifelse(test = length(xvc) > 1, yes = 'Multi_Hit', no = xvc)
                                }

                                x = ifelse(test = length(xad) > 0, yes = paste(xad, xvc, sep = ';'), no = xvc)
                                x = gsub(pattern = ';$', replacement = '', x = x)
                                x = gsub(pattern = '^;', replacement = '', x = x)
                                return(x)
                              } , value.var = 'Variant_Classification', fill = '', drop = FALSE)

  #convert to matrix
  data.table::setDF(oncomat)
  rownames(oncomat) = oncomat$Hugo_Symbol
  oncomat = as.matrix(oncomat[,-1, drop = FALSE])

  variant.classes = as.character(unique(subMaf[,Variant_Classification]))
  variant.classes = c('',variant.classes, 'Multi_Hit')
  names(variant.classes) = 0:(length(variant.classes)-1)

  #Complex variant classes will be assigned a single integer.
  vc.onc = unique(unlist(apply(oncomat, 2, unique)))
  vc.onc = vc.onc[!vc.onc %in% names(variant.classes)]
  names(vc.onc) = rep(as.character(as.numeric(names(variant.classes)[length(variant.classes)])+1), length(vc.onc))
  variant.classes2 = c(variant.classes, vc.onc)

  oncomat.copy <- oncomat
  #Make a numeric coded matrix
  for(i in 1:length(variant.classes2)){
    oncomat[oncomat == variant.classes2[i]] = names(variant.classes2)[i]
  }

  #If maf has only one gene
  if(nrow(oncomat) == 1){
    mdf  = t(matrix(as.numeric(oncomat)))
    rownames(mdf) = rownames(oncomat)
    colnames(mdf) = colnames(oncomat)
    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  }

  #convert from character to numeric
  mdf = as.matrix(apply(oncomat, 2, function(x) as.numeric(as.character(x))))
  rownames(mdf) = rownames(oncomat.copy)


  #If MAF file contains a single sample, simple sorting is enuf.
  if(ncol(mdf) == 1){
    sampleId = colnames(mdf)
    mdf = as.matrix(mdf[order(mdf, decreasing = TRUE),])
    colnames(mdf) = sampleId

    oncomat.copy = as.matrix(oncomat.copy[rownames(mdf),])
    colnames(oncomat.copy) = sampleId

    return(list(oncoMatrix = oncomat.copy, numericMatrix = mdf, vc = variant.classes))
  } else{
    #Sort by rows as well columns if >1 samples present in MAF
    #Add total variants per gene
    mdf = cbind(mdf, variants = apply(mdf, 1, function(x) {
      length(x[x != "0"])
    }))
    #Sort by total variants
    mdf = mdf[order(mdf[, ncol(mdf)], decreasing = TRUE), ]
    #colnames(mdf) = gsub(pattern = "^X", replacement = "", colnames(mdf))
    nMut = mdf[, ncol(mdf)]

    mdf = mdf[, -ncol(mdf)]

    mdf.temp.copy = mdf #temp copy of original unsorted numeric coded matrix

    mdf[mdf != 0] = 1 #replacing all non-zero integers with 1 improves sorting (& grouping)
    tmdf = t(mdf) #transposematrix
    mdf = t(tmdf[do.call(order, c(as.list(as.data.frame(tmdf)), decreasing = TRUE)), ]) #sort

    mdf.temp.copy = mdf.temp.copy[rownames(mdf),] #organise original matrix into sorted matrix
    mdf.temp.copy = mdf.temp.copy[,colnames(mdf)]
    mdf = mdf.temp.copy

    #organise original character matrix into sorted matrix
    oncomat.copy <- oncomat.copy[,colnames(mdf)]
    oncomat.copy <- oncomat.copy[rownames(mdf),]

    ## modified
    mhit.name = as.numeric(names(variant.classes[variant.classes == 'Multi_Hit']))
    ir.names = rownames(mdf)
    ic.names = colnames(mdf)
    mdf.2 = copy(mdf)
    oncomat.copy.2 = copy(oncomat.copy)
    for(ir in 1:nrow(mdf)){
	    for(ic in 1:ncol(mdf)){
		    if(mdf[ir, ic] == mhit.name){
			    ir.name = ir.names[ir]; ir.name
			    ic.name = ic.names[ic]; ic.name
			    vc = unlist(m@data[Tumor_Sample_Barcode == ic.name & Hugo_Symbol == ir.name, simple.vc.onco])
			    sel.name = levels(vc)[min(as.numeric(vc))]; sel.name
			    sel.num = as.numeric(names(variant.classes[variant.classes == sel.name])); sel.num
			    if(verbose == T){
				    cat('ir:', ir, '\tic', ic, '\tic.name', ic.name, '\tir.name', ir.name, '\tsel.name:', sel.name, "\tsel.num:", sel.num, '\n')
			    }
			    mdf.2[ir, ic] = sel.num
			    oncomat.copy.2[ir, ic] = sel.name
		    }
	    }
    }

    variant.classes.2 = variant.classes[variant.classes != 'Multi_Hit']
    lst = list(oncoMatrix = oncomat.copy.2, numericMatrix = mdf.2, vc = variant.classes.2, cnvc = cnv_events)
    save(lst, file='lst.rdata')
    lst
  }
}
