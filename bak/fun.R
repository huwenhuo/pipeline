filter_maf <- function(maf, vaf_cut, read_depth, exac=T, filters=T) {
  # filter on read depth and VAF
  maf.out <- maf[t_depth >= read_depth & ((t_var_freq > vaf_cut) | (t_var_freq > .02 & IMPACT_410 == TRUE))]
  # filter suspiciously long indel calls
  maf.out <- maf.out[nchar(Tumor_Seq_Allele1) < 50 & nchar(Tumor_Seq_Allele2) < 50]
  if (exac) {
    maf.out <- maf.out[(ExAC_AF < .0004 | is.na(ExAC_AF))]
  }
  if (filters) {
    if("FILTER" %in% names(maf.out)){
      maf.out <- maf.out[FILTER == "."]
    } else {
      cat("no FILTER column available in maf file\n")
    }
  }
  maf.out
}

#' Filter a maf file and produce a report of counts by Tumor_Sample_Barcode
#'
#' Each filter should correspond to a column in the report starting with "N_"
#' which counts the number of mutations passing this filter and all before it.
#' ie. the number should be no larger
#'
#' the final maf should correspond to the application of all filters
#'
#' @param maf
#' @param filters
#' @export filter_maf_report
filter_maf_report <- function(maf){

  ## code implementing each filter should be present for all subsequent filters
  ## in order to count the number of mutations pass each filter in turn
  ##
  ## use code from last combined filter to create filtered maf
  report <- maf[, list(
    N = .SD[,
            .N],
    N_minus_long_indels = .SD[nchar(Tumor_Seq_Allele1) < 50 & nchar(Tumor_Seq_Allele2) < 50,
                              .N],
    N_ExAC = .SD[nchar(Tumor_Seq_Allele1) < 50 & nchar(Tumor_Seq_Allele2) < 50 &
                   (ExAC_AF < .0004 | is.na(ExAC_AF)),
                 .N],
    N_coding = .SD[nchar(Tumor_Seq_Allele1) < 50 & nchar(Tumor_Seq_Allele2) < 50 &
                     (ExAC_AF < .0004 | is.na(ExAC_AF)) &
                     stringr::str_detect(Consequence,
                                         paste(c("synonymous_variant",
                                                 autospy::Nonsyn_Consequences),
                                               collapse = '|')),
                   .N]
    # ,N_nsyn = .SD[nchar(Tumor_Seq_Allele1) < 50 & nchar(Tumor_Seq_Allele2) < 50 &
    #                (ExAC_AF < .0004 | is.na(ExAC_AF)) &
    #                stringr::str_detect(Consequence,
    #                                    paste(
    #                                      autospy::Nonsyn_Consequences,
    #                                      collapse = '|')), .N]
  ),
  keyby = Tumor_Sample_Barcode]

  ############
  ## filter maf
  maf <- maf[nchar(Tumor_Seq_Allele1) < 50 & nchar(Tumor_Seq_Allele2) < 50 &
               (ExAC_AF < .0004 | is.na(ExAC_AF)) &
               stringr::str_detect(Consequence,
                                   paste(c("synonymous_variant",
                                           autospy::Nonsyn_Consequences),
                                         collapse = '|'))]

  ##
  column_definitions <- list("N" = "total number of mutations",
                            "N_minus_long_indels" = "...after indels length >50 removed",
                            "N_ExAC" = "...after removing common SNPs (ExAC_AF > 0.4%)",
                            "N_coding" = "...after requiring mutations are either synonymous_variant or in Nonsyn_Consequences"
                            #,"N_nsyn" = "...after requiring mutations have Nonsyn_Consequences"
  )

  ## check that all columns are defined
  if(!all(
    names(column_definitions) == setdiff(names(report), "Tumor_Sample_Barcode")
  )){
    cat("column descriptions:", paste(names(column_definitions)), "\n")
    cat("column names:", setdiff(names(report), "Tumor_Sample_Barcode"))
    stop("must define all columns")
  }


  ## check that there are the same number of mutations in the filtered maf as were reported passing the last filter
  if(!all(
    maf[, .N, keyby = Tumor_Sample_Barcode]$N == report[[ncol(report)]]
  )){
    stop("must have the same number of mutations in the filtered maf as were reported passing the last filter")
  }

  ## add diff columns to report
  report[, gsub("^N_", "fail_", names(report)[-1:-2]) := as.list(-diff(as.integer(.SD))),
         Tumor_Sample_Barcode]

  ## add total row to report
  report <- rbind(report,
                  c(Tumor_Sample_Barcode = "total",
                    as.list(colSums(report[, -1, with =F])))
  )

  list(maf = maf,
       report = copy(report),
       column_definitions = column_definitions)
}

plot_mutation_signatures <- function(sigs, pid, sid, primary = NA, fraction_threshold = .05) {

  dat <- as.data.table(sigs)
  #  setnames(dat, names(dat), c('sample_name','num_of_mutations', 1:30))
  setkeyv(dat, "Sample Name")
  dat.m <- melt(dat, id.vars= c('Sample Name', 'Number of Mutations'))
  setnames(dat.m, names(dat.m), c('Sample Name', 'Number of Mutations', 'signature', 'prob'))
  dat.m[, signature := as.integer(gsub("^Signature.", "", signature))]
  levels(dat.m$signature) <- c(levels(dat.m$signature), "31")
  dat.m[dat.m$prob <= fraction_threshold]$signature <- 31

  mmr_msi = colorRampPalette(c('#6BAED6', '#4292C6','#2171B5'))(5) # 6, 15, 20, 21, 26
  aid_apobec = rev(RColorBrewer::brewer.pal(8, 'Oranges'))[2:4] # 2, 9, 13
  t_c = rev(RColorBrewer::brewer.pal(8, 'Greens'))[2:4] # 5, 12, 16
  c_t = rev(RColorBrewer::brewer.pal(8, 'Purples'))[4:5] # 11, 23
  t_a = rev(RColorBrewer::brewer.pal(8, 'Reds'))[4:5] # 25, 27
  s = c(RColorBrewer::brewer.pal(12,'Set3')[c(1:3,7:12)], RColorBrewer::brewer.pal(8, 'Dark2')) # 1, 4, 7, 8, 14, 17, 18, 19, 22, 24, 28, 29, 30
  gg.cols = c(s[1],               #1
              aid_apobec[1],      #2
              s[2],               #3
              s[3],               #4
              t_c[1],             #5
              mmr_msi[1],         #6
              s[4],               #7
              s[5],               #8
              aid_apobec[2],      #9
              s[6],               #10
              c_t[1],             #11
              s[7],               #12
              aid_apobec[3],      #13
              s[8],               #14
              mmr_msi[2],         #15
              t_c[2],             #16
              s[9],               #17
              s[10],              #18
              s[11],              #19
              mmr_msi[3],         #20
              mmr_msi[4],         #21
              s[12],              #22
              c_t[2],             #23
              s[13],              #24
              t_a[1],             #25
              mmr_msi[5],         #26
              t_a[2],             #27
              s[14],              #28
              s[15],              #29
              s[16],              #30
              "darkgrey")         #other

  sig.order <- c(6,15,20,21,26,1,2,9,13,3,4,5,12,16,7,8,10,11,23,14,17,18,19,22,24,25,27,28,29,30,31)
  sig.names <- c("MMR/MSI (6)", "MMR/MSI (15)", "MMR/MSI (20)", "MMR/MSI (21)", "MMR/MSI (26)",
                 "Aging", "AID/APOBEC (2)", "AID/APOBEC (9)", "AID/APOBEC (13)", "BRCA1/2", "Smoking",
                 "Unknown T>C", "Unknown T>C", "Unknown T>C", "UV", "Unknown C>A", "POLE",
                 "Alkylating Agent", "Unknown C>T", "Unknown Hypermutated", "Unknown T>G", "Unknown C>A", "Unknown C>T",
                 "Aristolochic Acid", "Aflatoxin", "Unknown T>A", "Unknown T>A", "Unknown T>G", "Tobacco", "Unknown C>T",
                 "Other")
  sig.names <- sig.names[which(sig.order %in% dat.m$signature)]
  dat.m$`Sample Name` <- factor(dat.m$`Sample Name`, levels=unique(dat.m$`Sample Name`))
  dat.m[, sample_name := factor(stringr::str_extract(`Sample Name`, sid))]
  dat.m[, patient_name := factor(stringr::str_extract(`Sample Name`, pid))]

  if(!is.na(primary)){
    if(all(primary %in% unique(dat.m$`Sample Name`))){
      dat.m$`Sample Name` <- factor(dat.m$`Sample Name`, levels=c(primary, setdiff(unique(dat.m$`Sample Name`), primary)))
      dat.m$sample_name <- factor(dat.m$sample_name, levels = stringr::str_extract(levels(dat.m$`Sample Name`), sid))
    } else {
      warning(paste("Could not find primary file", primary))
    }
  }

  ## sum signatures labeled as 31 ("other") in to a single entry in the stacked bar
  dat.m <- dat.m[, list(prob = sum(prob)), list(sample_name, patient_name, `Number of Mutations`, signature)]

  sig.order <- sig.order[which(sig.order %in% dat.m$signature)]
  #  dat.m[, signature := factor(dat.m$signature, levels=sig.order, ordered = T)]

  dat.m[, signature := factor(dat.m$signature, levels=sig.order)]

  sig.plot <- ggplot(data=dat.m[order(signature)],
                     aes(x=sample_name,
                         y=prob,
                         fill=signature,
                         order=signature
                     )) +
    geom_bar(stat = 'identity', position='stack') +
    #    coord_flip() +
    scale_fill_manual(values=gg.cols[sig.order],
                      labels=sig.names,
                      name='Mutational Signature') +
    scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
    labs(x='Sample', y='Fraction of Mutations') +
    facet_grid(. ~ patient_name,
               space="free_x",
               scales="free_x",
               drop = TRUE) +
    theme_minimal(20) +
    theme(
      # plot.title=element_text(size=20, face = "bold"),
      axis.title=element_text(size=18, face = "bold"),
      # strip.text.y=element_text(size=15),
      axis.text.x=element_text(size=15, angle=45, face = "bold"),
      axis.text.y=element_text(size=15),
      legend.text=element_text(size=15),
      legend.title=element_text(size=15, face = "bold"),
      axis.ticks.x=element_blank(),
      axis.ticks.y=element_blank())

  if(uniqueN(dat.m$patient_name) == 1){
    sig.plot <- sig.plot + theme(
      strip.text.x=element_blank()
    )
  } else {
    sig.plot <- sig.plot + theme(
      strip.text.x = element_text(angle = 90)
    )
  }
  sig.plot <- sig.plot + geom_text(data = dat.m[signature == "31"],
                                   y = -0.03,
                                   aes(x=sample_name,
                                       label=`Number of Mutations`))
  sig.plot
}

plot_mutation_overlap <- function(maf, pid = NULL, log = F) {
	### create a binary matrix of mutations (columns: samples, rows: loci)
	dc <- dcast.data.table(maf,
			       TAG ~ Tumor_Sample_Barcode,
			       value.var = "Tumor_Sample_Barcode",
			       fun.aggregate = function(x) {
				       as.integer(length(x)>0)
			       },
			       fill = 0)
	mat <- as.matrix(dc[, -1, with = F])
	### take the pair-wise overlap of columns (Number of mutations in common)
	cm <- crossprod(mat)
	diag(cm) <- 0
	if(log == TRUE) {
		cm <- log10(cm+1)
		key.xlab = "log10( Nmutations )"
	} else {
		key.xlab = "Nmutations"
	}
	fill_palette <- colorRampPalette(c("white", "red"))(n = 299)
	if(!is.null(pid)){
		ipatient <- as.integer(factor(stringr::str_extract(rownames(cm), pid)))
		side_palette <- gg_color_hue(uniqueN(ipatient))[ipatient]
		gplots::heatmap.2(cm, trace = "none", col = fill_palette,
				  Rowv = NULL, Colv = NULL, lhei = c(2, 8), margins=c(10,10),
				  key.title = NA, key.ylab = NA, key.xlab = key.xlab,
				  RowSideColors = side_palette
				  )
	} else {
		gplots::heatmap.2(cm, trace = "none", col = fill_palette,
				  Rowv = NULL, Colv = NULL, lhei = c(2, 8), margins=c(10,10),
				  key.title = NA, key.ylab = NA, key.xlab = key.xlab
				  )
	}
}

make_mutation_overlap_plot <- function(maf, pid = NULL, log = F, out) {
  pdf(file=paste0(out, " <- mutation <- overlap", ".pdf"),
            width=16, height=16)
  plot_mutation_overlap(maf, pid, log)
    dev.off()
}

