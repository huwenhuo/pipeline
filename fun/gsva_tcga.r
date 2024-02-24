
gsva.score
intersect(colnames(gsva.score), k3.cls$rn)

dsn <- model.matrix(~k3.cls$k3.cls.3); head(dsn)
colnames(design) <- c("ALL", "PS")
head(design)
fit <- lmFit(gsva.score.sel, dsn)
fit <- eBayes(fit)
gsva.res <- topTable(fit, coef="PS", number=Inf)
gsva.res

summary(decideTests(fit, p.value=0.05))
gsva.res = as.data.table(gsva.res, keep.rownames=T)
gsva.res[, Significance := 'Not']
gsva.res[P.Value < 0.05, Significance := 'Sig']
gsva.res

gg = ggplot(gsva.res, aes(x=reorder(rn, P.Value), y=logFC, color = Significance, size=-log(P.Value))) + geom_point() +
	theme(axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=0)) + ggtitle('Immune cell infilteration \nby ssGSEA (Primary vs Secondary) ') + 
	ylab("logFC") + xlab("Immune cells") + coord_flip() + 
	geom_hline(yintercept=0, color=adjustcolor(2, .3))
ggsave(gg, file='res/immune_cell_deconv_ssGSEA_primary_vs_secondary.pdf', width=5, height=5)
sync('res')

gsva.score.sel.melt = as.data.table(melt(gsva.score.sel))
head(gsva.anno.sel.dt)
gsva.score.sel.melt = merge(gsva.score.sel.melt, gsva.anno.sel.dt, by.x = 'Var2', by.y = 'rn')
gsva.score.sel.melt
gg = ggplot(gsva.score.sel.melt, aes(x=type, y=value)) + geom_jitter(width=.2) +
	ylab('Score in ssGSEA analyiss') + facet_wrap(~Var1, scales = 'free')  + 
	stat_summary(fun.data = "mean_cl_boot", colour = "red", size = .2) + 
	scale_y_continuous(trans='log2')
ggsave(gg, file='res/dotplots_gsva_score_primary_secondary.pdf', width=12, height=12)
sync('res')
