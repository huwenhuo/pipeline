
tcga_sc_id = c('TCGA-BT-A0YX-01A', 'TCGA-BT-A20U-01A', 'TCGA-BT-A2LD-01A', 'TCGA-C4-A0F1-01A', 'TCGA-C4-A0F7-01A', 'TCGA-CU-A0YN-01A', 
  'TCGA-DK-A2I2-01A', 'TCGA-FD-A3B5-01A', 'TCGA-FD-A3N5-01A', 'TCGA-G2-A2ES-01A', 'TCGA-G2-A3IB-01A', 'TCGA-GC-A3I6-01A', 
  'TCGA-GD-A3OS-01A', 'TCGA-BT-A42E', 'TCGA-GU-A42Q', 'TCGA-FD-A43Y', 'TCGA-FD-A5BU', 'TCGA-K4-A4AC', 'TCGA-FD-A5BY', 'TCGA-PQ-A6FI', 
  'TCGA-GU-A766', 'TCGA-CU-A72E', 'TCGA-E7-A7XN', 'TCGA-XF-A8HE', 'TCGA-YC-A89H', 'TCGA-E7-A97P', 'TCGA-XF-A8HH', 'TCGA-ZF-A9RE', 
  'TCGA-ZF-AA4W', 'TCGA-4Z-AA80', 'TCGA-4Z-AA82', 'TCGA-4Z-AA89', 'TCGA-XF-A9SJ', 'TCGA-XF-A9T4', 'TCGA-XF-A9T8', 'TCGA-ZF-AA53', 
  'TCGA-XF-AAME', 'TCGA-XF-AAMH', 'TCGA-XF-AAMT', 'TCGA-XF-AAN2', 'TCGA-XF-AAN5', 'TCGA-ZF-A9RD', 'TCGA-ZF-A9RG', 'TCGA-BT-A20X-01A', 
  'TCGA-FD-A3B4-01A', 'TCGA-FD-A3N6-01A')

tcga_sc_id = sub("-01A", "", tcga_sc_id)
tcga_sc_id

tcga_basal_id = c('TCGA-DK-A3WW', 'TCGA-SY-A9G5', 'TCGA-DK-A2I4', 'TCGA-E7-A7XN', 'TCGA-XF-A9T8', 'TCGA-XF-A9T5', 'TCGA-ZF-AA4V', 'TCGA-BT-A0YX', 'TCGA-DK-AA6S', 
  'TCGA-UY-A9PB', 'TCGA-4Z-AA7W', 'TCGA-4Z-AA84', 'TCGA-BT-A3PJ', 'TCGA-XF-A8HD', 'TCGA-ZF-AA4W', 'TCGA-BT-A20J', 'TCGA-FD-A3SS', 'TCGA-XF-AAN2', 
  'TCGA-XF-A9SY', 'TCGA-DK-A3IU', 'TCGA-G2-A2ES', 'TCGA-XF-A8HE', 'TCGA-XF-A9T3', 'TCGA-DK-A3IN', 'TCGA-DK-AA6L', 'TCGA-FD-A3SN', 'TCGA-FD-A3SO', 
  'TCGA-FD-A5BX', 'TCGA-GC-A3RC', 'TCGA-XF-A9SI', 'TCGA-HQ-A5NE', 'TCGA-K4-A5RJ', 'TCGA-K4-A6FZ', 'TCGA-BT-A3PK', 'TCGA-DK-AA6R', 'TCGA-FD-A3N5', 
  'TCGA-E7-A7DV', 'TCGA-UY-A8OB', 'TCGA-XF-AAN5', 'TCGA-XF-A9T6', 'TCGA-XF-A9SJ', 'TCGA-UY-A78P', 'TCGA-G2-A2EJ', 'TCGA-DK-A2I6', 'TCGA-FD-A3B6', 
  'TCGA-G2-A2EF', 'TCGA-XF-A9SM', 'TCGA-GC-A3YS', 'TCGA-ZF-A9RF', 'TCGA-CF-A1HS', 'TCGA-G2-AA3C', 'TCGA-4Z-AA7Q', 'TCGA-CU-A3KJ', 'TCGA-GC-A6I1', 
  'TCGA-DK-AA74', 'TCGA-E5-A2PC', 'TCGA-4Z-AA86', 'TCGA-4Z-AA81', 'TCGA-ZF-AA58', 'TCGA-DK-A1AB', 'TCGA-UY-A78L', 'TCGA-E7-A97P', 'TCGA-GU-A766', 
  'TCGA-ZF-AA54', 'TCGA-K4-A5RH', 'TCGA-GU-A42Q', 'TCGA-ZF-A9RN', 'TCGA-BT-A42F', 'TCGA-GC-A3I6', 'TCGA-BT-A20O', 'TCGA-ZF-AA53', 'TCGA-C4-A0F1', 
  'TCGA-BT-A42E', 'TCGA-ZF-A9RD', 'TCGA-GU-AATQ', 'TCGA-FD-A3NA', 'TCGA-FT-A61P', 'TCGA-UY-A9PA', 'TCGA-FD-A3B3', 'TCGA-XF-AAN3', 'TCGA-ZF-AA4N', 
  'TCGA-CU-A0YN', 'TCGA-FD-A6TH', 'TCGA-BL-A5ZZ', 'TCGA-GV-A40E', 'TCGA-GC-A3WC', 'TCGA-K4-A4AC', 'TCGA-PQ-A6FI', 'TCGA-FD-A3B5', 'TCGA-FD-A5C1', 
  'TCGA-FD-A6TD', 'TCGA-LC-A66R', 'TCGA-K4-A5RI', 'TCGA-FD-A3SP', 'TCGA-GU-A764', 'TCGA-E7-A3X6', 'TCGA-FT-A3EE', 'TCGA-DK-A6B5', 'TCGA-DK-A2I2', 
  'TCGA-FD-A6TK', 'TCGA-DK-A1AE', 'TCGA-FD-A3N6', 'TCGA-FD-A5BT', 'TCGA-DK-AA6Q', 'TCGA-4Z-AA82', 'TCGA-BT-A2LD', 'TCGA-GD-A3OQ', 'TCGA-GU-A762', 
  'TCGA-CU-A72E', 'TCGA-ZF-AA56', 'TCGA-BL-A13I', 'TCGA-UY-A8OC', 'TCGA-XF-A9SX', 'TCGA-XF-AAMT', 'TCGA-BT-A20V', 'TCGA-G2-A3IB', 'TCGA-FD-A62N', 
  'TCGA-HQ-A5ND', 'TCGA-C4-A0F0', 'TCGA-GV-A3QG', 'TCGA-DK-A3IM', 'TCGA-BT-A20U', 'TCGA-BT-A20X', 'TCGA-PQ-A6FN', 'TCGA-FD-A3B7', 'TCGA-FD-A5BU', 
  'TCGA-FD-A62S', 'TCGA-DK-AA6M', 'TCGA-FD-A3B4', 'TCGA-ZF-AA5H', 'TCGA-FD-A43Y', 'TCGA-ZF-AA5N', 'TCGA-FD-A5BY', 'TCGA-DK-A3WX', 'TCGA-XF-A9T4', 
  'TCGA-XF-AAN4', 'TCGA-XF-AAMW', 'TCGA-FD-A3B8', 'TCGA-FD-A5BS', 'TCGA-XF-AAME', 'TCGA-DK-A3WY', 'TCGA-XF-AAN8');


tcga_sc_basal_id = intersect(tcga_sc_id, tcga_basal_sc_id)
tcga_sc_basal_id

save(tcga_sc_basal_id, tcga_sc_id, tcga_basal_id, file="tcga_sc_basal_id.RData")

## base 47
## mcconkey

## cells paper classifiers
require(data.table)
tcga.blca.marker = data.table(
		     symbol =  c('KRT20', 'PPARG', 'FOXA1', 'GATA3', 'SNX31', 'UPK1A', 'UPK2', 'FGFR3', 'PGM5', 'DES', 'C7', 'SFRP4', 'COMP', 'SGCD', 'ZEB1', 'ZEB2', 'SNAI1', 'TWIST1', 'CDH2', 'CLDN3', 'CLDN4', 'CLDN7', 'CD44', 'KRT6A', 'KRT5', 'KRT14', 'COL17A1', 'DSC3', 'GSDMC', 'TGM1', 'PI3', 'TP63', 'CD274', 'PDCD1LG2', 'IDO1', 'CXCL11', 'L1CAM', 'SAA1', 'MSI1', 'PLEKHG4B', 'GNG4', 'PEG10', 'RND2', 'APLP1', 'SOX2', 'TUBB2B', 'CRTAC1', 'CTSE', 'PADI3', 'MSN', 'NR3C1', 'SHH', 'BMP5') ,
		     Class = c(rep('Luminal markers', 8), rep('ECM and smooth muscle', 6), rep('EMT and Claudin markers', 8), rep('Basal markers', 5), 
			       rep('Squamous markers', 5), rep('Immune markers', 6), rep('Neuronal-differentation', 8), rep('DEG in CIS', 5), rep('Sonic Hedgehog', 2)))

row.names(res) -> gene.id
gene.id = data.table(symbol = sub(".*_", "", gene.id), 
		     ensembl = sub("_.*", "", gene.id))

save(gene.id, file='~/program/fun/gene.id.RData')
load(file='~/program/fun/gene.id.RData')

tcga.blca.marker = merge(tcga.blca.marker, gene.id, by.x='symbol', by.y='symbol', all.x=T)
tcga.blca.marker

## cell 2018 TCGA bladder cancer paper
mut.2017 = c('TP53', 'RB1', 'RHOB', 'PIK3CA', 'KDM6A', 'TSC1', 'ELF3', 'KMT2D', 'CREBBP', 'CDKN1A', 'EP300', 'ZFP36L1', 'ARID1A', 'STAG2', 'CDKN2A', 'HRAS', 'KRAS', 'FBXW7', 'ERCC2', 'ASXL2', 'RHOA', 'KMT2A', 'FGFR3', 'NFE2L2', 'KMT2C', 'PSIP1', 'KANSL1', 'C3orf70', 'FAT1', 'SPTAN1', 'RXRA', 'ZBTB7B', 'PTEN', 'ATM', 'KLF5', 'PARD3', 'CUL1', 'NRAS', 'SF3B1', 'GNA13', 'RBM10', 'ACTB', 'MBD1', 'CASP8', 'HIST1H3B', 'TAF11', 'ERBB2', 'NUP93', 'SF1', 'ERBB3', 'METTL3', 'SPN', 'MB21D2', 'SSH3', 'USP28', 'ASXL1', 'TMCO4', 'HES1', 'ZNF773')
mut.2014 = c('TP53', 'MLL2', 'ARID1A', 'KDM6A', 'PIK3CA', 'EP300', 'CDKN1A', 'RB1', 'ERCC2', 'FGFR3', 'STAG2', 'ERBB3', 'FBXW7', 'RXRA', 'ELF3', 'NFE2L2', 'TSC1', 'KLF5', 'TXNIP', 'FOXQ1', 'CDKN2A', 'RHOB', 'FOXA1', 'PAIP1', 'BTG2', 'HRAS', 'ZFP36L', 'RHOA', 'CCND3')
