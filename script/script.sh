## generate bed file 
awk '{if ($3 != "gene") print $0;}' ~/program/Homo_sapiens.GRCh37.87.chr.gtf | grep -v "^#" | gtfToGenePred /dev/stdin /dev/stdout | genePredToBed stdin annotation.bed

rsync -pogavur --chmod=Du=rwx,Dg=rx,Do=rx,Fu=rw,Fg=r,Fo=r source/ dest/

