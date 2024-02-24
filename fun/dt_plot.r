dt.plot = function(idt){
	if(class(idt) != 'data.table'){
		idt = as.data.table(idt, keep.rownames=T)
	}
	idt = melt(idt)
	gg = ggplot(idt, aes(x=variable, y=value)) + 
		xlab('')
	gg
}
