NimbleUserManual.pdf :	NimbleUserManual.Rnw
	Rscript -e "library(knitr);  knit('MCMC.Rnw'); knit2pdf('NimbleUserManual.Rnw')"
# should we remove the compiled .o/.so in model1 too?

NimbleUserManualNoEval.pdf :	NimbleUserManual.Rnw
	Rscript makeNimbleUserManualNoEval.R 

clean :
	rm -f NimbleUserManual.pdf NimbleUserManual.tex NimbleUserManual.aux NimbleUserManual.log NimbleUserManual.toc NimbleUserManual.out MCMC.tex

