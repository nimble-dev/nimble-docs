NimbleUserManual.pdf :	NimbleUserManual.Rnw
	Rscript -e "library(knitr);  knit('MCMC.Rnw'); knit2pdf('NimbleUserManual.Rnw')"
# should we remove the compiled .o/.so in model1 too?

NimbleUserManualNoEval.pdf :	NimbleUserManual.Rnw
	Rscript makeNimbleUserManualNoEval.R 
# makeNimbleUserManualNoEval.R doesn't seem to exist
# one option here is to programmatically (via sed, e.g)
# change all chunks to have eval=FALSE and then run knit/knit2pdf

clean :
	rm -f NimbleUserManual.pdf NimbleUserManual.tex NimbleUserManual.aux NimbleUserManual.log NimbleUserManual.toc NimbleUserManual.out MCMC.tex

