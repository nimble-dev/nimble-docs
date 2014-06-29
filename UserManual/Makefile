NimbleUserManual.pdf :	NimbleUserManual.Rnw
	Rscript -e "library(knitr);  knit('MCMC.Rnw'); knit2pdf('NimbleUserManual.Rnw')"
# should we remove the compiled .o/.so in model1 too?

NimbleUserManualNoEval.pdf :	NimbleUserManual.Rnw
	Rscript makeNimbleUserManualNoEval.R 

clean :
	rm -f NimbleUserManual.pdf NimbleUserManual.tex NimbleUserManual.aux NimbleUserManual.log MCMC.tex
# for some reason rm -f NimbleUserManual.{pdf,tex,aux,log} not removing anything on Chris' linux box...
