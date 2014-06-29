NimbleDeveloperManual.pdf :	NimbleDeveloperManual.Rnw
	Rscript -e "library(knitr); knit2pdf('NimbleDeveloperManual.Rnw')"
# should we remove the compiled .o/.so in model1 too?

clean :
	rm -f NimbleDeveloperManual.pdf NimbleDeveloperManual.tex NimbleDeveloperManual.aux NimbleDeveloperManual.log
# for some reason rm -f NimbleDeveloperManual.{pdf,tex,aux,log} not removing anything on Chris' linux box...
