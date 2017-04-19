all: index.html navigating_nimble_objects.html
.PHONY: all clean

%.html: %.Rmd
	Rscript -e "library(rmarkdown); rmarkdown:::render('$<')"

clean:
	git clean -dfX
