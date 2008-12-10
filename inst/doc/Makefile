LogConcDEAD: LogConcDEAD.pdf clean

LogConcDEAD.pdf: LogConcDEAD.tex 
	pdflatex LogConcDEAD; bibtex LogConcDEAD; pdflatex LogConcDEAD; pdflatex LogConcDEAD

clean: 
	rm -f *.log *.aux *-*.eps *.bak *.rel *.blg Rplots.pdf *-*.pdf qhull_out.txt *.out *.bbl *.tex

