BASENAME = manual

default:
	pdflatex  -shell-escape ${BASENAME}
###	bibtex    ${BASENAME}
### pdflatex  ${BASENAME}
	pdflatex  -shell-escape ${BASENAME}

.PHONY: clean

clean:
	rm -f *.aux *.log *.bbl *.blg *.brf *.cb *.ind *.idx *.ilg  \
	      *.inx *.dvi *.toc *lot *lof *.out *~ ~* spellTmp

