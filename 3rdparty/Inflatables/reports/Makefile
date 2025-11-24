all: 4_10_2019.html 6_8_2019.html

%.html: %.md
	pandoc --standalone --css=icg_report.css \
			--mathjax='https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_HTMLorMML' \
			--bibliography references.bib \
			$^ -o $@
