JUPYTER_BOOK ?= $(shell [ -x ./jupyterbook/bin/jupyter-book ] && echo ./jupyterbook/bin/jupyter-book || echo jupyter-book)

build:
	$(JUPYTER_BOOK) build src/
publish:
	ghp-import -n -p -f src/_build/html
show: build
	(cd src/_build/html/ && open index.html &)
latex:
	$(JUPYTER_BOOK) build src/ --builder pdflatex
clean:
	$(JUPYTER_BOOK) clean src/
