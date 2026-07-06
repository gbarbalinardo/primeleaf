PY ?= python3

.PHONY: help install test verify constants census paper clean

help:
	@echo "install     editable install with dev dependencies"
	@echo "test        run the fast test suite"
	@echo "verify      reproduce Appendix A.1 and A.2 in full"
	@echo "constants   print the paper's critical constants"
	@echo "census      reproduce the tree-zero census (about 10 minutes)"
	@echo "paper       build paper/prime-leaf-tree-theory.pdf"
	@echo "clean       remove build artifacts and caches"

install:
	$(PY) -m pip install -e ".[dev]"

test:
	$(PY) -m pytest

verify:
	$(PY) -m primeleaf verify-combinatorics
	$(PY) -m primeleaf verify-zeta

constants:
	$(PY) -m primeleaf constants

census:
	$(PY) -m primeleaf census --t-max 500 --out data/zeros_t500.txt

paper:
	cd paper && latexmk -pdf -interaction=nonstopmode prime-leaf-tree-theory.tex

clean:
	cd paper && latexmk -C
	rm -rf .pytest_cache primeleaf/__pycache__ tests/__pycache__ *.egg-info
