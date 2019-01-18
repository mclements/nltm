check:
	R CMD build .
	R CMD check --as-cran nltm_1.4.2.tar.gz
