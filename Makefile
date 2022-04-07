check:
	R CMD build .
	R CMD check --as-cran nltm_`grep Version DESCRIPTION | cut -b 10-15`.tar.gz

dev-check:
	R-devel CMD build .
	R-devel CMD check --as-cran --use-valgrind --no-stop-on-test-error nltm_`grep Version DESCRIPTION | cut -b 10-15`.tar.gz

dev-build:
	R-devel CMD build .
