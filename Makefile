# fastperm --- common package tasks. Run `make` (or `make help`) for the list.
# Everything routes through devtools so the whole dependency graph is one tool.
# Slow simulations, if any, run only when NOT_CRAN=true, which the relevant
# targets set for you. fastperm is pure R --- there is no src/ to compile, so
# the clean target removes only R build artifacts and the local pkgdown site.
# The site builds into pkgdown-site/ (git-ignored), not docs/, because docs/
# holds the design memo; GitHub Actions publishes the site from main.

RSCRIPT := Rscript

.DEFAULT_GOAL := help
.PHONY: help deps document test test-fast check build install site coverage clean

help: ## List the available targets
	@grep -hE '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) \
	  | awk 'BEGIN{FS=":.*?## "}{printf "  make %-11s %s\n", $$1, $$2}'

deps: ## Install development and package dependencies
	$(RSCRIPT) -e 'if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools"); devtools::install_dev_deps(upgrade="never")'

document: ## Regenerate NAMESPACE and the man/*.Rd files from roxygen
	$(RSCRIPT) -e 'devtools::document()'

test: ## Run the full test suite, including any slow simulations
	NOT_CRAN=true $(RSCRIPT) -e 'devtools::test()'

test-fast: ## Run only the fast tests (skip the simulations)
	NOT_CRAN=false $(RSCRIPT) -e 'devtools::test()'

check: document ## R CMD check via devtools (the gate before a change is done)
	$(RSCRIPT) -e 'devtools::check()'

build: document ## Build the source tarball
	$(RSCRIPT) -e 'devtools::build()'

install: document ## Install fastperm into the local library
	$(RSCRIPT) -e 'devtools::install(upgrade=FALSE)'

site: ## Build the pkgdown site into pkgdown-site/ (needs pkgdown)
	$(RSCRIPT) -e 'pkgdown::build_site()'

coverage: ## Report test coverage (needs covr)
	NOT_CRAN=true $(RSCRIPT) -e 'print(covr::package_coverage())'

clean: ## Remove build artifacts and the local pkgdown site
	rm -rf ..Rcheck *.Rcheck *.tar.gz doc Meta pkgdown-site
