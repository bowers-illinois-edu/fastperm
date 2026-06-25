# dev/ --- development notes and scratch (not shipped)

This directory holds development material that must not enter the package
tarball or the pkgdown site: design and theory notes, session handoffs, and
reproducible scratch/spike scripts that validate an idea but are not package
code. It is `.Rbuildignore`d, so `R CMD build` and `R CMD check` skip it.

The rule: everything here is either a non-shipping document or
throwaway-but-reproducible scratch. Once a script becomes real functionality
it moves to `R/` with tests in `tests/testthat/`; it does not stay here. The
package root is reserved for the standard files (DESCRIPTION, NAMESPACE,
NEWS.md) plus a few build-ignored infra files (Makefile, _pkgdown.yml) --- not
loose working files.

The Route B saddlepoint design and the riposte/fastperm cross-package handoff
currently live in the riposte repo as `FASTPERM_INTEGRATION.md`; its natural
long-term home is here.
