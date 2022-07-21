# Makefile for epydemic-signals
#
# Copyright (C) 2021 Simon Dobson
#
# This file is part of epydemic-signals, an experiment in epidemic processes.
#
# epydemic-signals is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# epydemic-signals is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with epydemic-signals. If not, see <http://www.gnu.org/licenses/gpl.html>.

# The name of our package on PyPi
PACKAGENAME = epydemic-signals

# The version we're building
VERSION = 0.1.1


# ----- Sources -----

# Source code
SOURCES_SETUP_IN = setup.py.in
SOURCES_CODE = \
	epydemic_signals/__init__.py \
	epydemic_signals/timeddict.py \
	epydemic_signals/sir_healing.py \
	epydemic_signals/sir_healing_one.py \
	epydemic_signals/hitting_healing.py \
	epydemic_signals/signal.py \
	epydemic_signals/signalgenerator.py \
	epydemic_signals/signaldynamics.py \
	epydemic_signals/stochasticsignaldynamics.py \
	epydemic_signals/synchronoussignaldynamics.py \
	epydemic_signals/compartment_signal.py \
	epydemic_signals/progress_signal.py \
	epydemic_signals/infection_boundary_signal.py
SOURCES_TESTS = \
	test/__init__.py \
	test/test_timeddict.py \
	test/test_signal.py \
	test/test_notebook.py \
	test/test_progress.py \
	test/test_boundary.py \
	test/progresssignalinvariants.py \
	test/test_stochasticsignalgenerator.py \
	test/test_synchronoussignalgenerator.py

SOURCES_DOC_CONF = doc/conf.py
SOURCES_DOC_BUILD_DIR = doc/_build
SOURCES_DOC_BUILD_HTML_DIR = $(SOURCES_DOC_BUILD_DIR)/html
SOURCES_DOC_ZIP = $(PACKAGENAME)-doc-$(VERSION).zip
SOURCES_DOCUMENTATION = \
	doc/index.rs
SOURCES_DIAGRAMS =
SOURCES_BIBLIOGRAPHY = doc/bibliography.bib

# Extras for building diagrams etc
SOURCES_UTILS =

# Extras for the build and packaging system
SOURCES_EXTRA = \
	README.rst \
	LICENSE \
	HISTORY \
	CONTRIBUTING.rst
SOURCES_GENERATED = \
	MANIFEST \
	setup.py

# Distribution files
DIST_SDIST = dist/$(PACKAGENAME)-$(VERSION).tar.gz
DIST_WHEEL = dist/$(PACKAGENAME)-$(VERSION)-py3-none-any.whl

# ----- Tools -----

# Base commands
PYTHON = python3
TOX = tox
COVERAGE = coverage
PIP = pip
TWINE = twine
GPG = gpg
GIT = git
ETAGS = etags
VIRTUALENV = $(PYTHON) -m venv
ACTIVATE = . $(VENV)/bin/activate
TR = tr
CAT = cat
SED = sed
RM = rm -fr
CP = cp
CHDIR = cd
ZIP = zip -r

# Files that are locally changed vs the remote repo
# (See https://unix.stackexchange.com/questions/155046/determine-if-git-working-directory-is-clean-from-a-script)
GIT_DIRTY = $(shell $(GIT) status --untracked-files=no --porcelain)

# Root directory
ROOT = $(shell pwd)

# Requirements for running the library and for the development venv needed to build it
VENV = venv3
REQUIREMENTS = requirements.txt
DEV_REQUIREMENTS = dev-requirements.txt

# Requirements for setup.py
# Note we elide dependencies to do with backporting the type-checking
PY_REQUIREMENTS = $(shell $(SED) -e '/^typing_extensions/d' -e 's/^\(.*\)/"\1",/g' $(REQUIREMENTS) | $(TR) '\n' ' ')

# Constructed commands
RUN_TESTS = $(TOX)
RUN_COVERAGE = $(COVERAGE) erase && $(COVERAGE) run -a setup.py test && $(COVERAGE) report -m --include '$(PACKAGENAME)*'
RUN_SETUP = $(PYTHON) setup.py
RUN_SPHINX_HTML = PYTHONPATH=$(ROOT) make html
RUN_TWINE = $(TWINE) upload dist/*


# ----- Top-level targets -----

# Default prints a help message
help:
	@make usage

# Build the tags file
tags:
	$(ETAGS) -o TAGS $(SOURCES_CODE) $(SOURCES_TESTS)

# Run tests for all versions of Python we're interested in
test: env Makefile setup.py
	$(ACTIVATE) && $(RUN_TESTS)

# Run coverage checks over the test suite
coverage: env
	$(ACTIVATE) && $(RUN_COVERAGE)

# Build the API documentation using Sphinx
.PHONY: doc
doc: env $(SOURCES_DOCUMENTATION) $(SOURCES_DOC_CONF)
	$(ACTIVATE) && $(CHDIR) doc && $(RUN_SPHINX_HTML)

# Build a development venv from the requirements in the repo
.PHONY: env
env: $(VENV)

$(VENV):
	$(VIRTUALENV) $(VENV)
	$(CAT) $(REQUIREMENTS) $(DEV_REQUIREMENTS) >$(VENV)/requirements.txt
	$(ACTIVATE) && $(PIP) install -U pip wheel && $(CHDIR) $(VENV) && $(PIP) install -r requirements.txt

# Build a source distribution
sdist: $(DIST_SDIST)

# Build a wheel distribution
wheel: $(DIST_WHEEL)

# Upload a source distribution to PyPi
upload: commit sdist wheel
	$(GPG) --detach-sign -a dist/$(PACKAGENAME)-$(VERSION).tar.gz
	$(ACTIVATE) && $(RUN_TWINE)

# Update the remote repos on release
commit: check-local-repo-clean
	$(GIT) push origin master
	$(GIT) tag -a v$(VERSION) -m "Version $(VERSION)"
	$(GIT) push origin v$(VERSION)

.SILENT: check-local-repo-clean
check-local-repo-clean:
	if [ "$(GIT_DIRTY)" ]; then echo "Uncommitted files: $(GIT_DIRTY)"; exit 1; fi

# Clean up the distribution build
clean:
	$(RM) $(SOURCES_GENERATED) $(SOURCES_DIST_DIR) $(PACKAGENAME).egg-info dist $(SOURCES_DOC_BUILD_DIR) $(SOURCES_DOC_ZIP) dist build

# Clean up everything, including the computational environment (which is expensive to rebuild)
reallyclean: clean
	$(RM) $(VENV)


# ----- Generated files -----

# Manifest for the package
MANIFEST: Makefile
	echo  $(SOURCES_EXTRA) $(SOURCES_GENERATED) $(SOURCES_CODE) | $(TR) ' ' '\n' >$@

# The setup.py script
setup.py: $(SOURCES_SETUP_IN) $(REQUIREMENTS) Makefile
	$(CAT) $(SOURCES_SETUP_IN) | $(SED) -e 's|VERSION|$(VERSION)|g' -e 's|REQUIREMENTS|$(PY_REQUIREMENTS)|g' >$@

# The source distribution tarball
$(DIST_SDIST): $(SOURCES_GENERATED) $(SOURCES_CODE) Makefile
	$(ACTIVATE) && $(RUN_SETUP) sdist

# The binary (wheel) distribution
$(DIST_WHEEL): $(SOURCES_GENERATED) $(SOURCES_CODE) Makefile
	$(ACTIVATE) && $(RUN_SETUP) bdist_wheel


# ----- Usage -----

define HELP_MESSAGE
Available targets:
   make test         run the test suite for all Python versions we support
   make coverage     run coverage checks of the test suite
   make tags         build the TAGS file
   make doc          build the API documentation using Sphinx
   make env          create a development virtual environment
   make sdist        create a source distribution
   make wheel        create binary (wheel) distribution
   make upload       upload distribution to PyPi
   make commit       tag current version and and push to master repo
   make clean        clean-up the build
   make reallyclean  clean up the virtualenv as well

endef
export HELP_MESSAGE

usage:
	@echo "$$HELP_MESSAGE"
