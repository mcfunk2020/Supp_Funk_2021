PROCESSED_DATA_DIR ?= processed_data

all: $(subst .mk,,$(wildcard *.mk))

install:
	Rscript -e 'remotes::install_local("sir")'

include *.mk

copy_processed_data: analysis_in_vivo
	- mkdir -p $(PROCESSED_DATA_DIR)
	find -L results_human/analysis_in_vivo -type f -name "*.tsv" -exec bash -c 'gzip < "$$0" > "'$(PROCESSED_DATA_DIR)'/$${0##*/}.gz"' {} \;
