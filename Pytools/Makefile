#PYT=python3

PYTOOLS_SRC = \
read_json.py \
setup.py \
west_clientserver.py \
west_fetch_input.py \
west_function3d.py \
west_read_bse.py \
west_utils.py

include ../west_make.inc

all: \
title \
install_timestamp

install_timestamp: $(PYTOOLS_SRC)
	$(PYT) -m pip install --user .
	@touch $@

clean:
	- rm -rf west.egg-info build dist
	- rm -f files.txt install_timestamp

title:
	@echo
	@echo "##############"
	@echo "### Pytools ##"
	@echo "##############"
	@echo
