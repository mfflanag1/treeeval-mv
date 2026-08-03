PYTHON ?= python3
export PYTHONPATH := $(CURDIR)/src

.PHONY: test calibrate search verify

test:
	$(PYTHON) -m unittest discover -s tests -v

calibrate:
	$(PYTHON) scripts/run_calibration.py

search:
	$(PYTHON) scripts/run_small_search.py

verify: test calibrate search
