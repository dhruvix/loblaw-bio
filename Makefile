PYTHON := .venv/bin/python
PIP := .venv/bin/pip
STREAMLIT := .venv/bin/streamlit

.PHONY: setup pipeline dashboard

setup:
	test -d .venv || python3 -m venv .venv
	$(PIP) install --upgrade pip
	$(PIP) install -r requirements.txt

pipeline: setup
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py

dashboard: setup
	$(STREAMLIT) run dashboard.py
