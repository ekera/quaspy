VENV = venv
PYTHON = $(VENV)/bin/python3
PIP = $(VENV)/bin/pip3

all: $(VENV) build

$(VENV):
	python3 -m venv $(VENV)
	$(PIP) install --upgrade build twine

build: $(VENV)
	$(PYTHON) -m pip install --upgrade build
	$(PYTHON) -m build

install:
	$(PIP) install dist/quaspy-*.whl

uninstall:
	$(PIP) uninstall -y quaspy

reinstall:
	$(PIP) uninstall -y quaspy
	$(PIP) install dist/quaspy-*.whl

pre-publish-check:
	$(PYTHON) -m twine check dist/*

publish:
	$(PYTHON) -m twine upload dist/*

publish-test:
	$(PYTHON) -m twine upload -r testpypi dist/*

clean:
	rm -rf dist
	rm -rf src/quaspy.egg-info

clean-all: clean
	rm -rf $(VENV)
