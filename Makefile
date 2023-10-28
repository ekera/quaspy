all:
	python3 -m pip install --upgrade build
	python3 -m build

install:
	pip3 install dist/quaspy-*.whl

uninstall:
	pip3 uninstall -y quaspy

reinstall:
	pip3 uninstall -y quaspy
	pip3 install dist/quaspy-*.whl

clean:
	rm -rf dist
	rm -rf src/quaspy.egg-info
