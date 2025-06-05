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

pre-publish-check:
	python3 -m twine check dist/*

publish:
	python3 -m twine upload dist/*

publish-test:
	python3 -m twine upload -r testpypi dist/*

clean:
	rm -rf dist
	rm -rf src/quaspy.egg-info
