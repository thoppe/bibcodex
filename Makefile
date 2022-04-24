package_name = bibcodex
coverage_omit="*/pubmed_parser*"

test:
	python -m pytest -s tests/

coverage:
	coverage run --omit=$(coverage_omit) --source=$(package_name) -m pytest tests/
	coverage report -m
	coverage html
	xdg-open htmlcov/index.html

lint:
	black bibcodex tests --line-length 80
	flake8 bibcodex --ignore=E203

clean:
	rm -rvf cover $(package_name).egg-info/ htmlcov dist *~

dist_test:
	rm -rvf dist
	python setup.py sdist
	twine upload -r test dist/*

dist_production:
	rm -rvf dist
	python setup.py sdist
	twine upload dist/*
