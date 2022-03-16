lint:
	black bibcodex --line-length 80
	flake8 bibcodex --ignore=E203
