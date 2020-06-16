# Makefile for AutoCV

cv: generate render

DOCKER_USERNAME = poldrack
current_dir = $(shell pwd)

clean:
	- rm autocv_template.aux
	- rm autocv_template.log
	- rm autocv_template.out
	- rm autocv_template.pdf
	- rm distinctions.tex
	- rm editorial.tex
	- rm education.tex
	- rm employment.tex
	- rm funding.tex
	- rm header.tex
	- rm memberships.tex
	- rm patents.tex
	- rm presentations.tex
	- rm pubs.tex
	- rm scholar.log
	- rm service.tex
	- rm talks.tex
	- rm teaching.tex
	- rm default.profraw

	
cv-local:
	python make_cv.py

# code to check environment variables
# from https://stackoverflow.com/questions/4728810/makefile-variable-as-prerequisite

guard-%:
	@ if [ "${${*}}" = "" ]; then \
    echo "Environment variable $* not set"; \
    exit 1; \
    fi

# use build arg to force installation of latest package from github
docker-build: guard-DOCKER_USERNAME
	docker build --build-arg DUMMY=`date +%s` -t $(DOCKER_USERNAME)/latex-python .

docker-deploy: docker-login docker-upload

docker-login: guard-DOCKER_USERNAME guard-DOCKER_PASSWORD
	docker login --username=$(DOCKER_USERNAME) --password=$(DOCKER_PASSWORD)

docker-upload: guard-DOCKER_USERNAME
	docker push $(DOCKER_USERNAME)/latex-python

shell: guard-DOCKER_USERNAME
	docker run --entrypoint /bin/bash -it -v $(current_dir):/data  $(DOCKER_USERNAME)/latex-python

render: guard-DOCKER_USERNAME
	docker run -it -v $(current_dir):/data  $(DOCKER_USERNAME)/latex-python

generate:
	docker run -it -v $(current_dir):/data  $(DOCKER_USERNAME)/latex-python python make_cv.py

pypi:
	python3 -m pip install --user --upgrade setuptools wheel
	python3 -m pip install --user --upgrade twine
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/*

test:
	python -m pytest --cov=autocv
