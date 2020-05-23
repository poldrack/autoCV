# Makefile for AutoCV

cv: generate render

DOCKER_USERNAME = poldrack
current_dir = $(shell pwd)

cv-local:
	python make_cv.py

# code to check environment variables
# from https://stackoverflow.com/questions/4728810/makefile-variable-as-prerequisite

guard-%:
	@ if [ "${${*}}" = "" ]; then \
    echo "Environment variable $* not set"; \
    exit 1; \
    fi

docker-build: guard-DOCKER_USERNAME
	docker build -t $(DOCKER_USERNAME)/latex-python .

docker-deploy: docker-login docker-upload

docker-login: guard-DOCKER_USERNAME guard-DOCKER_PASSWORD
	docker login --username=$(DOCKER_USERNAME) --password=$(DOCKER_PASSWORD)

docker-upload: guard-DOCKER_USERNAME
	docker push $(DOCKER_USERNAME)/latex-python

shell: guard-DOCKER_USERNAME
	docker run -it -v $(current_dir):/data  $(DOCKER_USERNAME)/latex-python

shell: guard-DOCKER_USERNAME
	docker run -it -v $(current_dir):/data  $(DOCKER_USERNAME)/latex-python

render: guard-DOCKER_USERNAME
	docker run -it -v $(current_dir):/data  $(DOCKER_USERNAME)/latex-python

generate:
	docker run -it -v $(current_dir):/data  $(DOCKER_USERNAME)/latex-python python make_cv.py

