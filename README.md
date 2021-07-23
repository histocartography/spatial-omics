# spatial-omics
[![Build Status](https://travis.ibm.com/art-zurich/spatial-omics.svg?token=bmUqdLriQp1g3yv7TJC6&branch=master)](https://travis.ibm.com/art-zurich/spatial-omics)
[![GitHub Pages](https://img.shields.io/badge/docs-sphinx-blue)](https://pages.github.ibm.com/art-zurich/spatial-omics)

Blueprint for a python package to create reproducible assets in a collaboration friendly manner.

## Usage

Click the "Use this template" button (on GitHub Enterprise) to create a new repo (if needed, dashes are recommended over underscores for this name).

Your package is installable for anyone with access to your repo.

### Adapting the spatialOmics for you: TODOs

Of course, clone your repo.

- [ ] [.travis.yml](.travis.yml) notifications: delete the section if undesired, else set up a slack channel and set the token.
- [ ] If you do not want to build documentation same as last step but also delete the `docs` folder. Else consider deleting the example markdown file.
- [ ] change the name "spatialOmics" to your desired package name. If really needed, underscores are valid. Note that the repo name is independent. Don't forget to change it in [docs/conf.py](docs/conf.py) and [docs/index.md](docs/index.md).
- [ ] Update author information in [setup.py](setup.py) and [docs/conf.py](docs/conf.py).
- [ ] If you decide against some checks, remove them from .travis.yml and the respective tools from the development requirements ([setup.py](setup.py) "dev" extras and [dev_requirements.txt](dev_requirements.txt)).
- [ ] Set up Travis, including a github token in case of pushing docs.
- [ ] If not building a docker image, remove the Dockerfile, the .travis directory and the "Docker" build stage in .travis.yml. Else you will have to implement .travis/deploy.sh at some point.
- [ ] Make this README file your own, and update the different banners.

Happy developing and documenting! 
### Install
```sh
# assuming you have a ssh key set up on GitHub
pip install "git+ssh://git@github.ibm.com/USERNAME_OR_ORGANIZATION/NEW_REPOSITORY.git@master"
```

see the [VCS](#vcs) paragraph how to handle requirements from GitHub or other version control system (VCS)

### Suggested setup for development

Create a `virtualenv`:

```sh
python -m venv venv
```

Activate it:

```sh
source venv/bin/activate
```

Install the package (including required packages) as "editable", meaning changes to your code do not require another installation to update.

```sh
pip install -e ".[dev]"
# or `pip install -e ".[vcs,dev]"`  # if you rely on other packages from github
```

To ensure exact versions:
```sh
pip install -r requirements.txt
pip install -r dev_requirements.txt
pip install -e .
```