To run the notebooks containing case studies:

1. Create and activate the case study conda environment:

```sh
conda env create -f environment.yml
source activate explosig-connect-case-study-env
```

2. Pull git submodules to download [ARDNMF](https://github.com/lrgr/ardnmf):

```sh
cd ..
git submodule update --init --recursive
git submodule foreach git pull origin master
cd -
```

3. Start jupyter:

```sh
jupyter notebook
```