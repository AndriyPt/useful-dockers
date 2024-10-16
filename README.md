# useful-dockers
Set of useful docker images which would be build on Docker Hub

Based on https://github.com/agolovynskyi1/latex-ukrainian-diser

## Install packages

texstudio
texlive all packages
scalable-cyrfonts-tex
pandoc
mupdf

## Software

https://pypi.org/project/scikit-fem/

### ParaView

```bash
wget -O ParaView.tar.gz "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.13&type=binary&os=Linux&downloadFile=ParaView-5.13.0-MPI-Linux-Python3.10-x86_64.tar.gz"
sudo mkdir /opt/paraview
sudo tar -xvzf ./ParaView.tar.gz -C /opt/paraview
```

### Jax-Fem

https://github.com/deepmodeling/jax-fem


#### Installation of jax-fem dependencies

Install conda

```bash
curl -O https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh
bash Anaconda3-2024.06-1-Linux-x86_64.sh -b -p /opt/anaconda
```

```bash
conda env create -f environment.yml
conda activate jax-fem-env
```

```bash
python3 -m pip install --upgrade pip
python3 -m pip install jax
```

### NGSolve

Install `add-apt-repository` by running:
```bash
sudo apt install software-properties-common
``` 

First make sure to activate the 'universe' repository:
```bash
sudo apt-add-repository universe
```

Now activate the PPA (use ppa:ngsolve/nightly to get nightly builds instead of monthly releases):
```bash
sudo add-apt-repository ppa:ngsolve/ngsolve
sudo apt-get update
```

Now you can install NGSolve like any other Ubuntu package with
```bash
sudo apt-get install ngsolve
```

## Useful links

FEM Solvers

https://www.firedrakeproject.org/demos/helmholtz.py.html

https://github.com/NGSolve/ngsolve

https://github.com/deepmodeling/jax-fem

NGSolve BEM extension https://weggler.github.io/ngbem/intro.html 

## Working articles

Thermophysical and mechanical properties of biological tissues as a function of temperature: a systematic literature review 

https://www.tandfonline.com/doi/full/10.1080/02656736.2022.2028908#d1e2189

q_m => Specific metabolic rates of major organs and tissues across adulthood: evaluation by mechanistic model of resting energy expenditure1

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2980962/

Tissue physical properties ETH https://itis.swiss/virtual-population/tissue-properties/database/heat-generation-rate/

