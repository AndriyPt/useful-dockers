# useful-dockers
Set of useful docker images which would be build on Docker Hub

Based on https://github.com/agolovynskyi1/latex-ukrainian-diser

## Install packages

texstudio
textlive all packages
scalable-cyrfonts-tex
pandoc
mupdf

## Software

https://pypi.org/project/scikit-fem/

ParaView

```bash
wget -O ParaView.tar.gz "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.13&type=binary&os=Linux&downloadFile=ParaView-5.13.0-MPI-Linux-Python3.10-x86_64.tar.gz"
sudo mkdir /opt/paraview
sudo tar -xvzf ./ParaView.tar.gz -C /opt/paraview
```

https://github.com/deepmodeling/jax-fem


## Installation of jax-fem dependencies

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


