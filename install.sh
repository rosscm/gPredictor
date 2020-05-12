#!/bin/bash

# Script to install dependencies/packages to run CRISPR guide analysis pipeline
# IMPORTANT: currently only supports Mac OS

# Get Homebrew if not already available on system
SCRIPT_DIR=$(pwd)
BREW_EXEC=$(which brew)
if [ -f "$BREW_EXEC" ]; then
  echo "Homebrew is already installed"
else
  /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
fi

brew update && brew upgrade

### ====================================
# Install brew packages
for pkg in python@2 git-lfs wget cmake; do
    if brew list -1 | grep -q "^${pkg}\$"; then
        echo "Package '$pkg' is installed"
    else
        echo "Package '$pkg' is not installed"
        brew install $pkg
    fi
done

# Install miniconda2 (can install anaconda instead if preferred)
## Both Anaconda and Miniconda use Conda as the package manager. The difference
## is that Miniconda only comes the package management system. So when you install
## it, there is just the management system and not coming with a bundle of
## pre-installed packages like Anaconda does.
echo '** INSTALLING MINICONDA2 / ANACONDA2 **'
CONDA_BASE=$(conda info --base)
if [ -d "$CONDA_BASE" ]; then
    echo "Miniconda2 / Anaconda2 exists"
else
  wget -c https://repo.anaconda.com/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
  chmod +x Miniconda2-latest-MacOSX-x86_64.sh
  ./Miniconda2-latest-MacOSX-x86_64.sh
  rm Miniconda2-latest-MacOSX-x86_64.sh
fi

# Install python modules within softwarae-specific virtual environments
if pip list | grep -F virtualenv; then
    echo "Package 'virtualenv' is installed"
else
    pip install virtualenv
fi

mkdir virtualenvs && cd virtualenvs

# Create Aziumth virtual environment
virtualenv azimuth
source azimuth/bin/activate
pip install pandas numpy scikit-learn==0.17.1
deactivate

# Create inDelphi virtual environment
virtualenv indelphi
source indelphi/bin/activate
pip install pandas scikit-learn==0.18.1 scipy numpy
deactivate

### ====================================
# Install CRISPR pipeline software
SOFTWARE_DIR=$SCRIPT_DIR/software

## CRISPRO
echo '** INSTALLING CRISPRO **'
### Create CRISPRO virtual environment
conda update -n base conda
env=$CONDA_BASE/envs/crispro
if [ -d "$env" ]; then
    echo "CRISPRO conda environment exists"
else
  conda env create -f $SCRIPT_DIR/virtualenvs/crispro_env.yml
fi

cd $SOFTWARE_DIR
CRISPRO=crispro-master
if [ -d "$CRISPRO" ]; then
    echo "$CRISPRO exists"
else
  conda activate crispro
  wget -c https://gitlab.com/bauerlab/crispro/-/archive/master/crispro-master.zip
  unzip crispro-master.zip
  rm crispro-master.zip
  cd crispro-master
  python setup.py install
  conda deactivate
fi

## ProTiler
echo '** INSTALLING ProTiler **'
cd $SOFTWARE_DIR
PROTILER=ProTiler-1.0.0
if [ -d "$PROTILER" ]; then
    echo "$PROTILER exists"
else
  git clone https://github.com/MDhewei/ProTiler-1.0.0.git
  cd ProTiler-1.0.0
  python setup.py install
fi

## inDelphi
echo '** INSTALLING inDelphi **'
cd $SOFTWARE_DIR
INDELPHI=inDelphi-model
if [ -d "$INDELPHI" ]; then
    echo "$INDELPHI exists"
else
  git clone https://github.com/maxwshen/inDelphi-model.git
fi

## Azimuth
echo '** INSTALLING Azimuth **'
## IMPORTANT: modifications needed to get this package running (Mac OSX issue)
### Solve 'RuntimeError: Python is not installed as a framework.' as follows
FILE=~/.matplotlib/matplotlibrc
if [ -f "$FILE" ]; then
    echo "Matplotlib backend exists"
else
  touch ${FILE}
  echo "backend: TkAgg" >> $FILE
fi

Azimuth=Azimuth
if pip list | grep -F "${Azimuth}"; then
    echo "Package '$Azimuth' is installed"
else
    echo "Package '$Azimuth' is not installed"
    pip install $Azimuth --user
fi

cd $SOFTWARE_DIR
if [ -d "$Azimuth" ]; then
    echo "$Azimuth exists"
else
  git clone https://github.com/MicrosoftResearch/Azimuth.git
  cd Azimuth
  python setup.py install
  ### Generating new model .pickle files - issue with incompatiblity with different
  ### versions of scikitlearn (https://github.com/MicrosoftResearch/Azimuth/issues/20)
  cd azimuth
  python model_comparison.py
  ### Run unit tests
  nosetests
fi

## FORECasT
echo '** INSTALLING FORECasT **'
cd $SOFTWARE_DIR
FORECAST=SelfTarget
if [ -d "$FORECAST" ]; then
    echo "$FORECAST exists"
else
  # Create a Python 3 virtual environment and activate it
  git clone https://github.com/felicityallen/SelfTarget.git
  cd SelfTarget/selftarget_pyutils
  pip install -e .
  cd ../indel_prediction
  pip install -e .
  cd ../
  ## Compile predictor
  cd indel_analysis/indelmap
  cmake . -DINDELMAP_OUTPUT_DIR=/usr/local/bin
  make && make install
  export INDELGENTARGET_EXE=/usr/local/bin/indelgentarget
fi

cd $SCRIPT_DIR
