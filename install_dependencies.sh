# this file will install all dependencies needed for PyRocket
# run 'bash install_dependencies.sh' in your Ubuntu terminal

# ensure that apt is up to date
sudo apt-get update

# install C compiler
sudo apt-get install gcc

# install FORTRAN compiler
sudo apt-get install gfortran

# install gmsh
sudo apt-get install gmsh

# install pip
sudo apt-get install python3-pip

# install python libraries
pip3 install scipy
pip3 install thermo
pip3 install matplotlib
pip3 install fipy
pip3 install rocketcea