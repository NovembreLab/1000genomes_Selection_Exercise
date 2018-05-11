https://mybinder.org/v2/gh/jhmarcus/1000genomes_Selection_Exercise/master

# Workshop: 1000 Genomes Selection Exercise

A workshop on detecting positive selection.  The workshop introduces PBS and iHS approaches.  

It usable and instructive as is, but it is still in development. 

To get started run:
```
git clone https://github.com/NovembreLab/1000genomes_Selection_Exercise
cd 1000genomes_Selection_Exercise
```
Then
```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
chmod 755 Miniconda2-latest-Linux-x86_64.sh
./Miniconda2-latest-Linux-x86_64.sh
```
Agree to everything that the installation of conda asks you to do then run:
```
source ~/.bashrc
conda env create -f selection_exercise_environment.yml
```
And then: 
```
source activate selex
````


