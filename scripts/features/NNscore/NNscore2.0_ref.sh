#!/bin/bash

## you should first install vina (http://vina.scripps.edu/download.html) and mgltools (http://mgltools.scripps.edu/downloads)
## you should assign the path of vina (vina_executable) in nnscore/NNScore2module.py
## then, add the path of nnscore in your PYTHONOATH. (~/data/scripts/features/NNscore/nnscore)

# it should be noted that only the python2 is supported because the scripts from mgltools are not compatibale with python3.


export PYTHONPATH=/home/shenchao/python_module2/lib/python2.7/site-packages
module purge
module load anaconda2
export PYTHONPATH=/home/shenchao/AI_based_SFs/software/NNscore/NNscore_2.0:${PYTHONPATH}
python NNscore2.0_ref.py

