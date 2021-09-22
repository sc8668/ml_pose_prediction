## Some representative scripts.

files  | descriptions 
---- | ----- 
./features |     
./features/E3FP/e3fp_features.py  | A demo to generate the E3FP fingerprint. 
./features/ECIF_and_ELEM/ecif_calculate.py  | The script to calculate the ECIF and ELEM features.
./features/ECIF_and_ELEM/PDB_Atom_Keys.csv  | The file recording the atom information of the protein for the calculation of ECIF.
./features/ECIF_and_ELEM/run_ecif.py  | A script to call ecif_calculate.py to calculate the ECIF features.
./features/ECIF_and_ELEM/run_elem.py  | A script to call ecif_calculate.py to calculate the ELEM features.
./features/NNscore/nnscore  | The module to calculate the NNscore features.
./features/NNscore/NNscore2.0_ref.py  | The script the calculate the NNscore features with the use of nnscore module.
./features/NNscore/vina  | A executable file of Vina.
./ml |  
./ml/test.sh   | A demo to reproduce the training and testing of the models in this study.
./ml/\*.py  | Some scripts utilized in the model training and evaluation.
./utils |  
./utils/3-fold_CCV/clustering_modified.py  | The script for 3-fold clustered cross validation.
./utils/3-fold_CCV/run_cluster.sh   | A demo to use the clustering_modified.py. 
./utils/3-fold_CCV/gnina_clust_info.txt  | A demo input file needed for the clustering.
./utils/3-fold_CCV/subset*   | A demo result of the clustering.
./utils/myalign2.py  | A script for the alignment of different protein with the structalign utility in Schrödinger.
./utils/run_prime-mmgbsa.py  | A script to perform Prime MM/GBSA in Schrödinger.
./utils/xscore_rescore.py    | A script to use X-Score for rescoring. 
