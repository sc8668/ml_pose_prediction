### The results are based on the poses provided by CASF2016 (http://www.pdbbind-cn.org/)

***lig_id = "xxx_native" represents the native poses.***          
     
     

						
files  | description   
---- | ----- 
rmsd_statistics_casf.csv | RMSD value of each pose. These RMSD values are directly obtained from the original CASF benchmark.
rmsd_statistics_casfmy.csv | The RMSD values here are calcalated using OpenBabel.
xscore_casf.csv.bz2 | the scores of X-Score ("average_score" is used as the final scores).
NNscore2_casf.csv.bz2 | the features of NNscore ("vina_affinity","vina_gauss_1","vina_gauss_2","vina_repulsion","vina_hydrophobic",and "vina_hydrogen"  are used as Vina features, and "vina_affinity" is used as the Vina score).
---- | -----
e3fp_casf.csv.bz2  | the features of E3FP.
ecif_casf.csv.bz2  | the features of ECIF.
elem_casf.csv.bz2 | the features of ELEM.
---- | ----- 
example/ | the example to assess the docking power in this study.
example/rmsd/ | the directory to store the RMSD values for each protein obtained from the original CASF benchmark.
example/rmsdmy/ | the directory to store the RMSD values for each protein calculated by OpenBabel.
example/scores/ | the directory to store the scores calculated by different SFs/models.	"xscore","vina","vina_xgb2","nnscore-vina_xgb2","nnscore_xgb2","elem_xgb2","ecif+vina_xgb2","ecif_xgb2","e3fp_xgb2" are calculated by ourselves, while the others are obtained from the original CASF benchmark.
example/CoreSet.dat | a file to store the PDB ID.
example/docking_power_withresample.py |a script the assess the docking power	

												
#### example
` ` `
python docking_power_withresample.py -c ./CoreSet.dat -s ./scores/X-Score 
-p positive -r ./rmsd -o X-Score.csv -l 2.0 -i 1000 -resample
` ` `
