### The results are based on the re-docked poses generated by Surflex-Dock and the native poses.


***lig_id = "1_000x" represents the native poses***



files  | description   
---- | ----- 
rmsd_statistics.csv | RMSD value of each pose.
rmsd_statistics_surflex.csv | when using Surflex-Dock for rescoring, it will minimize the poses, so the RMSDs of those native poses not always equal to 0.
surflex_scores.csv  | the scores from Surflex-Dock ("Total_Score" is finally utilized).
prime-mmgbsa.csv.bz2 | the scores of Prime-MM/GBSA ("r_psp_MMGBSA_dG_Bind" is the final scores).
xscore.csv.bz2 | the scores of X-Score ("average_score" is used as the final scores).
NNscore2.csv.bz2 | the features of NNscore ("vina_affinity","vina_gauss_1","vina_gauss_2","vina_repulsion","vina_hydrophobic",and "vina_hydrogen" are used as Vina features, and "vina_affinity" is used as the Vina score).
---- | -----
e3fp.csv.bz2  | the features of E3FP.
ecif.csv.bz2 | the features of ECIF.
elem.csv.bz2 | the features of ELEM.





