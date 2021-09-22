# The impact of cross-docked poses on performance of machine learning classifier for protein-ligand binding pose prediction

<div align=center>
<img src="https://github.com/sc8668/ml_pose_prediction/blob/main/For_Table_of_Contents_Use_Only.jpg" width="500px" height="500px">
</div>

features.zip   				 --- The feature files utilized in this study.

models.zip					 --- Some models trained in this study.

scripts.zip    				 --- Some representative scripts.

example.zip					 --- A demo to use trained models to predict the binding poses of the in-house ligands

PDBbind-ReDocked.tar.bz2     --- PDBbind-ReDocked dataset.
			e.g.	--- 1a28/  --- the directory for "1a28".
						--- 1a28_prot/ 
							--- 1a28_p.pdb   the prepared protein file.
							--- 1a28_l.sdf   the prepared ligand file.
							--- 1a28_l.mol2   the prepared ligand file.
						--- 1a28_surflex/   --- the directory for Surflex docking results.
							--- each pose has been move the individual directory (i.e. 1_0000~1_0019)
							--- 1_000x is the rescoring results by Surflex-Dock (not directly scoring, the ligand will minimized).
					
PDBbind-CrossDocked-Core.tar.bz2  --- PDBbind-CrossDocked-Core dataset.
			--- dataset/ --- here the complexes with the same cluster are in the same directory.
				--- 1/ --- the cluster 1.
					--- 3ui7/ --- the directory for "3ui7".
						--- 3ui7_prot/ --- the directory contains the prepared protein and ligand for 3ui7.
			--- docking/ --- the directory for docking results.	
				--- 1_3ui7_3uuo  --- for cluster 1, the ligand of 3uuo is docked into the protein of 3ui7.
								 --- If the latter two IDs are same, re-docking is conducted; elsewise, cross-docking is conducted.
					--- 1_3ui7_3uuo_surflex  --- the directory for Surflex docking results.
									--- each pose has been move the individual directory (i.e. 1_0000~1_0019)
					--- 1_3ui7_3uuo_glideSP  --- the directory for Glide SP docking results. 
									--- each pose has been move the individual directory (i.e. 1_0000~1_0019)
					--- 1_3ui7_3uuo_vina     --- the directory for Autodock Vina docking results. 
									--- each pose has been move the individual directory (i.e. 1_0001~1_0020)		
									
PDBbind-CrossDocked-Refined.tar.bz2  --- PDBbind-CrossDocked-Refined dataset.
			--- dataset/ --- here the complexes with the same cluster are in the same directory.
				--- 1/ --- the cluster 1.
					--- 5d21/ --- the directory for "5d21".
						--- 5d21_prot/ --- the directory contains the prepared protein and ligand for 5d21.
			--- docking/ --- the directory for Surflex-Dock docking results.				
				--- here the docking poses for a single complex (e.g. 999_3ip5_3ip9) has been merged. (e.g. 999_3ip5_3ip9_surflex.mol2)

									
									
