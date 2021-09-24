###  The information of the datasets utilized in this study.
<br/>  
<br/>	    
<br/>	    
	    
	 
                                  Re-docked poses                                   Cross-docked poses

|Dataset| Complexes | Poses  | Positives  |  Negatives |  Complexes | Poses  | Positives | 	Negatives|
| ------ | :------: | :------: | :------: | :------: | :------: | :-----: | :------: | :------: |
| PDBbind-ReDocked	| 4,057	| 83,876	| 39,978 | 	43,898	| /	| /	| /	| /| 
| PDBbind-ReDocked-Refined	| 3,767	| 77,922	| 37,114	| 40,808	| /	| /	| /	| /| 
| PDBbind-ReDocked-Core| 	290	| 5,954 (5,664)<sup>a</sup>	| 2,864 (2,574)	| 3,090	| /	| /	| /	| /
| CASF-Docking	| 285<sup>b</sup>	| 22,777 (22,492)	| 5,494 (5,209)	| 17,283	| /	| /	| /	| /| 
| PDBbind-CrossDocked-Core-s	| 285	| 5,551	| 2,565	| 2,986	| 1,058	| 20,859	| 5,872	| 14,987| 
| PDBbind-CrossDocked-Core-g	| 282	| 4,795	| 1,596	| 3,199	| 1,030	| 17,814	| 3,768	| 14,046| 
| PDBbind-CrossDocked-Core-v	| 285	| 5,693	| 301	| 5,392	| 1,058	| 21,145	| 740	| 20,405| 
| PDBbind-CrossDocked-Refined	| 3,767	| 77,839 (74,072)	| 37,028 (33,261)	| 40,811	| 90,002	| 1,874,433	| 1,499,702	| 374,731| 
| PDBbind-CrossDocked-Refined<sup>c</sup>	| 3,767	| 77,839 (74,072)	| 37,028 (33,261)	| 40,811	| 90,002	| 1,731,351	| 1,428,161	| 303,190| 

<sup>**a**</sup>: The number in bracket refers to the number after removing the crystal poses.     
<sup>**b**</sup>: The core set of original PDBbind 2016 has 290 complexes belonging to 58 clusters, while only 285 are remained when constructing the CASF because there is a duplicated cluster.     
<sup>**c**</sup>: The set eliminates the cross-native poses.      
<br/> 
<br/>  
*	#### PDBbind-ReDocked dataset.
Each PDB has its own directory.    
e.g.       
1a28/  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory for "1a28".        
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/1a28_prot/     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/1a28_p.pdb   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the prepared protein file of "1a28".     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/1a28_l.sdf   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the prepared ligand file of "1a28".      
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/1a28_l.mol2  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the prepared ligand file of "1a28".      
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/1a28_surflex/   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory for Surflex docking results.     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; each pose has been moved into its individual directory (i.e. 1_0000~1_0019)     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; "1_000x" is the rescoring results by Surflex-Dock (the ligand will be minimized when scoring).     
<br/> 					
*	#### PDBbind-CrossDocked-Core dataset.
This set has two directories, one for protein and ligand (**dataset/**), and the other for docking results(**docking/**).                
e.g.       
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**dataset/** &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; here the complexes with the same cluster are in the same directory.             
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /1/ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the cluster 1.                 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /3ui7/ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the directory for "3ui7".               
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /3ui7_prot/ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory contains the prepared protein and ligand for "3ui7".                   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**docking/**  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory for docking results.	                 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /1_3ui7_3uuo/  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; for cluster 1, the ligand of "3uuo" is docked into the protein of "3ui7".       
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; If the latter two IDs are same, re-docking is conducted; elsewise, cross-docking is conducted.               
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/1_3ui7_3uuo_surflex/  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory for Surflex docking results.             
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; each pose has been moved to its individual directory (i.e. 1_0000~1_0019)              
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /1_3ui7_3uuo_glideSP/  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory for Glide SP docking results.           
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; each pose has been moved to its individual directory (i.e. 1_0000~1_0019)                 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /1_3ui7_3uuo_vina/     &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory for Autodock Vina docking results.             
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; each pose has been moved to its individual directory (i.e. 1_0001~1_0020)		             
<br/> 									
*	#### PDBbind-CrossDocked-Refined dataset.      
This set also has two directories, one for protein and ligand (**dataset/**), and the other for docking results(**docking/**).         
e.g.       
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**dataset/** &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; here the complexes with the same cluster are in the same directory.     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/1/ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the cluster 1.     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;/5d21/ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory for "5d21".     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; /5d21_prot/ &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory contains the prepared protein and ligand for "5d21".     
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;**docking/** &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; the directory for Surflex-Dock docking results.       
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;here the docking poses for a single complex (e.g. 999_3ip5_3ip9) has been merged. (e.g. 999_3ip5_3ip9_surflex.mol2)     




   
