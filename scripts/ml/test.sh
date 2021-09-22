#!/bin/bash


##random splitting
python xgb_training_and_predct.py -input_file ~/data/features/PDBbind-ReDocked/NNscore2.csv.bz2 \
-ref_file ~/data/features/PDBbind-ReDocked/rmsd_statistics.csv \
--rand 349 \
-f nnscore -vm random -max_evals 50 -o nnscore_random.csv -p -gpu 1

##3-fold CV
for i in $(seq 1 3) 
do
mkdir -p ccv${i}
cd ccv${i}
python ../xgb_training_and_predct.py -input_file ~/data/features/PDBbind-ReDocked/NNscore2.csv.bz2 \
-ref_file ~/data/features/PDBbind-ReDocked/rmsd_statistics.csv \
-pdb_prefix ~/data/scripts/utils/3-fold_CCV/subset -ccv_id ${i} \
-f nnscore -vm 3-CCV -max_evals 50 -o nnscore_ccv${i}.csv -p -gpu 1
cd ..
done

##refined-core splitting
python xgb_training_and_predct.py -input_file ~/data/features/PDBbind-ReDocked/NNscore2.csv.bz2 \
-ref_file ~/data/features/PDBbind-ReDocked/rmsd_statistics.csv \
--core_file ~/data/features/core_pdbid.txt \
-f nnscore -vm refined_core -max_evals 50 -o nnscore_refined-core.csv -p -gpu 1

##train on PDBbind-CrossDocked-Refined
python xgb_training_and_predct.py -input_file ~/data/features/PDBbind-CrossDocked-Refined/NNscore2_refinedcrossdock.csv.bz2 \
-ref_file ~/data/features/PDBbind-CrossDocked-Refined/rmsd_statistics.csv.bz2 \
-f nnscore -vm training_only -max_evals 50 -gpu 1

##  test on CASF-docking
python xgb_predict.py -input_file ~/data/features/CASF-Docking/NNscore2_casf.csv.bz2 \
-ref_file ~/data/features/CASF-Docking/rmsd_statistics_casf.csv \
-model ./final_best_model.pkl -f nnscore -t casf -gpu 1 -o xgb_out_casf.csv

python metrics2_resample.py -input_file ./xgb_out_casf.csv \
-ref_file ~/data/features/CASF-Docking/rmsd_statistics_casf.csv \
-t casf -d nnscore -resample -i 1000

##  test on PDBbind-CrossDocked-Core-g
python xgb_predict.py -input_file ~/data/features/PDBbind-CrossDocked-Core-g/NNscore2_casfcrossdockg.csv.bz2 \
-ref_file ~/data/features/PDBbind-CrossDocked-Core-g/glide_rmsd_statistics.csv \
-model ./final_best_model.pkl -f nnscore -t crosscore -gpu 1 -o xgb_out_crosscore.csv

##obtain the metrics focused on just the cross-docked poses
python metrics2_resample.py -input_file ./xgb_out_crosscore.csv \
-ref_file ~/data/features/PDBbind-CrossDocked-Core-g/glide_rmsd_statistics.csv \
-t crosscore-cross -d nnscore -resample -i 1000



