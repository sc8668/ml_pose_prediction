## An example to use trained models to predict the binding poses of the in-house ligands


Requirements：       
      python  3.6    
      openbabel 3.0.0    
     rdkit 2019.03.1    
      xgboost  0.90    
      sklearn  0.23.2    
      vina  1.1.2    
      mgltools  1.5.6   


#### 1. generate the features
**The MGLTOOLS_DIR, BABEL_DIR, VINA_EXEC, SCRIPT_DIR in generate_features.py should be first assigned.**

` ` `
python generate_features.py -p ~/input/3b27_p.pdb -l ~/input/3b27_docking.mol2 -f all -o out -bo
` ` `

#### 2. predict

` ` `
python xgb_predict2.py -input_file ./out_ecif.csv.bz2 
-input_file2 ./out_vina.csv.bz2 
-model ~/models/all-docked/ecif+vina_withrank/final_best_model.pkl 
-f ecif+vina+rank -gpu 0 -o xgb_out.csv -en 1
` ` `

#### 3. then, the users can select the best pose according to the predicted scores.





