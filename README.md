# Stress Predictiction using Machine Learning

There are four folders, each containing codes for respective datasets.

For each code you will need to download the cvxEDA algorithm (for feature extraction from EDA)
https://www.mathworks.com/matlabcentral/fileexchange/53326-cvxeda
and place the cvxEDA.m file in each folder

The detailed description of codes on how, or what, to run, is given below.


# Dataset-1 (WESAD)

Download the data at https://ubicomp.eti.uni-siegen.de/home/datasets/icmi18/
The data is in .pickle format so you need to run the python code in the file named "pkl_to_mat.py" to convert the data from "pickle" format to .mat format. Of course, this file needs to be executed in Python.
Then, in MATLAB run the file named "run_2_ext_feat" to extract features. You can adjust the window shift (located in the starting lines of the code). Make sure to change the filename accordingly.
After the features are extracted, run the file named "run_2_loso_classify" for LOSO classification using ML classifiers. Make sure to load the relevant feature file, for example if you have set 30 second shift then the file name would be "features_30.mat"

The files with names starting with "feat_ext_..." are functions for feature extraction of different modalities. 
The file named "load_wesad_data" is for loading the data. 


# Dataset-2 (UGI)

Download the data at https://github.com/italha-d/Stress-Predict-Dataset

You need to run three different files to extract features.
"run_2_ext_feat_interview"  "run_2_ext_feat_relax" and "run_2_ext_feat_stroop" are three files that must be executed to extract features from all three phases of this study. You can adjust the window shift (should be same in all files) and then run the file named "run_2_loso_classify" for LOSO classification. Make sure to load the relevant features, for example, in case you have set 5 second shift during feature extraction then the feature set is named "features_relax_5.mat". Similarly, for the Stroop and Relax cases.

In this dataset there are only 4 modalities, but for the sake of comparison with other datasets, only the EDA data is processed.

# Dataset-3 (UTD)

Dataset-3 can be downloaded at https://personal.utdallas.edu/~nourani/Bioinformatics/Biosensor_Data/

Run the file named "run_2_ext_feat" for feature extraction. You can adjust the window shift. Make sure to change the filename.
Run "run_2_loso_classify" for LOSO classification. Make sure to load the relevant feature file.

# Dataset-4 (CLAS)

Dataset-4 can be downloaded at https://www.wwwsensornetworkslab.com/clas (Needs an EULA form)

Run the file named "run_2_ext_feat" for feature extraction. You can adjust the window shift. Make sure to change the filename.
Run "run_2_loso_classify" for LOSO classification. Make sure to load the relevant feature file.
