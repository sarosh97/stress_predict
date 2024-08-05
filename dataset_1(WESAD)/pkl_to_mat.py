import numpy, scipy.io
import pickle, sys


subjects = ["S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","S13","S14","S15","S16","S17"]
base_file = "D:\\Research Papers and Data\\Wearable sensing\\WESAD\\"

for i in subjects:
    file_name = base_file + i + "\\" + i
    source_name = file_name + ".pkl"
    dest_name = file_name + ".mat"
	
    a=pickle.load( open( source_name, "rb" ), encoding='latin1')

    scipy.io.savemat(dest_name, mdict={'pickle_data': a})