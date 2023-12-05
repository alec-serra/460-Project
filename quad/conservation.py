import os
import numpy as np


FILES_DIR = '/media/alec/0cc34806-c13b-42cf-8529-9b13b4562d0f/Setup 1/2d_Tree_data_/'

def fuzzy_search(string):
	os.chdir(FILES_DIR)
	directory = os.listdir()
	res = []
	for fname in directory:
		if string in fname:
			res.append(fname)
	return res

def parse_file(handle):
	data = np.loadtxt(handle)
	return data

files = sorted(fuzzy_search('.dat'))

pc_to_m = 3.2407792896664e-17 
Solar_to_kg = 5.0289921396853e-31
G =  6.674e-11 * (pc_to_m)**3/(Solar_to_kg) 

res = []
for f in files:
	
	#states = parse_file(f)
	#print(states[:, 2:4])
	print(energy(states))
	exit(-1)

