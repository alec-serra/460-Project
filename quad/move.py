# import the PIL library for image manipulation
from PIL import Image
# import the imageio library for creating animations
import imageio
import numpy as np
import os 
# define the size of the image in pixels
width = 2500
height = 2500

# define the function that generates black and white pixels given an array of X, Y positions
def generate_pixels(positions):
	# create a new image object with mode '1' for black and white
	img = Image.new('1', (width, height))
	
	# get the pixel data of the image
	pixels = img.load()
	
	# loop over the positions
	for x, y in positions:
		int_x = int(x)
		int_y = int(y)
		if int_x > 0 and int_x < width and int_y > 0 and int_y < height:
			# set the pixel at (x, y) to white (value 1)
			pixels[int_x, int_y] = 1#(100//x, 200//y, 200/(x*y))
	# return the image object
	return img

# define the list of arrays of X, Y positions for animation
# each array is a list of tuples of (x, y) coordinates
# the coordinates are zero-based and start from the top-left corner
FILES_DIR = '/home/alec/Backup/2d_Tree_data_/'

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
	res = np.zeros((len(data), 2))
	for ndx, d in enumerate(data):
		res[ndx] = d[0:2]
	return res


def gen_images():
	file_list = sorted(fuzzy_search('.dat'))
	positions = [parse_file(p) for p in file_list]
	max_x = max([ max(p[:, 0]) for p in positions])
	min_x = min([ min(p[:, 0]) for p in positions])
	max_y = max([ max(p[:, 1]) for p in positions])
	min_y = min([ min(p[:, 1]) for p in positions])
	print(min_x, max_x, min_y, max_y)
	# create an empty list to store the images for animation
	images = []
	# loop over the list of positions
	for f, p in zip(file_list, positions):
		# generate an image with the given positions
		#Reshape data to be in picture
		p[:, 0] += (max_x - min_x)/2
		p[:, 0] *= width/(max_x - min_x)
		
		p[:, 1] += (max_y - min_y)/2
		p[:, 1] *= width/(max_y - min_y)
		#print(min(p[:, 0]), max(p[:, 0]),
		#		  min(p[:, 1]), max(p[:, 1]))
		ndx = f.split('_')[1]
		print(ndx)
		img = generate_pixels(p)
		# append the image to the list
		img.save('/home/alec/Backup/2d_Tree_Movie/{0}.png'.format(ndx))
def make_movie():
	# Define the folder path that contains the PNG images
	folder_path = "/home/alec/Backup/2d_Tree_Movie/"
	# Get the list of PNG files in the folder
	png_files = [f for f in os.listdir(folder_path) if f.endswith(".png")]
	# Sort the files by name
	png_files.sort()
	# Create a list of images
	images = []
	# Loop through the files and append them to the list
	for file in png_files:
		# Read the image file
		image = imageio.imread(os.path.join(folder_path, file))
		# Append the image to the list
		images.append(image)
	# Define the output file name
	output_file = "/home/alec/Backup/2d_Tree_Movie/output.gif"
	# Write the images to the output file
	imageio.mimsave(output_file, images, fps=20)
	for p in png_files:
		print(p)
		os.remove('/home/alec/Backup/2d_Tree_Movie/{0}'.format(p))
gen_images()
make_movie()
