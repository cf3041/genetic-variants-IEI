import numpy as np
import pandas as pd
import os, glob
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import re
from matplotlib.colors import LinearSegmentedColormap

LEN_HEX = '#991D1F'

def main():

	skeleton_in = './5DXU_chimX_skeleton.cxc'
	out_file_len = './5DXU_len_resistant.cxc'

	with open(skeleton_in, 'r') as skel_fp:
		skel_lines = skel_fp.readlines()

	out_file_len_fp = open(out_file_len, 'w')
	out_file_len_fp.writelines(skel_lines)
	
	# Annotate specific variants
	out_file_len_fp.write('color /B:567 %s\n' % LEN_HEX) # PIK3R1 K567G
	out_file_len_fp.write('color /B:573 %s\n' % LEN_HEX) # PIK3R1 L573P
	out_file_len_fp.write('color /B:570 %s\n' % LEN_HEX) # PIK3R1 L570P
	out_file_len_fp.write('color /B:583 %s\n' % LEN_HEX) # PIK3R1 Ex14-SA

	
	out_file_len_fp.write('delete /B:431-442\n')
	out_file_len_fp.write('view matrix camera 0.83597,-0.4801,0.26579,95.871,0.22682,-0.13873,-0.964,-256.49,0.4997,0.86617,-0.0070714,6.3601\n')
	out_file_len_fp.close()


if __name__=="__main__":
    main()
