import numpy as np
import pandas as pd
import os, glob
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import re

FDR_THRESH = 0.3
P_THRESH = 0.05
# ENRICH_HEX = '#00FF00'
# DEPLETE_HEX = '#FF0000'
ENRICH_HEX = '#00796B'
DEPLETE_HEX = '#D32F2F'

# Good annotation options:
# Black #000000
# Bright cyan #00FFFF
# Bright magenta #FF00FF
# Deep blue #0000FF

def main():

	plot_GOF = True
	plot_LOF = True

	in_txt = './DP_vs_DN_Maxcyte_X2.sgrna_summary.txt' 
	in_lib = 'PIK3CD_R1_ABE8e_tiling_library.xlsx'
	skeleton_in = './5DXU_chimX_skeleton.cxc'
	out_file = './5DXU_chimX_perspectives.cxc'

	with open(skeleton_in, 'r') as skel_fp:
		skel_lines = skel_fp.readlines()

	out_fp = open(out_file, 'w')
	out_fp.writelines(skel_lines)

	mageck_df = pd.read_csv(in_txt, delimiter='\t')
	lib_df = pd.read_excel(in_lib)

	LFC_arr = mageck_df['LFC'].to_numpy()

	max_LFC = LFC_arr.max()
	min_LFC = LFC_arr.min()
	assert(min_LFC < 0)
	
	colored_dict_3CD = {}
	colored_dict_3R1 = {}
	for curr_index, curr_row in mageck_df.iterrows():

		# Only color significant edits
		# if not curr_row['FDR'] < FDR_THRESH:
		if not curr_row['p.twosided'] < P_THRESH:
			continue

		curr_LFC = curr_row['LFC']

		# Get edit made by this sgRNA
		curr_edit = lib_df[lib_df['DICTIONARY'] == curr_row['sgrna']]['Amino acid edits'].to_numpy().astype(str)[0]

		if curr_edit == 'nan':
			continue

		if curr_row['Gene'] == 'PIK3CD':
			curr_chain = 'A'
		elif curr_row['Gene'] == 'PIK3R1':
			curr_chain = 'B'
		else:
			continue

		aas_to_colors = np.array([re.findall(r'\d+', item) for item in curr_edit.split(';') if item and not item.startswith('Exon')]).flatten().astype(int)		

		for aa_num in aas_to_colors:

			if curr_chain == 'A':
				if aa_num in colored_dict_3CD:
					if np.abs(colored_dict_3CD[aa_num]) > np.abs(curr_LFC):
						continue
					else:
						colored_dict_3CD[aa_num] = curr_LFC
			elif curr_chain == 'B':
				if aa_num in colored_dict_3R1:
					if np.abs(colored_dict_3R1[aa_num]) > np.abs(curr_LFC):
						continue
					else:
						colored_dict_3R1[aa_num] = curr_LFC

			if curr_LFC > 0 and plot_GOF:
				curr_max_hex = ENRICH_HEX
				curr_norm_LFC = curr_LFC / max_LFC

				curr_color = interpolate_color(curr_norm_LFC, curr_max_hex)
				out_fp.write('color /%s:%d %s\n' % (curr_chain, aa_num, curr_color))

			elif curr_LFC < 0 and plot_LOF:
				curr_max_hex = DEPLETE_HEX
				curr_norm_LFC = np.abs(curr_LFC) / np.abs(min_LFC)

				curr_color = interpolate_color(curr_norm_LFC, curr_max_hex)
				out_fp.write('color /%s:%d %s\n' % (curr_chain, aa_num, curr_color))


	# Write perspective coordinates
	out_fp.write('view matrix camera 0.83597,-0.4801,0.26579,88.694,0.22682,-0.13873,-0.964,-230.46,0.4997,0.86617,-0.0070714,6.5511\n')
	out_fp.write('turn x 20\n')
	out_fp.write('color /A:1021 #C4AC84\n') # PIK3CD: E1021K
	out_fp.write('color /B:564 #3B6EB3\n') # PIK3R1: N564K
	out_fp.write('save 5DXU_patient_samples.png width 3000 height 3000 transparentBackground true\n')

	out_fp.close()


# Convert hex color to RGB (normalized between 0 and 1)
def hex_to_rgb(hex_color):
	hex_color = hex_color.lstrip('#')
	return tuple(int(hex_color[i:i+2], 16) / 255.0 for i in (0, 2, 4))


# Convert RGB (0-1 range) to hex
def rgb_to_hex(rgb):
	return f"#{int(rgb[0] * 255):02X}{int(rgb[1] * 255):02X}{int(rgb[2] * 255):02X}"


def interpolate_color(value, upper_limit_hex):

	# Ensure the value is between 0 and 1
	value = max(0, min(1, value))

	# RGB for white (#FFFFFF) and the upper limit color (provided as hex)
	white_rgb = (1, 1, 1)
	upper_limit_rgb = hex_to_rgb(upper_limit_hex)

	# Interpolating between white and the upper limit color
	r = white_rgb[0] + value * (upper_limit_rgb[0] - white_rgb[0])
	g = white_rgb[1] + value * (upper_limit_rgb[1] - white_rgb[1])
	b = white_rgb[2] + value * (upper_limit_rgb[2] - white_rgb[2])

	# Return the interpolated color as a hex value
	return rgb_to_hex((r, g, b))


if __name__=="__main__":
    main()