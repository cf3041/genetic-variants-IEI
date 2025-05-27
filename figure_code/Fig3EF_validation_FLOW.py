import numpy as np
import pandas as pd
import os, glob
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import re
from matplotlib.colors import LinearSegmentedColormap

AKT_HEX = '#D987A3'
S6_HEX = '#009E73'

def main():

	Fig3_df = pd.read_excel('./Fig3_plot_df.xlsx')

	skeleton_in = './5DXU_chimX_skeleton.cxc'
	out_file_AKT = './5DXU_AKT_validation.cxc'
	out_file_S6 = './5DXU_S6_validation.cxc'

	with open(skeleton_in, 'r') as skel_fp:
		skel_lines = skel_fp.readlines()

	out_file_AKT_fp = open(out_file_AKT, 'w')
	out_file_AKT_fp.writelines(skel_lines)

	out_file_S6_fp = open(out_file_S6, 'w')
	out_file_S6_fp.writelines(skel_lines)
	
	colored_dict_3CD = {}
	colored_dict_3R1 = {}
	stim_df = Fig3_df[Fig3_df['Condition']=='Stimulated']
	AKT_max = stim_df['pAKT 473'].max()
	len_df = Fig3_df[Fig3_df['Condition']=='Leniolisib']
	S6_max = len_df['pS6'].max()

	sgRNA_IDs = Fig3_df['sgRNA'].unique().astype(str)
	for sgRNA_ID in sgRNA_IDs:

		if sgRNA_ID == 'NT-sgRNA':
			continue

		curr_df = Fig3_df[Fig3_df['sgRNA']==sgRNA_ID]
		assert(np.where(curr_df['Condition']=='Stimulated')[0].size == 3)
		assert(np.where(curr_df['Condition']=='Leniolisib')[0].size == 3)

		curr_AKT_mag = curr_df[curr_df['Condition']=='Stimulated']['pAKT 473'].mean()
		curr_S6_mag = curr_df[curr_df['Condition']=='Leniolisib']['pS6'].mean()

		assert(curr_df['sgRNA_map'].unique().size==1)
		curr_sgRNA_map = curr_df['sgRNA_map'].to_numpy()[0]

		curr_gene = curr_sgRNA_map.split(':')[0]
		curr_edits = curr_sgRNA_map.split(' ')[1].split('_')

		if curr_gene == 'PIK3CD':
			curr_chain = 'A'
		elif curr_gene == 'PIK3R1':
			curr_chain = 'B'
		else:
			bp()

		aas_to_colors = [int(re.search(r'\d+', curr_edit).group()) for curr_edit in curr_edits]

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

			curr_AKT_norm_mag = curr_AKT_mag / AKT_max
			curr_AKT_color = interpolate_color(curr_AKT_norm_mag, AKT_HEX)
			out_file_AKT_fp.write('color /%s:%d %s\n' % (curr_chain, aa_num, curr_AKT_color))

			curr_S6_norm_mag = curr_S6_mag / S6_max
			curr_S6_color = interpolate_color(curr_S6_norm_mag, S6_HEX)
			out_file_S6_fp.write('color /%s:%d %s\n' % (curr_chain, aa_num, curr_S6_color))


	# Write perspective coordinates

	out_file_AKT_fp.write('view matrix camera 0.83597,-0.4801,0.26579,88.694,0.22682,-0.13873,-0.964,-230.46,0.4997,0.86617,-0.0070714,6.5511\n')
	out_file_AKT_fp.write('turn z 90\n')
	out_file_AKT_fp.write('turn x -10\n')
	out_file_AKT_fp.write('turn y 10\n')
	out_file_AKT_fp.write('save AKT_view1.png width 3000 height 3000 transparentBackground true\n')
	out_file_S6_fp.write('view matrix camera 0.83597,-0.4801,0.26579,88.694,0.22682,-0.13873,-0.964,-230.46,0.4997,0.86617,-0.0070714,6.5511\n')
	out_file_S6_fp.write('turn z 90\n')
	out_file_S6_fp.write('turn x -10\n')
	out_file_S6_fp.write('turn y 10\n')
	out_file_S6_fp.write('save S6_view1.png width 3000 height 3000 transparentBackground true\n')
	
	out_file_AKT_fp.close()
	out_file_S6_fp.close()

	print('Max AKT = %s' % str(AKT_max))
	print('Max S6 = %s' % str(S6_max))
	plot_custom_colorbar(AKT_HEX, [0, AKT_max], ['0', str(AKT_max)], './AKT_colorbar.pdf', colorbar_label='AKT Magnitude')
	plot_custom_colorbar(S6_HEX, [0, S6_max], ['0', str(S6_max)], './S6_colorbar.pdf', colorbar_label='S6 Magnitude')


def plot_custom_colorbar(color_hex, tick_positions, tick_labels, save_pn, colorbar_label="Colorbar"):
    """
    Plot a custom colorbar from 0 to a predefined color with custom tick marks.

    Parameters:
        color_hex (str): Hex color for the end of the colorbar (e.g., "#FF5733").
        tick_positions (list of float): Positions for tick marks (normalized between 0 and 1).
        tick_labels (list of str): Labels for the tick marks.
        colorbar_label (str): Label for the colorbar. Defaults to "Colorbar".
    """
    # Convert hex to RGB (normalized to 0-1 range)
    color_rgb = tuple(int(color_hex[i:i+2], 16) / 255.0 for i in (1, 3, 5))

    # Create a custom colormap from 0 (white) to the specified color
    cmap = LinearSegmentedColormap.from_list("custom_cmap", [(1, 1, 1), color_rgb])

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 1))

    # Create a gradient for the colorbar
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    ax.imshow(gradient, aspect="auto", cmap=cmap)

    # Remove axis ticks and labels for the image
    ax.set_axis_off()

    # Add a colorbar
    cbar = fig.colorbar(
        plt.cm.ScalarMappable(cmap=cmap),
        ax=ax,
        orientation="horizontal",
        fraction=0.046,
        pad=0.04,
    )

    # Set custom ticks and labels
    # cbar.set_ticks(tick_positions)
    # cbar.set_ticklabels(tick_labels)

    # Add label to the colorbar
    # cbar.set_label(colorbar_label)

    # Show the plot
    fig.savefig(save_pn, bbox_inches='tight')
    plt.close()


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
