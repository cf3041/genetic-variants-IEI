import numpy as np
import os, glob
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import pandas as pd
import re
from matplotlib.colors import LinearSegmentedColormap

custom_cmap_colors = [(0.827, 0.184, 0.184), # LOF  '#D32F2F'
					  (1, 1, 1),			 # White at transition
					  (0.0, 0.474, 0.419)]   # GOF '#00796B'   

def main():

	sgRNA_IDs = np.loadtxt('./validation_sgRNA_IDs.txt', dtype=str)
	sgRNA_map = load_sgRNA_map('../sgRNA_ID_to_label_dict/sgRNA_ID_map.txt')
	# screen_data = pd.read_excel('./DP_vs_DN_Maxcyte_X2.sgrna_summary.xlsx')
	screen_data_D1 = pd.read_excel('./DP_vs_DN_Maxcyte_X2_by_donor.xlsx', sheet_name='D1')
	screen_data_D2 = pd.read_excel('./DP_vs_DN_Maxcyte_X2_by_donor.xlsx', sheet_name='D2')
	screen_data_D3 = pd.read_excel('./DP_vs_DN_Maxcyte_X2_by_donor.xlsx', sheet_name='D3')

	validation_data_D1 = screen_data_D1[np.isin(screen_data_D1['sgrna'], sgRNA_IDs)]
	validation_data_D2 = screen_data_D2[np.isin(screen_data_D2['sgrna'], sgRNA_IDs)]
	validation_data_D3 = screen_data_D3[np.isin(screen_data_D3['sgrna'], sgRNA_IDs)]

	# Ensure same order
	assert(np.all(np.equal(validation_data_D1['sgrna'],validation_data_D2['sgrna'])))
	assert(np.all(np.equal(validation_data_D1['sgrna'],validation_data_D3['sgrna'])))

	LFC_arr = np.vstack([np.vstack([validation_data_D1['LFC'],validation_data_D2['LFC']]), validation_data_D3['LFC']]).T

	# Calculate the row sums
	row_sums = np.sum(LFC_arr, axis=1)

	# Get the sorted row indices (highest to lowest row sum)
	sorted_indices = np.argsort(row_sums)[::-1]

	# Create figure and axis
	fig, ax = plt.subplots(figsize=(6, 10))  # Adjust figure size as needed

	# Display the heatmap using imshow
	# heatmap = ax.imshow(LFC_arr[sorted_indices, :], cmap='PiYG', aspect='auto', vmin=-4, vmax=4)

	# Custom cmap
	custom_cmap = LinearSegmentedColormap.from_list('Fig2_heatmap', list(zip([0, 0.5, 1], custom_cmap_colors)))
	heatmap = ax.imshow(LFC_arr[sorted_indices, :], cmap=custom_cmap, aspect='auto', vmin=-4, vmax=4)

	# Add colorbar
	cbar = fig.colorbar(heatmap, ax=ax, orientation='vertical', shrink=0.8)
	cbar.set_label('LFC', fontsize=12)

	# Set axis labels
	ax.set_xlabel('Donor', fontsize=12)
	# ax.set_ylabel('sgRNA', fontsize=12)

	# Customize tick positions and labels
	ax.set_xticks(np.arange(LFC_arr.shape[1]))
	ax.set_yticks(np.arange(LFC_arr.shape[0]))
	ax.set_xticklabels(['Donor 1', 'Donor 2', 'Donor 3'], fontsize=10)
	edit_labels = np.array([sgRNA_map[curr_ID] for curr_ID in validation_data_D1['sgrna'].to_numpy()[sorted_indices].astype(str)])
	ax.set_yticklabels(edit_labels, fontsize=6)

	# Rotate x-axis labels for better readability
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

	# Add grid lines for a cleaner look
	ax.set_xticks(np.arange(-0.5, LFC_arr.shape[1], 1), minor=True)
	ax.set_yticks(np.arange(-0.5, LFC_arr.shape[0], 1), minor=True)
	ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
	ax.tick_params(which="minor", bottom=False, left=False)

	pd.DataFrame(LFC_arr[sorted_indices, :], index=edit_labels, columns=['Donor 1', 'Donor 2', 'Donor 3']).to_excel('Fig2B_raw_data.xlsx')

	fig.savefig('./validation_sgRNAs_heatmap.png', bbox_inches='tight', dpi=800)
	fig.savefig('./validation_sgRNAs_heatmap.pdf', bbox_inches='tight')
	plt.close()

def load_sgRNA_map(in_file):

	out_dict = {}

	with open(in_file, 'r') as file:
		for line in file:
			key, value = line.strip().split("; ", 1)
			out_dict[key] = value

	return out_dict

if __name__=="__main__":
    main()
