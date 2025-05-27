import numpy as np
import os, glob
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import pandas as pd
import re
import seaborn as sns
import warnings
import scipy.stats as stats
import scipy.cluster.hierarchy as sch
from matplotlib.colors import LinearSegmentedColormap

# warnings.simplefilter("error", FutureWarning)

custom_cmap_colors = [(0, 0, 1),     # Blue at 0
         			  (1, 1, 1),     # White at transition
         			  (1, 1, 1),     # White at transition
          			  (1, 0, 0)]     # Red at 1

def main():

	screen_data = pd.read_excel('./DP_vs_DN_Maxcyte_X2.sgrna_summary.xlsx')
	FLOW_data = pd.read_excel('./BEv3_009_Flow_Analysis_withS6.xlsx')
	sgRNA_IDs = np.loadtxt('./validation_sgRNA_IDs.txt', dtype=str)
	sgRNA_map = load_sgRNA_map('../sgRNA_ID_to_label_dict/sgRNA_ID_map.txt')

	column_names = ['sgRNA', 'Substitution', 'LFC', 'Condition', 'pS6', 'pAKT 308', 'pAKT 473', 'p Value']
	plot_df = pd.DataFrame(columns=column_names)

	assert(np.all(FLOW_data['Condition_x']==FLOW_data['Condition_y']))
	assert(np.all(FLOW_data['Treatment_x']==FLOW_data['Treatment_y']))

	for curr_FLOW_ind, curr_FLOW_row in FLOW_data.iterrows():

		if 'FMO' in curr_FLOW_row['Condition_x'] or 'Unstained' in curr_FLOW_row['Condition_x']: # Skip unstained and FMO controls
			# print('Skipping FMO/unstained')
			continue

		if curr_FLOW_row['Condition_x'] == 'LacZ':
			
			curr_sgRNA_name = 'NT-sgRNA'
			curr_aa_substituion = 'NT-sgRNA'
			curr_LFC = 0
			curr_p_val = 0
		
		else:

			curr_gene = curr_FLOW_row['Condition_x'].split(':')[0]
			curr_aa_sub = curr_FLOW_row['Condition_x'].split(':')[1][1:]

			if curr_gene == 'PIK3CD':
				if curr_aa_sub == 'Glu1025Gly;Ser1026Gly;' or curr_aa_sub == 'Phe392Pro;Cys391Arg;':
					# print('Skipping %s' % curr_aa_sub)
					continue

			subset1 = screen_data[(screen_data['Gene']==curr_gene) & (screen_data['Amino acid edits (3-9 window)'].str.contains(curr_aa_sub, case=False))]
			correct_sgRNA_row = subset1.iloc[np.where(np.isin(subset1['sgrna'], sgRNA_IDs))[0]]
			
			if correct_sgRNA_row.shape[0] == 2:
				if curr_aa_sub == 'Tyr688Cys;':
					correct_sgRNA_row = correct_sgRNA_row[correct_sgRNA_row['sgrna'] == 'PIK3R1_1053']
				elif curr_aa_sub == 'Tyr688Cys;Ser689Gly;':
					correct_sgRNA_row = correct_sgRNA_row[correct_sgRNA_row['sgrna'] == 'PIK3R1_1055']
				else:
					bp()
			
			assert(correct_sgRNA_row.shape[0] == 1)

			curr_sgRNA_name = correct_sgRNA_row['sgrna'].to_numpy().astype(str)[0]
			curr_aa_substituion = correct_sgRNA_row['Amino acid edits (3-9 window)'].to_numpy().astype(str)[0]
			curr_LFC = correct_sgRNA_row['LFC'].to_numpy()[0]
			curr_p_val = -np.log10(correct_sgRNA_row['p.twosided'].to_numpy()[0])
		
		curr_condition = curr_FLOW_row['Treatment_x']
		curr_pS6 = curr_FLOW_row['Remove Doublets/Singlets/Lymphocytes/pS6 pos | Freq. of Parent']
		# curr_pS6 = curr_FLOW_row['Remove Doublets/Singlets/Lymphocytes | Median (pS6_235_236)']
		curr_pAKT_308 = curr_FLOW_row['Remove Doublets/Singlets/Lymphocytes | Median (pAKT_308)']
		curr_pAKT_472 = curr_FLOW_row['Remove Doublets/Singlets/Lymphocytes | Median (pAKT_473)']

		temp_df = pd.DataFrame([[curr_sgRNA_name,curr_aa_substituion,curr_LFC,curr_condition,curr_pS6,curr_pAKT_308,curr_pAKT_472,curr_p_val]], columns=column_names)
		if plot_df.empty:
			plot_df = temp_df
		else:
			plot_df = pd.concat([plot_df, temp_df], ignore_index=True)

	# Set a professional style
	sns.set_theme(style="whitegrid")  # Choose a grid-based professional style
	sns.set_context("paper", font_scale=1.5)  # Adjust context for publication-quality visuals

	# Ensure the 'Condition' column follows this order
	plot_df['Condition'] = pd.Categorical(plot_df['Condition'], categories=['Unstimulated', 'Stimulated', 'Leniolisib'], ordered=True)

	# Order by magnitude of stimulation
	stimulated_df = plot_df[plot_df['Condition'] == 'Stimulated']
	sorted_sgRNA = stimulated_df.sort_values(by='pAKT 473', ascending=False)['sgRNA']
	sgRNA_order = pd.unique(sorted_sgRNA.to_numpy()).astype(str)
	plot_df['sgRNA'] = pd.Categorical(plot_df['sgRNA'], categories=sgRNA_order, ordered=True)
	plot_df = plot_df.sort_values('sgRNA')

	plot_df_stim = plot_df[plot_df['Condition']=='Stimulated']

	plot_df = add_mapped_column(plot_df, sgRNA_map)
	plot_df_stim = add_mapped_column(plot_df_stim, sgRNA_map)

	plot_heatmap(plot_df_stim, './heatmap_validation_sgRNAs.pdf')
	plot_heatmap(plot_df_stim, './heatmap_validation_sgRNAs_no_map.pdf', with_cmap=False)

def plot_heatmap(df, save_pn, with_cmap=True):

	sgRNA_labels = df['sgRNA_map'].to_numpy().astype(str)
	# raw_mat = df[['LFC', 'pS6', 'pAKT 308', 'pAKT 473']].to_numpy()
	raw_mat = df[['LFC', 'pS6', 'pAKT 473']].to_numpy()

	# Normalize each column
	# display_mat = raw_mat / raw_mat.max(axis=0)
	sgRNA_inds = np.where(sgRNA_labels=='NT-sgRNA')[0]
	
	# Norm method 1
	# norm1 = 5*(raw_mat[:,0] / np.abs(raw_mat[:,0]).max())
	# norm2 = raw_mat[:, 1] / np.mean(raw_mat[sgRNA_inds,:][:,1])
	# norm3 = raw_mat[:, 2] / np.mean(raw_mat[sgRNA_inds,:][:,2])
	
	# Norm method 2
	# norm1 = raw_mat[:,0]  / np.abs(raw_mat[:,0]).max()
	# norm2 = raw_mat[:, 1] / np.abs(raw_mat[:,1]).max()
	# norm3 = raw_mat[:, 2] / np.abs(raw_mat[:,2]).max()

	# Norm method 3: project pFLOW data onto interval [0, 1]
	norm1 = raw_mat[:,0] / np.abs(raw_mat[:,0]).max()
	norm2 = (raw_mat[:,1]-raw_mat[:,1].min()) / (raw_mat[:,1].max()-raw_mat[:,1].min())
	norm3 = (raw_mat[:,2]-raw_mat[:,2].min()) / (raw_mat[:,2].max()-raw_mat[:,2].min())

	norm2_NT_mean = norm2[sgRNA_inds].mean()
	norm3_NT_mean = norm3[sgRNA_inds].mean()

	display_mat = np.vstack([np.vstack([norm1, norm2]), norm3]).T

	row_order = sch.leaves_list(sch.linkage(display_mat, method='ward'))

	# Create figure and axis
	fig, ax = plt.subplots(1, 3, figsize=(6, 10))  # Adjust figure size as needed

	# Load manually ordered heatmap
	if True:
		# pd.DataFrame(display_mat[row_order,:], index=sgRNA_labels[row_order]).to_excel('./ordered_display_mat.xlsx', header=False)
		temp_df = pd.read_excel('./ordered_display_mat.xlsx', header=None, index_col=0)

		del display_mat
		del row_order
		del sgRNA_labels

		display_mat = temp_df.to_numpy()
		row_order = np.arange(display_mat.shape[0])
		sgRNA_labels = temp_df.index.to_numpy().astype(str)

	# Display the heatmap using imshow
	NT_white_spread = 0.05
	heatmap1 = ax[0].imshow(display_mat[row_order, :][:, 0].reshape(row_order.size, 1), cmap='PRGn', aspect='auto', vmin=-0.5, vmax=0.5)
	cmap2 = LinearSegmentedColormap.from_list('heatmap2_cm', list(zip([0, norm2_NT_mean-NT_white_spread, norm2_NT_mean+NT_white_spread, 1], custom_cmap_colors)))
	heatmap2 = ax[1].imshow(display_mat[row_order, :][:, 1].reshape(row_order.size, 1), cmap=cmap2, aspect='auto', vmin=0, vmax=1)
	cmap3 = LinearSegmentedColormap.from_list('heatmap3_cm', list(zip([0, norm3_NT_mean-NT_white_spread, norm3_NT_mean+NT_white_spread, 1], custom_cmap_colors)))
	heatmap3 = ax[2].imshow(display_mat[row_order, :][:, 2].reshape(row_order.size, 1), cmap=cmap3, aspect='auto', vmin=0, vmax=1)

	# Add colorbar
	if with_cmap:
		cbar1 = fig.colorbar(heatmap1, ax=ax, shrink=0.4, location='bottom')
		cbar1.set_label('Screen LFC', fontsize=12)
		cbar2 = fig.colorbar(heatmap2, ax=ax, shrink=0.4, location='bottom')
		cbar2.set_label('pS6 Magnitude', fontsize=12)
		cbar3 = fig.colorbar(heatmap3, ax=ax, shrink=0.4, location='bottom')
		cbar3.set_label('pAKT Magnitude', fontsize=12)

	# Set axis labels
	# ax.set_xlabel('Donor', fontsize=12)
	# ax.set_ylabel('sgRNA', fontsize=12)

	# Customize tick positions and labels
	ax[0].set_xticks([0])
	ax[1].set_xticks([0])
	ax[2].set_xticks([0])

	ax[0].set_yticks(np.arange(display_mat.shape[0]))
	ax[1].set_yticks([])
	ax[2].set_yticks([])
	
	ax[0].set_xticklabels(['LFC'], fontsize=10)
	ax[1].set_xticklabels(['pS6'], fontsize=10)
	ax[2].set_xticklabels(['pAKT (pS473)'], fontsize=10)
	ax[0].set_yticklabels(sgRNA_labels[row_order], fontsize=4)

	# Rotate x-axis labels for better readability
	plt.setp(ax[0].get_xticklabels(), rotation=45, ha="right")
	plt.setp(ax[1].get_xticklabels(), rotation=45, ha="right")
	plt.setp(ax[2].get_xticklabels(), rotation=45, ha="right")

	# Add grid lines for a cleaner look
	ax[0].grid(False)
	ax[1].grid(False)
	ax[2].grid(False)
	# ax.set_xticks(np.arange(-0.5, display_mat.shape[1], 1), minor=True)
	# ax.set_yticks(np.arange(-0.5, display_mat.shape[0], 1), minor=True)
	# ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
	# ax.tick_params(which="minor", bottom=False, left=False)
	
	fig.savefig(save_pn, bbox_inches='tight', dpi=800)
	plt.close()

def add_mapped_column(in_df, sgRNA_map, col_name='sgRNA_map'):

	add_column = np.full(in_df.shape[0], '', dtype='<U100')
	count = 0
	for index, row in in_df.iterrows():

		if row['sgRNA'] == 'NT-sgRNA':
			add_val = 'NT-sgRNA'
		else:
			assert(row['sgRNA'] in sgRNA_map)

			add_val = sgRNA_map[row['sgRNA']]

		add_column[count] = add_val

		count += 1

	in_df[col_name] = add_column

	return in_df.copy()


def load_sgRNA_map(in_file):

	out_dict = {}

	with open(in_file, 'r') as file:
		for line in file:
			key, value = line.strip().split("; ", 1)
			out_dict[key] = value

	return out_dict

if __name__=="__main__":
    main()
