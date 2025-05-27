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
AA_ABBREVS = {
	'Ala' : 'A',
	'Arg' : 'R',
	'Asn' : 'N',
	'Asp' : 'D',
	'Cys' : 'C',
	'Glu' : 'E',
	'Gln' : 'Q',
	'Gly' : 'G',
	'His' : 'H',
	'Ile' : 'I',
	'Leu' : 'L',
	'Lys' : 'K',
	'Met' : 'M',
	'Phe' : 'F',
	'Pro' : 'P',
	'Ser' : 'S',
	'Thr' : 'T',
	'Trp' : 'W',
	'Tyr' : 'Y',
	'Val' : 'V',
	'Ter' : '*'
}

custom_cmap_colors = [(0, 0, 1),     # Blue at 0
         			  (1, 1, 1),     # White at transition
         			  (1, 1, 1),     # White at transition
          			  (1, 0, 0)]     # Red at 1

def main():

	FLOW_data = pd.read_excel('./BEv3_009point2_analysis.xlsx', sheet_name='CD3_ALL')
	pFLOW_data = pd.read_excel('./BEv3_009_Flow_Analysis.xlsx')
	TCF_data = pd.read_excel('TCF_all 121924.xlsx')
	S6p_data = pd.read_excel('./pS6 pos new 121924.xlsx')
	
	TCFp_label = 'Lymphocytes/Single Cells/Live/TCFhi | Freq. of Parent'
	TOX_label = 'Lymphocytes/Single Cells/Live | Mean (Comp-APC-A :: null)'
	AKT_label = 'Remove Doublets/Singlets/Lymphocytes | Median (pAKT_473)'
	S6p_label = 'Remove Doublets/Singlets/Lymphocytes/pS6 pos | Freq. of Parent'
	PD1p_label = 'Lymphocytes/Single Cells/Live/PD1 POS | Freq. of Parent'
	CTLA4p_label = 'Lymphocytes/Single Cells/Live/CTLA4 POS | Freq. of Parent'

	column_names = ['Sample', 'Treatment', 'TCF p', 'TOX', 'AKT', 'S6 p', 'PD1 p', 'CTLA-4 p']
	plot_df = pd.DataFrame(columns=column_names)

	sample_IDs = FLOW_data['Sample ID'].unique().astype(str)

	for sample_ID in sample_IDs:

		if 'FMO' in sample_ID: # Skip unstained and FMO controls for now
			# print('Skipping FMO')
			continue

		if 'Glu1025Gly;Ser1026Gly' in sample_ID or 'Phe392Pro;Cys391Arg;' in sample_ID:
			# print('Skipping %s' % sample_ID)
			continue

		if sample_ID == 'LacZ':
			
			curr_sample = 'NT-sgRNA'

		else:

			if sample_ID[-1] == ';':
				temp_ID = sample_ID[:-1]
			else:
				temp_ID = sample_ID

			curr_sample = temp_ID.replace(';', '_')

			if 'Exon' in curr_sample:
				curr_sample = curr_sample.replace('Exon', 'Ex')

			aa_abbreviations = re.findall(r'[A-Z][a-z][a-z]', curr_sample)

			for curr_abbreviation in aa_abbreviations:

				assert(curr_abbreviation in AA_ABBREVS)
				curr_sample = curr_sample.replace(curr_abbreviation, AA_ABBREVS[curr_abbreviation])

		curr_treatment = 'Stimulated'
		
		FLOW_rows = FLOW_data[(FLOW_data['Sample ID']==sample_ID) & (FLOW_data['Treatment']=='Stimulated')]
		assert(FLOW_rows.shape[0] == 3)

		pFLOW_rows = pFLOW_data[(pFLOW_data['Condition']==sample_ID) & (pFLOW_data['Treatment']=='Stimulated')]
		assert(pFLOW_rows.shape[0] == 3)

		TCF_rows = TCF_data[(TCF_data['Sample ID']==sample_ID) & (TCF_data['Treatment']=='Stimulated')]
		assert(TCF_rows.shape[0] == 3)

		S6p_rows = S6p_data[(S6p_data['Condition']==sample_ID) & (S6p_data['Treatment']=='Stimulated')]
		assert(S6p_rows.shape[0] == 3)
		

		temp_df = pd.DataFrame( {
				'Sample' : [curr_sample]*3,
				'Treatment' : [curr_treatment]*3,
				'TCF p' : TCF_rows[TCFp_label].to_numpy(), 
				'TOX' : FLOW_rows[TOX_label].to_numpy(), 
				'AKT' : pFLOW_rows[AKT_label].to_numpy(), 
				'S6 p' : S6p_rows[S6p_label].to_numpy(), 
				'PD1 p' : FLOW_rows[PD1p_label].to_numpy(), 
				'CTLA-4 p' : FLOW_rows[CTLA4p_label].to_numpy()
			})

		if plot_df.empty:
			plot_df = temp_df
		else:
			plot_df = pd.concat([plot_df, temp_df], ignore_index=True)

	# Set a professional style
	sns.set_theme(style="whitegrid")  # Choose a grid-based professional style
	sns.set_context("paper", font_scale=1.5)  # Adjust context for publication-quality visuals

	plot_heatmap(plot_df, './heatmap_phopheno_no_map.pdf', with_cmap=False)
	plot_heatmap(plot_df, './heatmap_phopheno.pdf')
	

def plot_heatmap(df, save_pn, with_cmap=True):

	sgRNA_labels = df['Sample'].to_numpy().astype(str)
	NT_sgRNA_inds = np.where(sgRNA_labels=='NT-sgRNA')[0]
	# column_names = ['TCF p', 'TOX', 'AKT', 'S6 p', 'PD1 p', 'CTLA-4 p']
	column_names = ['S6 p', 'AKT', 'PD1 p', 'CTLA-4 p', 'TOX', 'TCF p'] # Define column order for plotting
	raw_mat = df[column_names].to_numpy()

	# Normalize each column
	display_mat = np.full(raw_mat.shape, -1, dtype=raw_mat.dtype)
	NT_means = np.zeros(raw_mat.shape[1])
	for curr_col_ind in range(raw_mat.shape[1]):
		curr_col = raw_mat[:,curr_col_ind]
		display_mat[:,curr_col_ind] = (curr_col-curr_col.min()) / (curr_col.max()-curr_col.min())


		NT_means[curr_col_ind] = display_mat[:,curr_col_ind][NT_sgRNA_inds].mean()

		assert(display_mat[:,curr_col_ind].max()==1)
		assert(display_mat[:,curr_col_ind].min()==0)

	assert(np.where(display_mat==-1)[0].size==0)

	# Load manually ordered heatmap from Figure 3
	if True:

		temp_df = pd.read_excel('../Fig3_validation_FLOW_correlations/ordered_display_mat.xlsx', header=None, index_col=0)

		Fig3_labels = temp_df.index.to_numpy().astype(str)
		assert(np.all(np.isin(sgRNA_labels, Fig3_labels)))
		assert(sgRNA_labels.size==Fig3_labels.size)

		temp_indexes = np.unique(Fig3_labels, return_index=True)[1]
		unique_Fig3_labels = [Fig3_labels[index] for index in sorted(np.unique(Fig3_labels, return_index=True)[1])]

		row_order = []
		for unique_label in unique_Fig3_labels:

			curr_inds_to_add = np.where(sgRNA_labels==unique_label)[0]
			assert(curr_inds_to_add.size == 3)

			row_order.extend(curr_inds_to_add.tolist())			

		row_order = np.array(row_order, dtype=int)

		assert(row_order.size == sgRNA_labels.size)
		assert(np.all(np.char.equal(sgRNA_labels[row_order], Fig3_labels)))

	else:

		row_order = sch.leaves_list(sch.linkage(display_mat, method='ward'))


	# Create figure and axis
	fig, ax = plt.subplots(1, raw_mat.shape[1], figsize=(6, 10))  # Adjust figure size as needed

	# Display the heatmap using imshow
	NT_white_spread = 0.05
	heatmap_list = []
	for curr_col_ind in range(raw_mat.shape[1]):

		cmap_name = 'custom_cmap' + str(curr_col_ind)
		curr_cmap = LinearSegmentedColormap.from_list(cmap_name, list(zip([0, NT_means[curr_col_ind]-NT_white_spread, NT_means[curr_col_ind]+NT_white_spread, 1], custom_cmap_colors)))
		heatmap_list.append(ax[curr_col_ind].imshow(display_mat[row_order, :][:, curr_col_ind].reshape(row_order.size, 1), cmap=curr_cmap, aspect='auto', vmin=0, vmax=1))

	# Add colorbar
	if with_cmap:
		
		for curr_ind, curr_heatmap in enumerate(heatmap_list):

			curr_cbar = fig.colorbar(curr_heatmap, ax=ax, shrink=0.4, location='bottom')
			curr_cbar.set_label(column_names[curr_ind], fontsize=12)


	for curr_col_ind in range(raw_mat.shape[1]):

		ax[curr_col_ind].set_xticks([0])
		ax[curr_col_ind].set_xticklabels([column_names[curr_col_ind]], fontsize=8)
		ax[curr_col_ind].set_yticks([])
		ax[curr_col_ind].grid(False)

	ax[0].set_yticks(np.arange(display_mat.shape[0]))
	ax[0].set_yticklabels(sgRNA_labels[row_order], fontsize=4)

	# Rotate x-axis labels for better readability
	# plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

	df.to_excel('./Fig4A_raw_data.xlsx', index=False)
	pd.DataFrame(display_mat[row_order, :], index=sgRNA_labels[row_order], columns=column_names).to_excel('Fig4A_normalized_data.xlsx')

	fig.savefig(save_pn, bbox_inches='tight', dpi=800)
	plt.close()


if __name__=="__main__":
    main()
