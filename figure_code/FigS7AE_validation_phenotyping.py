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

def main():

	curr_sheet = 'CD3_ALL'
	# curr_sheet = 'CD8'
	# curr_sheet = 'CD4'

	FLOW_data = pd.read_excel('./BEv3_009point2_analysis.xlsx', sheet_name=curr_sheet)

	TCF_data = pd.read_excel('TCF_all 121924.xlsx')
	TCFp_label = 'Lymphocytes/Single Cells/Live/TCFhi | Freq. of Parent'

	Ki67_data = pd.read_excel('./BEv3_009point2_analysis.xlsx', sheet_name='Ki67_TCF')
	Ki67_label = 'Lymphocytes/Single Cells/Live/Ki67 hi | Freq. of Parent'

	title_text = 'CD3'
	TOX_label = 'Lymphocytes/Single Cells/Live | Mean (Comp-APC-A :: null)'
	LAG3_label = 'Lymphocytes/Single Cells/Live | Median (Comp-BV605-A :: null)'
	PD1_label = 'Lymphocytes/Single Cells/Live | Median (Comp-PE-A :: null)'
	CTLA4_label = 'Lymphocytes/Single Cells/Live | Median (Comp-PE-Cy7-A :: null)'
	CTLA4p_label = 'Lymphocytes/Single Cells/Live/CTLA4 POS | Freq. of Parent'
	PD1p_label = 'Lymphocytes/Single Cells/Live/PD1 POS | Freq. of Parent'

	# elif curr_sheet == 'CD8':
	# 	title_text = 'CD8'
	# 	TOX_label = 'Lymphocytes/Single Cells/Live/Q1: Comp-BV480-A- , Comp-BV421-A+ | Mean (Comp-APC-A :: null)'
	# 	LAG3_label = 'Lymphocytes/Single Cells/Live/Q1: Comp-BV480-A- , Comp-BV421-A+ | Median (Comp-BV605-A :: null)'
	# 	PD1_label = 'Lymphocytes/Single Cells/Live/Q1: Comp-BV480-A- , Comp-BV421-A+ | Median (Comp-PE-A :: null)'
	# 	CTLA4_label = 'Lymphocytes/Single Cells/Live/Q1: Comp-BV480-A- , Comp-BV421-A+ | Median (Comp-PE-Cy7-A :: null)'
	# 	CTLA4p_label = 'Lymphocytes/Single Cells/Live/Q1: Comp-BV480-A- , Comp-BV421-A+/CTLA4 POS | Freq. of Parent'
	# 	PD1p_label = 'Lymphocytes/Single Cells/Live/Q1: Comp-BV480-A- , Comp-BV421-A+/PD1 POS | Freq. of Parent'
	# elif curr_sheet == 'CD4':
	# 	title_text = 'CD4'
	# 	TOX_label = 'Lymphocytes/Single Cells/Live/Q3: Comp-BV480-A+ , Comp-BV421-A- | Mean (Comp-APC-A :: null)'
	# 	LAG3_label = 'Lymphocytes/Single Cells/Live/Q3: Comp-BV480-A+ , Comp-BV421-A- | Median (Comp-BV605-A :: null)'
	# 	PD1_label = 'Lymphocytes/Single Cells/Live/Q3: Comp-BV480-A+ , Comp-BV421-A- | Median (Comp-PE-A :: null)'
	# 	CTLA4_label = 'Lymphocytes/Single Cells/Live/Q3: Comp-BV480-A+ , Comp-BV421-A- | Median (Comp-PE-Cy7-A :: null)'
	# 	CTLA4p_label = 'Lymphocytes/Single Cells/Live/Q3: Comp-BV480-A+ , Comp-BV421-A-/CTLA4 POS | Freq. of Parent'
	# 	PD1p_label = 'Lymphocytes/Single Cells/Live/Q3: Comp-BV480-A+ , Comp-BV421-A-/PD1 POS | Freq. of Parent'
	# else:
	# 	bp()

	
	column_names = ['Sample', 'Treatment', 'TOX', 'LAG-3', 'PD1', 'CTLA-4', 'CTLA-4 p', 'PD1 p', 'TCF p', 'Ki67 p']
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

		# Add stimulated
		curr_treatment = 'Stimulated'
		
		FLOW_rows = FLOW_data[(FLOW_data['Sample ID']==sample_ID) & (FLOW_data['Treatment']==curr_treatment)]
		assert(FLOW_rows.shape[0] == 3)

		TCF_rows = TCF_data[(TCF_data['Sample ID']==sample_ID) & (TCF_data['Treatment']==curr_treatment)]
		assert(TCF_rows.shape[0] == 3)

		Ki67_rows = Ki67_data[(Ki67_data['Sample ID']==sample_ID) & (Ki67_data['Treatment']==curr_treatment)]
		assert(Ki67_rows.shape[0] == 3)

		temp_df = pd.DataFrame({
				'Sample' : [curr_sample]*3,
				'Treatment' : [curr_treatment]*3,
				'TOX' : FLOW_rows[TOX_label].to_numpy(),
				'LAG-3' : FLOW_rows[LAG3_label].to_numpy(),
				'PD1' : FLOW_rows[PD1_label].to_numpy(),
				'CTLA-4' : FLOW_rows[CTLA4_label].to_numpy(),
				'CTLA-4 p' : FLOW_rows[CTLA4p_label].to_numpy(),
				'PD1 p' : FLOW_rows[PD1p_label].to_numpy(),
				'TCF p' : TCF_rows[TCFp_label].to_numpy(),
				'Ki67 p' : Ki67_rows[Ki67_label].to_numpy()
			})

		if plot_df.empty:
			plot_df = temp_df
		else:
			plot_df = pd.concat([plot_df, temp_df], ignore_index=True)

		
		# Add leniolisib
		curr_treatment = 'Leniolisib'
		
		FLOW_rows = FLOW_data[(FLOW_data['Sample ID']==sample_ID) & (FLOW_data['Treatment']==curr_treatment)]
		assert(FLOW_rows.shape[0] == 3)

		TCF_rows = TCF_data[(TCF_data['Sample ID']==sample_ID) & (TCF_data['Treatment']==curr_treatment)]
		assert(TCF_rows.shape[0] == 3)

		Ki67_rows = Ki67_data[(Ki67_data['Sample ID']==sample_ID) & (Ki67_data['Treatment']==curr_treatment)]
		assert(Ki67_rows.shape[0] == 3)

		temp_df = pd.DataFrame({
				'Sample' : [curr_sample]*3,
				'Treatment' : [curr_treatment]*3,
				'TOX' : FLOW_rows[TOX_label].to_numpy(),
				'LAG-3' : FLOW_rows[LAG3_label].to_numpy(),
				'PD1' : FLOW_rows[PD1_label].to_numpy(),
				'CTLA-4' : FLOW_rows[CTLA4_label].to_numpy(),
				'CTLA-4 p' : FLOW_rows[CTLA4p_label].to_numpy(),
				'PD1 p' : FLOW_rows[PD1p_label].to_numpy(),
				'TCF p' : TCF_rows[TCFp_label].to_numpy(),
				'Ki67 p' : Ki67_rows[Ki67_label].to_numpy()
			})

		plot_df = pd.concat([plot_df, temp_df], ignore_index=True)

	# Set a professional style
	sns.set_theme(style="whitegrid")  # Choose a grid-based professional style
	sns.set_context("paper", font_scale=1.5)  # Adjust context for publication-quality visuals

	# Ensure the 'Condition' column follows this order
	plot_df['Treatment'] = pd.Categorical(plot_df['Treatment'], categories=['Stimulated', 'Leniolisib'], ordered=True)

	# Load order from figure 3
	Fig3_order = np.load('./Fig3_bar_order.npy')
	sgRNA_labels = plot_df['Sample'].to_numpy().astype(str)
	assert(np.all(np.isin(sgRNA_labels, Fig3_order)))
	assert(np.unique(sgRNA_labels).size==np.unique(Fig3_order).size)

	unique_Fig3_labels = [Fig3_order[index] for index in sorted(np.unique(Fig3_order, return_index=True)[1])]

	row_order = []
	for unique_label in unique_Fig3_labels:

		curr_inds_to_add = np.where(sgRNA_labels==unique_label)[0]
		assert(curr_inds_to_add.size == 6)

		row_order.extend(curr_inds_to_add.tolist())			

	row_order = np.array(row_order, dtype=int)

	assert(row_order.size == sgRNA_labels.size)

	plot_df = plot_df.iloc[row_order]

	

	plot_df_stim = plot_df[plot_df['Treatment']=='Stimulated']

	# Perform Dunnett test
	TOX_p_vals, TOX_conds = calculate_Dunnett_from_df(plot_df_stim, 'TOX')
	CTLA4_p_vals, CTLA4_conds = calculate_Dunnett_from_df(plot_df_stim, 'CTLA-4 p')
	PD1_p_vals, PD1_conds = calculate_Dunnett_from_df(plot_df_stim, 'PD1 p')
	TCF_p_vals, TCF_conds = calculate_Dunnett_from_df(plot_df_stim, 'TCF p')
	Ki67_p_vals, Ki67_conds = calculate_Dunnett_from_df(plot_df_stim, 'Ki67 p')
	
	assert(np.all(np.char.equal(TOX_conds, CTLA4_conds)))
	assert(np.all(np.char.equal(TOX_conds, PD1_conds)))
	assert(np.all(np.char.equal(TOX_conds, TCF_conds)))
	assert(np.all(np.char.equal(TOX_conds, Ki67_conds)))
	
	pd.DataFrame({'TOX p-value' : TOX_p_vals, 'CTLA-4 p-value' : CTLA4_p_vals,
				  'PD1 p-value' : PD1_p_vals, 'TCF p-value' : TCF_p_vals,
				  'Ki67 p-value' : Ki67_p_vals}, index=TOX_conds).to_excel('./Dunnett_results.xlsx')

	bp()

	plot_df.to_excel('./FigS7_raw_data.xlsx', index=False)

	plot_bar_graph(plot_df, 'Sample', 'TOX', 'Treatment', '', 'TOX MFI', title_text, './barplot_%s_TOX.pdf' % curr_sheet)
	# plot_bar_graph(plot_df, 'Sample', 'LAG-3', 'Treatment', '', 'LAG-3 MFI', title_text, './barplot_%s_LAG3.pdf' % curr_sheet)
	# plot_bar_graph(plot_df, 'Sample', 'PD1', 'Treatment', '', 'PD1 MFI', title_text, './barplot_%s_PD1.pdf' % curr_sheet)
	# plot_bar_graph(plot_df, 'Sample', 'CTLA-4', 'Treatment', '', 'CTLA-4 MFI', title_text, './barplot_%s_CTLA4.pdf' % curr_sheet)
	plot_bar_graph(plot_df, 'Sample', 'CTLA-4 p', 'Treatment', '', 'CTLA-4 Percent Positive', title_text, './barplot_%s_CTLA4p.pdf' % curr_sheet)
	plot_bar_graph(plot_df, 'Sample', 'PD1 p', 'Treatment', '', 'PD1 Percent Positive', title_text, './barplot_%s_PD1p.pdf' % curr_sheet)

	plot_bar_graph(plot_df, 'Sample', 'TCF p', 'Treatment', '', 'TCF Percent Positive', title_text, './barplot_%s_TCFp.pdf' % curr_sheet)
	plot_bar_graph(plot_df, 'Sample', 'Ki67 p', 'Treatment', '', 'Ki67 Percent Positive', title_text, './barplot_%s_Ki67p.pdf' % curr_sheet)

def calculate_Dunnett_from_df(in_df, col, cond_name = 'Sample', control_name = 'NT-sgRNA'):
	
	all_conds = in_df[cond_name].unique().astype(str)

	sample_list = []
	for curr_cond in all_conds:

		curr_vals = in_df[in_df[cond_name]==curr_cond][col].to_numpy()

		if curr_cond == control_name:
			control_list = curr_vals.tolist()
		else:
			sample_list.append(curr_vals.tolist())
	
	dunnett_results = stats.dunnett(*sample_list, control=control_list)

	return dunnett_results.pvalue, all_conds[all_conds != control_name]


def plot_bar_graph(df, x_val, y_val, hue_val, xlabel, ylabel, title, save_pn):

	fig, ax = plt.subplots(figsize=(18, 4))  # Adjust size for clarity
	sns.barplot(df, x=x_val, y=y_val, hue=hue_val, ax=ax, errorbar='sd', edgecolor='black', capsize=0.1)  # Add error bars and edges for better aesthetics
	sns.stripplot(data=df, x=x_val, y=y_val, hue=hue_val, dodge=True, ax=ax, marker="o", size=3, palette='dark:black', jitter=True, legend=False)
	ax.set_xlabel(xlabel, fontsize=14, labelpad=10)  # X-axis label
	ax.set_ylabel(ylabel, fontsize=14, labelpad=10)  # Y-axis label
	ax.set_title(title, fontsize=14)
	legend = ax.legend(title=hue_val, fontsize=12, title_fontsize=13, loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
	legend.get_frame().set_edgecolor('black')
	plt.xticks(rotation=45, horizontalalignment='right')
	fig.savefig(save_pn, bbox_inches='tight', dpi=800)
	plt.close()

def plot_scatterplot(df, x_val, y_val, hue_val, xlims, save_pn):

	df = df[df['LFC'] > xlims[0]]

	fig, ax = plt.subplots(figsize=(8, 6))
	# sns.scatterplot(df, x=x_val, y=y_val, hue=hue_val, size=hue_val, sizes=(20, 200), hue_norm=(0, 7), ax=ax)
	# sns.scatterplot(data=df, x=x_val, y=y_val, hue=hue_val, size=hue_val, sizes=(20, 200), hue_norm=(0, 7), ax=ax, palette='Blues', edgecolor="w", linewidth=0.5)
	sns.scatterplot(data=df, x=x_val, y=y_val, size=hue_val, sizes=(20, 200), hue_norm=(0, 7), ax=ax, color='blue', edgecolor="w", linewidth=0.5)
	ax.set_xlabel(x_val, fontsize=14, labelpad=10)  # X-axis label
	ax.set_ylabel(y_val, fontsize=14, labelpad=10)  # Y-axis label

	# Fit a linear regression line and overlay it
	sns.regplot(data=df, x=x_val, y=y_val, scatter=False, ax=ax, color='black', line_kws={'linewidth': 2}, truncate=True)

	# Calculate the regression statistics (slope, intercept, r-value, p-value, std err)
	slope, intercept, r_value, p_value, std_err = stats.linregress(df[x_val], df[y_val])

	# Annotate R² and p-value on the plot
	r_squared = r_value**2
	ax.annotate(f'R² = {r_squared:.2f}\np-value = {p_value:.2e}', 
	            xy=(0.05, 0.95), xycoords='axes fraction', ha='left', va='top', 
	            fontsize=12, color='black', backgroundcolor='white', weight='bold')

	legend = ax.legend(title=r'-log$_{10}$(p)', fontsize=12, title_fontsize=13, loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
	legend.get_frame().set_edgecolor('black')
	# ax.set_xlim(xlims)
	ax.grid(True, linestyle='--', alpha=0.7)
	ax.tick_params(axis='both', which='major', labelsize=12)
	fig.savefig(save_pn, bbox_inches='tight', dpi=800)
	plt.close()

def plot_heatmap(df, save_pn):

	sgRNA_labels = df['sgRNA'].to_numpy().astype(str)
	# raw_mat = df[['LFC', 'pS6', 'pAKT 308', 'pAKT 473']].to_numpy()
	raw_mat = df[['LFC', 'pS6', 'pAKT 473']].to_numpy()

	# Normalize each column
	# display_mat = raw_mat / raw_mat.max(axis=0)
	norm1 = 5*(raw_mat[:,0] / np.abs(raw_mat[:,0]).max())
	sgRNA_inds = np.where(sgRNA_labels=='NT-sgRNA')[0]
	norm2 = raw_mat[:, 1] / np.mean(raw_mat[sgRNA_inds,:][:,1])
	norm3 = raw_mat[:, 2] / np.mean(raw_mat[sgRNA_inds,:][:,2])

	display_mat = np.vstack([np.vstack([norm1, norm2]), norm3]).T

	row_order = sch.leaves_list(sch.linkage(display_mat, method='ward'))

	# Create figure and axis
	fig, ax = plt.subplots(figsize=(6, 10))  # Adjust figure size as needed

	# Display the heatmap using imshow
	heatmap = ax.imshow(display_mat[row_order, :], cmap='PRGn', aspect='auto', vmin=-3, vmax=3)

	# Add colorbar
	cbar = fig.colorbar(heatmap, ax=ax, orientation='vertical', shrink=0.8)
	cbar.set_label('Magnitude', fontsize=12)

	# Set axis labels
	# ax.set_xlabel('Donor', fontsize=12)
	ax.set_ylabel('sgRNA', fontsize=12)

	# Customize tick positions and labels
	ax.set_xticks(np.arange(display_mat.shape[1]))
	ax.set_yticks(np.arange(display_mat.shape[0]))
	# ax.set_yticks([])
	# ax.set_xticklabels(['LFC', 'pS6', 'pAKT 308', 'pAKT 473'], fontsize=10)
	ax.set_xticklabels(['Scaled LFC', 'Norm pS6', 'Norm pAKT 473'], fontsize=10)
	ax.set_yticklabels(sgRNA_labels[row_order], fontsize=4)

	# Rotate x-axis labels for better readability
	plt.setp(ax.get_xticklabels(), rotation=45, ha="right")

	# Add grid lines for a cleaner look
	ax.grid(False)
	# ax.set_xticks(np.arange(-0.5, display_mat.shape[1], 1), minor=True)
	# ax.set_yticks(np.arange(-0.5, display_mat.shape[0], 1), minor=True)
	# ax.grid(which="minor", color="gray", linestyle='-', linewidth=0.5)
	# ax.tick_params(which="minor", bottom=False, left=False)

	fig.savefig(save_pn, bbox_inches='tight', dpi=800)
	plt.close()

if __name__=="__main__":
    main()
