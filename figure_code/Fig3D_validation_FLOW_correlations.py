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

	# plot_df.to_excel('./Fig3_plot_df.xlsx')

	# plot_heatmap(plot_df_stim, './heatmap_validation_sgRNAs.pdf')
	# plot_heatmap(plot_df_stim, './heatmap_validation_sgRNAs.png')

	# Perform Dunnett test
	AKT_p_vals, AKT_conds = calculate_Dunnett_from_df(plot_df_stim, 'pAKT 473')
	S6_p_vals, S6_conds = calculate_Dunnett_from_df(plot_df_stim, 'pS6')

	assert(np.all(np.char.equal(AKT_conds, S6_conds)))
	bp()
	pd.DataFrame({'pAKT p-value' : AKT_p_vals, 'S6 p-value' : S6_p_vals}, index=AKT_conds).to_excel('./Dunnett_results.xlsx')

	plot_bar_graph(plot_df, 'sgRNA_map', 'pS6', 'Condition', './barplot_sgRNA_vs_S6.pdf', '', 'pS6 Percent Positive')
	plot_bar_graph(plot_df, 'sgRNA_map', 'pAKT 308', 'Condition', './barplot_sgRNA_vs_pAKT_308.pdf', '','pAKT(308) MFI')
	plot_bar_graph(plot_df, 'sgRNA_map', 'pAKT 473', 'Condition', './barplot_sgRNA_vs_pAKT_473.pdf', '', 'pAKT(473) MFI')

	plot_scatterplot(plot_df_stim, 'LFC', 'pS6', 'p Value', [-0.2, 5], './scatterplot_LFC_vs_S6.pdf')
	plot_scatterplot(plot_df_stim, 'LFC', 'pAKT 308', 'p Value', [-0.2, 5], './scatterplot_LFC_vs_pAKT_308.pdf')
	plot_scatterplot(plot_df_stim, 'LFC', 'pAKT 473', 'p Value', [-0.2, 5], './scatterplot_LFC_vs_pAKT_473.pdf')

def calculate_Dunnett_from_df(in_df, col, cond_name = 'sgRNA_map', control_name = 'NT-sgRNA'):
	
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

def plot_bar_graph(df, x_val, y_val, hue_val, save_pn, xlabel, ylabel):

	fig, ax = plt.subplots(figsize=(18, 4))  # Adjust size for clarity
	sns.barplot(df, x=x_val, y=y_val, hue=hue_val, ax=ax, errorbar='sd', edgecolor='black', capsize=0.1)  # Add error bars and edges for better aesthetics
	sns.stripplot(data=df, x=x_val, y=y_val, hue=hue_val, dodge=True, ax=ax, marker="o", size=3, palette='dark:black', jitter=True, legend=False)
	ax.set_xlabel(xlabel, fontsize=14, labelpad=10)  # X-axis label
	ax.set_ylabel(ylabel, fontsize=14, labelpad=10)  # Y-axis label
	legend = ax.legend(title=hue_val, fontsize=12, title_fontsize=13, loc='center left', bbox_to_anchor=(1, 0.5), frameon=True)
	legend.get_frame().set_edgecolor('black')
	plt.xticks(rotation=45, horizontalalignment='right', fontsize=10)
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

	sgRNA_labels = df['sgRNA_map'].to_numpy().astype(str)
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
	# ax.set_ylabel('sgRNA', fontsize=12)

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
