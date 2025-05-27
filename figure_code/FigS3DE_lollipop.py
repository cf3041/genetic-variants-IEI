import numpy as np
import os, glob
from pdb import set_trace as bp
import matplotlib.pyplot as plt
import pandas as pd
import re

# PIK3CD exon annotations from: https://useast.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000171608;r=1:9629889-9729114;t=ENST00000377346
PIK3CD_exon_mapping = {
	3 : [1, 48],
	4 : [48, 125],
	5 : [125, 201],
	6 : [201, 261],
	7 : [261, 311],
	8 : [311, 341],
	9 : [341, 415],
	10 : [415, 448],
	11 : [448, 491],
	12 : [491, 508],
	13 : [508, 564],
	14 : [564, 605],
	15 : [605, 653],
	16 : [653, 686],
	17 : [686, 746],
	18 : [746, 784],
	19 : [784, 810],
	20 : [810, 866],
	21 : [866, 907],
	22 : [907, 956],
	23 : [956, 1000],
	24 : [1000, 1044]
}

# PIK3R1 exon annotations from https://useast.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000145675;r=5:68215740-68301821;t=ENST00000521381
PIK3R1_exon_mapping = {
	2 : [1, 113],
	3 : [113, 144],
	4 : [144, 169],
	5 : [169, 213],
	6 : [213, 280],
	7 : [280, 307],
	8 : [307, 341],
	9 : [341, 374],
	10 : [374, 434],
	11 : [434, 476],
	12 : [476, 524],
	13 : [524, 583],
	14 : [583, 606],
	15 : [606, 663],
	16 : [663, 724]
}

PIK3CD_domains = {
'PIK3CD': [
	(16, 105, 'PI3K-ABD (Adaptor-Binding Domain)'),
	(187, 278, 'PI3K-RBD (Ras-Binding Domain)'),
	(287, 312, 'Disordered Region'),
	(319, 476, 'C2 PI3K-type Domain'),
	(497, 674, 'PIK Helical Domain'),
	(745, 1027, 'PI3K/PI4K Catalytic Domain'),
	(751, 757, 'G-loop'),
	(890, 898, 'Catalytic Loop'),
	(909, 935, 'Activation Loop')
	]
}

PIK3R1_domains = {
'PIK3R1': [
	(3, 79, 'SH3 Domain'),
	(80, 108, 'Disordered Region'),
	(84, 98, 'Proline-Rich Region'),
	(113, 301, 'Rho-GAP Domain'),
	(333, 428, 'SH2 Domain 1'),
	(429, 623, 'Inter-SH2 Domain (iSH2)'),
	(624, 718, 'SH2 Domain 2')
	]
}

P_THRESH = 0.05

def main():

	# Load the Excel files
	AKT_data = pd.read_excel('./DP_vs_DN_Maxcyte_X2.sgrna_summary.xlsx')
	prolif_data = pd.read_excel('./restim_vs_unselected.sgrna_summary.xlsx')

	AKT_pos_arr = np.full(fill_value=-10, dtype=int, shape=AKT_data.shape[0])
	prolif_pos_arr = np.full(fill_value=-10, dtype=int, shape=prolif_data.shape[0])

	for curr_AKT_ind, curr_AKT_row in AKT_data.iterrows():

		curr_gene = curr_AKT_row['Gene']
		curr_edit = curr_AKT_row['Amino acid edits (3-9 window)']

		curr_prolif_ind = np.where(prolif_data['sgrna']==curr_AKT_row['sgrna'])[0]
		assert(curr_prolif_ind.size == 1)

		curr_prolif_row = prolif_data.iloc[curr_prolif_ind]

		curr_prolif_pval = curr_prolif_row['p.twosided'].to_numpy()[0]

		if curr_gene != 'PIK3CD' and curr_gene != 'PIK3R1': # Only plot PIK3CD and PIK3R1
			continue

		if not isinstance(curr_edit, str): # Ignore empty window edits
			continue


		if curr_AKT_row['p.twosided'] >= P_THRESH and curr_prolif_pval >= P_THRESH: # Only plot edits significant in one or the other screen
			continue

		if 0 in np.array(curr_AKT_row['control_count'].split('/')).astype(float): # Ignore any control count = 0 in pAKT data
			continue

		if 0 in np.array(curr_prolif_row['control_count'].to_numpy().astype(str)[0].split('/')).astype(float): # Ignore any control count = 0 in proliferation data
			continue


		if 'Exon' in curr_edit: # Splice edits

			exon_num = int(curr_edit.split('Exon')[1].split(':')[0])
			exon_sign = curr_edit.split('Exon')[1].split(':')[1][0]

			if exon_sign == '-':
				exon_pos = 0
			elif exon_sign == '+':
				exon_pos = 1
			else:
				bp()

			if curr_gene == 'PIK3CD':
				AKT_pos_arr[curr_AKT_ind] = PIK3CD_exon_mapping[exon_num][exon_pos]
				prolif_pos_arr[curr_prolif_ind] = PIK3CD_exon_mapping[exon_num][exon_pos]
			elif curr_gene == 'PIK3R1':
				AKT_pos_arr[curr_AKT_ind] = PIK3R1_exon_mapping[exon_num][exon_pos]
				prolif_pos_arr[curr_prolif_ind] = PIK3R1_exon_mapping[exon_num][exon_pos]
			else:
				bp()

		else: # Coding edits

			all_positions = np.array(re.findall(r'\d+', curr_edit)).astype(int)

			if all_positions.size == 0:
				pass

			else:
				AKT_pos_arr[curr_AKT_ind] = all_positions.mean()
				prolif_pos_arr[curr_prolif_ind] = all_positions.mean()

	AKT_data['Amino Acid Position'] = AKT_pos_arr
	prolif_data['Amino Acid Position'] = prolif_pos_arr

	filtered_AKT_data = AKT_data[AKT_data['Amino Acid Position'] > 0]
	filtered_prolif_data = prolif_data[prolif_data['Amino Acid Position'] > 0]

	# Generate plots for PIK3CD and PIK3R1
	plot_gene_with_translucent_colors('PIK3CD', filtered_AKT_data, filtered_prolif_data, PIK3CD_domains, 'PIK3CD_lollipop_AKT_vs_prolif.pdf')
	plot_gene_with_translucent_colors('PIK3R1', filtered_AKT_data, filtered_prolif_data, PIK3R1_domains, 'PIK3R1_lollipop_AKT_vs_prolif.pdf')


# Function to generate lollipop plots with annotated domains
def plot_gene_with_translucent_colors(gene, data1, data2, domains, file_name):

	subset1 = data1[data1['Gene'] == gene]
	subset2 = data2[data2['Gene'] == gene]
	paler_colors_translucent = [
	    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', 
	    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5'
	]

	fig, ax = plt.subplots(figsize=(14, 6))
	ax.vlines(subset1['Amino Acid Position'], 0, subset1['LFC'], color='gray', alpha=0.6)
	ax.vlines(subset2['Amino Acid Position'], 0, subset2['LFC'], color='gray', alpha=0.6)

	plt.scatter(subset1['Amino Acid Position'], subset1['LFC'], s=abs(subset1['LFC']) * 10, c='Blue', alpha=0.8, zorder=3)
	plt.scatter(subset2['Amino Acid Position'], subset2['LFC'], s=abs(subset2['LFC']) * 10, c='Red', alpha=0.8, zorder=3)

	for i, (start, end, domain_name) in enumerate(domains.get(gene, [])):
		ax.axvspan(start, end, color=paler_colors_translucent[i % len(paler_colors_translucent)], alpha=0.2, label=domain_name)

	ax.axhline(0, color='black', linewidth=0.8)
	ax.set_title(f'Lollipop Plot of Significant Variants for {gene}')
	ax.set_ylabel('Log Fold Change (LFC)')
	ax.set_xlabel('Amino Acid Position')
	ax.grid(alpha=0.3)
	handles, labels = ax.get_legend_handles_labels()
	by_label = dict(zip(labels, handles))
	ax.legend(by_label.values(), by_label.keys(), loc='center left', bbox_to_anchor=(1, 0.5), title='Protein Domains')
	plt.tight_layout()
	plt.savefig(file_name, bbox_inches='tight', format='pdf')
	plt.close()


if __name__=="__main__":
    main()
