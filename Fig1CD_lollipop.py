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

def main():

	# Load the Excel file
	file_path = '/Users/cf3041/Library/CloudStorage/OneDrive-ColumbiaUniversityIrvingMedicalCenter/Base Editing v3 - Clinical Assay Screens to pair with HACE (PIK3)/BEv3_008_pAKT_pS6_Screen_r3 (with MaxCyte electroporation)/Analysis with CJF/MAGeCK/DP_vs_DN_Maxcyte_X2.sgrna_summary.xlsx'
	data = pd.read_excel(file_path)
	data.head()  # View the structure of the data

	# Filter for PIK3CD and PIK3R1 and exclude rows with zero control counts
	filtered_data = data[(data['Gene'].isin(['PIK3CD', 'PIK3R1']))].copy()            

	# Exclude rows with zero control counts
	filtered_data['Control Counts'] = filtered_data['control_count'].str.split('/').apply(
	    lambda x: [float(count) for count in x if count.strip().isdigit()]
	)

	filtered_data = filtered_data[filtered_data['Control Counts'].apply(lambda x: all(c > 0 for c in x))]
	# edit_arr = filtered_data['Amino acid edits (3-9 window)'].to_numpy().astype(str)
	# gene_arr = filtered_data['Gene'].to_numpy().astype(str)
	splice_variant_arr = np.full(fill_value=False, dtype=bool, shape=filtered_data.shape[0])
	pos_arr = np.full(fill_value=-10, dtype=int, shape=filtered_data.shape[0])

	for curr_num, (_, curr_row) in enumerate(filtered_data.iterrows()):

		curr_gene = curr_row['Gene']
		curr_edit = curr_row['Amino acid edits (3-9 window)']
		curr_mut = curr_row['Mutation category']



		if not isinstance(curr_mut, str):
			if np.isnan(curr_mut):
				continue

		if 'Splice' in curr_mut:

			assert('Exon' in curr_edit)

			splice_variant_arr[curr_num] = True

			# exon_removed = re.sub(r'Exon\d+:\+?\d+;?', '', curr_edit)
	
			exon_num = int(curr_edit.split('Exon')[1].split(':')[0])
			exon_sign = curr_edit.split('Exon')[1].split(':')[1][0]

			if exon_sign == '-':
				exon_pos = 0
			elif exon_sign == '+':
				exon_pos = 1
			else:
				bp()

			if curr_gene == 'PIK3CD':
				pos_arr[curr_num] = PIK3CD_exon_mapping[exon_num][exon_pos]
			elif curr_gene == 'PIK3R1':
				pos_arr[curr_num] = PIK3R1_exon_mapping[exon_num][exon_pos]
			else:
				bp()

			continue

			# print(curr_edit + ' identified exon #' + '%d'%exon_num)

		all_positions = np.array(re.findall(r'\d+', curr_edit)).astype(int)

		if all_positions.size == 0:
			pass

		else:

			if 'Exon' in curr_edit:
				assert('Intron' in curr_mut)
			
			else:
				pos_arr[curr_num] = all_positions.mean()
				

	filtered_data['Amino Acid Position'] = pos_arr
	filtered_data['Splice Variant'] = splice_variant_arr
	# filtered_data['Significant'] = (filtered_data['LFC'] > 1.8) | (filtered_data['LFC'] < -3)
	filtered_data['Significant'] = filtered_data['p.twosided'] < 0.05

	filtered_data_pos = filtered_data[filtered_data['Amino Acid Position'] > 0]

	# Generate plots for PIK3CD and PIK3R1
	plot_gene_with_translucent_colors('PIK3CD', filtered_data_pos, PIK3CD_domains, 'PIK3CD_lollipop_plot_translucent_colors.pdf')
	plot_gene_with_translucent_colors('PIK3R1', filtered_data_pos, PIK3R1_domains, 'PIK3R1_lollipop_plot_fleshed_out_domains.pdf')


# Function to generate lollipop plots with annotated domains
def plot_gene_with_translucent_colors(gene, data, domains, file_name):

	subset = data[data['Gene'] == gene]
	paler_colors_translucent = [
	    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', 
	    '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5'
	]

	fig, ax = plt.subplots(figsize=(14, 6))
	ax.vlines(subset['Amino Acid Position'], 0, subset['LFC'], color='gray', alpha=0.6)

	# ax.scatter(subset['Amino Acid Position'], subset['LFC'], 
	#            s=abs(subset['LFC']) * 10, 
	#            c=subset['Significant'].map({True: 'red', False: 'gray'}), 
	#            alpha=0.8, zorder=3, marker='o')

	plt.scatter(subset[subset['Splice Variant'] == True]['Amino Acid Position'], 
            subset[subset['Splice Variant'] == True]['LFC'],
            s=abs(subset[subset['Splice Variant'] == True]['LFC']) * 10, 
            c=subset[subset['Splice Variant'] == True]['Significant'].map({True: 'red', False: 'gray'}), 
            alpha=0.8, zorder=3, marker='x')

	# Plot points where Splice Variant is False with 'o' marker
	plt.scatter(subset[subset['Splice Variant'] == False]['Amino Acid Position'], 
	            subset[subset['Splice Variant'] == False]['LFC'],
	            s=abs(subset[subset['Splice Variant'] == False]['LFC']) * 10, 
	            c=subset[subset['Splice Variant'] == False]['Significant'].map({True: 'red', False: 'gray'}), 
	            alpha=0.8, zorder=3, marker='o')


	for i, (start, end, domain_name) in enumerate(domains.get(gene, [])):
		ax.axvspan(start, end, color=paler_colors_translucent[i % len(paler_colors_translucent)], alpha=0.2, label=domain_name)

	ax.axhline(0, color='black', linewidth=0.8)
	ax.set_title(f'Lollipop Plot of Missense Variants for {gene} with Domains')
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