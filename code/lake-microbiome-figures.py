import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import operator
import plotly.express as px
import seaborn as sns
from matplotlib.patches import Patch  # Import Patch for custom legend

# Setting the output folder for the figures
output_folder = Path("OutputFigures")
output_folder.mkdir(exist_ok=True)  # Create the folder if it does not exist


data = pd.read_csv("coverm_431_MAGS_metagenomes_reads_count.csv")
data=data.T
data.columns = data.iloc[0]

# drop the first row
data= data.iloc[1:]
data

#sort the data by the samples
data_sorted = data.sort_index(ascending=True)


#Figure 1
fig, ax = plt.subplots(figsize=(10, 10.5))

# Iterate over each row in the DataFrame and plot the boxplot
for i, row in enumerate(data_sorted.iterrows(), 1):
    positions = [i]
    widths = [0.6]
    ax.boxplot(row[1], positions=positions, widths=widths)

# Set the x-axis labels using the DataFrame index
labels = data_sorted.index
ax.set_xticks(range(1, len(labels) + 1))
ax.set_xticklabels(labels, rotation=45, ha='right')

# Set the plot
plt.ylabel('The Relative Abundance of Genomes')
#have more outlier in the middle season
plt.savefig(output_folder / 'Boxplot_distribution_b.png', dpi=300)

# Normalization
total_row=data.sum(axis=1)
data = data.div(total_row, axis=0)

data_sorted = data.sort_index(ascending=True)


#Figure 2
fig, ax = plt.subplots(figsize=(10, 10.5))

# Iterate over each row in the DataFrame and plot the boxplot
for i, row in enumerate(data_sorted.iterrows(), 1):
    positions = [i]
    widths = [0.6]
    ax.boxplot(row[1], positions=positions, widths=widths)

# Set the x-axis labels using the DataFrame index
labels = data_sorted.index
ax.set_xticks(range(1, len(labels) + 1))
ax.set_xticklabels(labels, rotation=45, ha='right')

# Set the plot
plt.ylabel('The Relative Abundance of Genomes')
#have more outlier in the middle season
plt.savefig(output_folder / 'Boxplot_distribution.png', dpi=300)


#Figure 3
#for each month, count the total
july=data_sorted.iloc[0:3].sum()
aug=data_sorted.iloc[3:9].sum()
sep=data_sorted.iloc[9:11].sum()
october=data_sorted.iloc[11:16].sum()

# fig, ax = plt.subplots()
fig, ax = plt.subplots(figsize=(10,5))


ax.boxplot(july)
ax.boxplot(aug,positions=[2], widths=0.6)
ax.boxplot(sep,positions=[3], widths=0.6)
ax.boxplot(october,positions=[4], widths=0.6)

labels = ['July(3 samples)','August(6 samples)','September(2 samples)','October(5 samples)']
ax.set_xticklabels(labels, fontsize=16)  # Increase font size for x-axis labels

ax.set_ylabel('The relative abundance of genomes', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=13)  # Increase font size for tick labels
plt.savefig(output_folder /'month_abundance.png', dpi=300)
#august has more outliers than all the other months

#sorted genome each month by descedning order
sorted_ser7 = july.sort_values(ascending=False)
sorted_ser8 = aug.sort_values(ascending=False)
sorted_ser9 = sep.sort_values(ascending=False)
sorted_ser10 = october.sort_values(ascending=False)

different_month=[]
dd11 =sorted_ser7[0:10]
dd21 =sorted_ser8[0:10]
dd31 =sorted_ser9[0:10]
dd41 =sorted_ser10[0:10]


fig, (ax1_, ax2_,ax3_,ax4_) = plt.subplots(4, 1,figsize=(20, 20))
ax1_.pie(dd11,autopct='%1.2f%%')
for i in dd11.index:
    different_month.append(i)
ax1_.legend(dd11.index,title = "types of genome:",bbox_to_anchor=(1.1, 1.05))

# ax2.pie(sd2)
ax2_.pie(dd21,autopct='%1.2f%%')
for i in dd21.index:
    different_month.append(i)
ax2_.legend(dd21.index,title = "types of genome:",bbox_to_anchor=(1.1, 1.05))

# ax3.pie(sd3)
ax3_.pie(dd31,autopct='%1.2f%%')
for i in dd31.index:
    different_month.append(i)
ax3_.legend(dd31.index,title = "types of genome:",bbox_to_anchor=(1.1, 1.05))

ax4_.pie(dd41,autopct='%1.2f%%')
for i in dd41.index:
    different_month.append(i)

# dd11.name = 'July'
# df_top10_month=dd11.to_frame()
# dd21.name = 'August'
# df_top10_month=df_top10_month.join(dd21, how='outer')
# dd31.name = 'September'
# df_top10_month=df_top10_month.join(dd31, how='outer')
# dd41.name = 'October'
# df_top10_month=df_top10_month.join(dd41, how='outer')
# df_top10_month= df_top10_month.astype(np.float32)

# #Figure 4

# fig, ax = plt.subplots(figsize=(12, 10))  # Adjusted figure size to give more space

# df_top10_month_binary = df_top10_month.fillna(0)
# df_binary = df_top10_month_binary.where(df_top10_month_binary == 0, other=1)
# df_binary = df_binary.loc[df_binary.sum(axis=1).sort_values(ascending=False).index]
# df_binary = df_binary.T

# # Create a heatmap with binary values
# sns.heatmap(df_binary, annot=True, cmap='coolwarm', vmin=0, vmax=1,
#             cbar=False, square=True, linewidths=0.5, linecolor='gray')

# # Custom legend
# legend_handles = [Patch(facecolor='blue', label='0: does not exist in top10'),
#                   Patch(facecolor='red', label='1: exists in top10')]
# ax.legend(handles=legend_handles, title='Legend', loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

# plt.subplots_adjust(right=0.8)  # Adjust subplot parameters to give more space on the right for the legend
# plt.tight_layout()  # Optional: can use plt.tight_layout() to adjust layouts automatically

# # Remove the x-axis title and adjust labels
# ax.set_xlabel('')
# ax.set_ylabel('Months')

# # Save and show the plot
# plt.savefig(output_folder / 'Binary Map of genome vs. months.png', dpi=300, bbox_inches='tight')

# Load the phylum data
phylum_data = pd.read_csv('MAG_taxonomy.tsv', sep='\t')[['Genome', 'Phylum']]

# Assume 'july', 'aug', 'sep', 'october' are Series objects with genomes as indices and values are some form of abundance
# Replace these with your actual Series data
sorted_ser7 = july.sort_values(ascending=False).head(10)
sorted_ser8 = aug.sort_values(ascending=False).head(10)
sorted_ser9 = sep.sort_values(ascending=False).head(10)
sorted_ser10 = october.sort_values(ascending=False).head(10)

# Combine data into one DataFrame
df = pd.concat([sorted_ser7, sorted_ser8, sorted_ser9, sorted_ser10], axis=1)
df.columns = ['July', 'August', 'September', 'October']
df = df.join(phylum_data.set_index('Genome'), how='left')

# Fill NaN with 0 for binary heatmap
df_binary = df[['July', 'August', 'September', 'October']].notnull().astype(int)

# Sort the dataframe by Phylum
df['Phylum'] = df.index.map(phylum_data.set_index('Genome')['Phylum'])
df_sorted = df.sort_values(by=['Phylum', 'July', 'August', 'September', 'October'], ascending=[True, False, False, False, False])
df_binary_sorted = df_binary.loc[df_sorted.index]

# Create the heatmap with annotations
fig, ax = plt.subplots(figsize=(15, 12))
sns.heatmap(df_binary_sorted, cmap='coolwarm', cbar=False, linewidths=0.5, linecolor='gray', annot=True, ax=ax)

# Move the x-axis (months) to the top
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')

# Add horizontal lines to group the taxonomic groups
unique_phyla = df_sorted['Phylum'].unique()
previous_phylum = None
for idx, genome in enumerate(df_sorted.index):
    current_phylum = df_sorted.loc[genome, 'Phylum']
    if previous_phylum and previous_phylum != current_phylum:
        ax.axhline(y=idx, color='black', linestyle='--')
    previous_phylum = current_phylum

# Custom legend with descriptive labels
legend_handles = [Patch(facecolor='blue', label='Not in Top 10'),
                  Patch(facecolor='red', label='In Top 10')]
ax.legend(handles=legend_handles, title='Presence in Top 10', loc='upper left', bbox_to_anchor=(1.05, 1))

plt.tight_layout()
plt.savefig(output_folder / 'Binary Map of genome vs. months.svg', format='svg', dpi=300, bbox_inches='tight')
## change in Inkscape



#Figure 5
#organized them based on depth 5,10,15,23.5
five=data_sorted.loc[['2020-07-24_5m','2020-10-19_5m']].sum()
ten=data_sorted.loc[['2020-08-25_10m','2020-08-05_10m']].sum()
fifteen=data_sorted.loc[['2020-07-24_15m','2020-08-05_15m','2020-08-25_15m','2020-09-11_15m','2020-10-08_15m','2020-10-19_15m']].sum()
twenth2_5=data_sorted.loc[['2020-07-24_23.5m','2020-08-05_23.5m','2020-08-25_23.5m','2020-09-11_23.5m','2020-10-08_23.5m','2020-10-19_23.5m']].sum()

#fig, ax = plt.subplots()
fig, ax = plt.subplots(figsize=(10,5))
ax.boxplot(five)
ax.boxplot(ten,positions=[2], widths=0.6)
ax.boxplot(fifteen,positions=[3], widths=0.6)
ax.boxplot(twenth2_5,positions=[4], widths=0.6)

labels = ['5m(2 samples)','10m(2 samples)','15m(6 samples)','23.5m(6 samples)']
ax.set_xticklabels(labels, fontsize=16)  # Increase font size for x-axis labels

ax.set_ylabel('The relative abundance of genomes', fontsize=16)
plt.tick_params(axis='both', which='major', labelsize=13)

plt.savefig(output_folder /'depth_abundance.png', dpi=300)

sorted_five = five.sort_values(ascending=False)
sorted_ten = ten.sort_values(ascending=False)
sorted_fifteen = fifteen.sort_values(ascending=False)
sorted_23 = twenth2_5.sort_values(ascending=False)

different_depth=[]
five_ =sorted_five[0:10]
_ten =sorted_ten[0:10]
_fifteen =sorted_fifteen[0:10]
_23 =sorted_23[0:10]


fig, (ax1_, ax2_,ax3_,ax4_) = plt.subplots(4, 1,figsize=(20, 20))
ax1_.pie(five_,autopct='%1.2f%%')
for i in five_.index:
    different_depth.append(i)
ax1_.legend(five_.index,title = "types of genome:",bbox_to_anchor=(1.1, 1.05))

# ax2.pie(sd2)
ax2_.pie(_ten,autopct='%1.2f%%')
for i in _ten.index:
    different_depth.append(i)
ax2_.legend(_ten.index,title = "types of genome:",bbox_to_anchor=(1.1, 1.05))

# ax3.pie(sd3)
ax3_.pie(_fifteen,autopct='%1.2f%%')
for i in _fifteen.index:
    different_depth.append(i)
ax3_.legend(_fifteen.index,title = "types of genome:",bbox_to_anchor=(1.1, 1.05))

ax4_.pie(_23,autopct='%1.2f%%')
for i in _23.index:
    different_depth.append(i)

five_.name = 'depth5m'
df_top10=five_.to_frame()
_ten.name = 'depth10m'
df_top10=df_top10.join(_ten, how='outer')
_fifteen.name = 'depth15m'
df_top10=df_top10.join(_fifteen, how='outer')
_23.name = 'depth23.5m'
df_top10=df_top10.join(_23, how='outer')
df_top10= df_top10.astype(np.float32)

#Figure 6

# fig, ax = plt.subplots(figsize=(12, 10))  # Adjusted figure size for better layout

# # Assuming df_top10 is defined somewhere in your code, process it:
# df_top10_depth_binary = df_top10.fillna(0)
# df_depth_binary = df_top10_depth_binary.where(df_top10_depth_binary == 0, other=1)
# df_depth_binary = df_depth_binary.T  # Transpose the DataFrame

# # Create a heatmap with binary values
# sns.heatmap(df_depth_binary, annot=True, cmap='coolwarm', vmin=0, vmax=1,
#             cbar=False, square=True, linewidths=0.5, linecolor='gray')

# # Remove x-axis title
# ax.set_xlabel('')

# # Custom legend
# legend_handles = [Patch(facecolor='blue', label='0: does not exist in top10'),
#                   Patch(facecolor='red', label='1: exists in top10')]
# ax.legend(handles=legend_handles, title='Legend', loc='upper left', bbox_to_anchor=(1.05, 1), borderaxespad=0.)

# # Add axis labels and title
# ax.set_ylabel('Depths')

# # Save and show the plot
# plt.subplots_adjust(right=0.85)  # Adjust right margin to make room for the legend
# plt.savefig(output_folder / 'binarymap_depth.png', dpi=300, bbox_inches='tight')

# Join phylum data with top10 data
df_top10 = df_top10.join(phylum_data.set_index('Genome'), how='left')

# Fill NaN with 0 for binary heatmap
df_binary = df_top10.notnull().astype(int)

# Sort the dataframe by Phylum
df_top10['Phylum'] = df_top10.index.map(phylum_data.set_index('Genome')['Phylum'])
df_sorted = df_top10.sort_values(by=['Phylum', 'depth5m', 'depth10m', 'depth15m', 'depth23.5m'], ascending=[True, False, False, False, False])
df_binary_sorted = df_binary.loc[df_sorted.index]

# Create the heatmap with annotations
fig, ax = plt.subplots(figsize=(15, 12))
sns.heatmap(df_binary_sorted.iloc[:, :-1], cmap='coolwarm', cbar=False, linewidths=0.5, linecolor='gray', annot=True, ax=ax)


# Move the x-axis (depths) to the top
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')

# Add horizontal lines and annotations to group the taxonomic groups by phylum
unique_phyla = df_sorted['Phylum'].unique()
previous_phylum = None
for idx, genome in enumerate(df_sorted.index):
    current_phylum = df_sorted.loc[genome, 'Phylum']
    if previous_phylum and previous_phylum != current_phylum:
        ax.axhline(y=idx, color='black', linestyle='--')
    previous_phylum = current_phylum

# Annotate the last group
if previous_phylum:
    ax.text(-0.5, len(df_sorted) - 0.5, previous_phylum, verticalalignment='center', fontsize=10, color='black')

# Custom legend with descriptive labels
legend_handles = [Patch(facecolor='blue', label='Not In top10'),
                  Patch(facecolor='red', label='In top10')]
ax.legend(handles=legend_handles, title='Presence in Top 10', loc='upper left', bbox_to_anchor=(1.05, 1))

plt.tight_layout()
plt.savefig(output_folder /'binarymap_depth.svg', format='svg', dpi=300, bbox_inches='tight')

#Figure 7

from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.max_display_count = 10000  # Set the max display count to a high number
count_top10_month={}
for i in different_month:
    if i not in count_top10_month:
        count_top10_month[i]=1
    else:
        count_top10_month[i]+=1
count_top10_month=dict(sorted(count_top10_month.items(), key=operator.itemgetter(1),reverse=True))

#save it as fraction
import operator
from fractions import Fraction

count_top10_month = {}
for i in different_month:
    if i not in count_top10_month:
        count_top10_month[i] = 1
    else:
        count_top10_month[i] += 1

count_top10_month = dict(sorted(count_top10_month.items(), key=operator.itemgetter(1), reverse=True))

# Divide each value by 4 and convert to fractions
for key in count_top10_month:
    value = count_top10_month[key]
    fraction = Fraction(value, 4)
    count_top10_month[key] = str(fraction)
    
    
#here we are only focusing on top 10 for each depth
count_top10_depth={}
for i in different_depth:
    if i not in count_top10_depth:
        count_top10_depth[i]=1
    else:
        count_top10_depth[i]+=1
count_top10_depth=dict(sorted(count_top10_depth.items(), key=operator.itemgetter(1),reverse=True))

count_top10_depth = {}
for i in different_depth:
    if i not in count_top10_depth:
        count_top10_depth[i] = 1
    else:
        count_top10_depth[i] += 1

count_top10_depth = dict(sorted(count_top10_depth.items(), key=operator.itemgetter(1), reverse=True))

# Divide each value by 4 and convert to fractions
for key in count_top10_depth:
    value = count_top10_depth[key]
    fraction = Fraction(value, 4)
    count_top10_depth[key] = str(fraction)
    
all_dic = {}

# Convert the values in count_top10_month to fractions and add to the 'all_dic' dictionary
for i in count_top10_month:
    if i not in all_dic:
        fraction_str = count_top10_month[i]
        fraction = Fraction(fraction_str)
        all_dic[i] = fraction
    else:
        fraction_str = count_top10_month[i]
        fraction = Fraction(fraction_str)
        all_dic[i] += fraction

# Convert the values in count_top10_depth to fractions and add to the 'all_dic' dictionary
for j in count_top10_depth:
    if j not in all_dic:
        fraction_str = count_top10_depth[j]
        fraction = Fraction(fraction_str)
        all_dic[j] = fraction
    else:
        fraction_str = count_top10_depth[j]
        fraction = Fraction(fraction_str)
        all_dic[j] += fraction


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Load the taxonomy data from the TSV file
file_path = 'MAG_taxonomy.tsv'
taxonomy_data = pd.read_csv(file_path, sep='\t')

# Create a dictionary mapping genomes to their phylum
genome_phylum_dict = pd.Series(taxonomy_data.Phylum.values, index=taxonomy_data.Genome).to_dict()

# Assuming all_dic is a dictionary defined elsewhere in your code with genome IDs and some numeric values
all_dic = dict(sorted(all_dic.items(), key=lambda item: item[1], reverse=True))

# Convert the dictionary to a DataFrame
data = {'Genome': list(all_dic.keys()), 'Frequency': [float(value) for value in all_dic.values()]}
df = pd.DataFrame(data)

# Sort the DataFrame by 'Frequency' in descending order
df.sort_values(by='Frequency', ascending=False, inplace=True)

# Incorporate phylum information
df['Phylum'] = df['Genome'].map(genome_phylum_dict)

# Create the horizontal bar plot
plt.figure(figsize=(25, 15))
bar_plot = sns.barplot(x='Frequency', y='Genome', data=df)

# Set the plot title and axis labels with increased font sizes
plt.ylabel("Genome", fontsize=24)  # Increase font size for y-axis title
plt.yticks(fontsize=24)  # Increase font size for y-axis tick labels
plt.xlabel("")  # Remove the x-axis label
plt.xticks(fontsize=24)  # Increase font size for x-axis tick labels

# Add labels for the phylum to the bars
for index, row in df.iterrows():
    bar_plot.text(row['Frequency'] + 0.01, index, str(row['Phylum']), color='black', ha="left", va='center', fontsize=24)

# Save the plot ensuring that all layout elements fit well
plt.tight_layout()  # Adjust subplots to give some padding and prevent cut-off
plt.savefig(output_folder /'Sum_ratio.png', dpi=300)


#Figure 8
#measured by Pearson’s correlation
fig, ax = plt.subplots(figsize=(15,13))
data_sorted_=data_sorted.T.astype(np.float32)
corr=data_sorted_.corr()
sns.heatmap(corr, annot=True, cmap='coolwarm')
plt.savefig(output_folder / 'correlation matrix_Pearson.png',dpi=300)

#Figure 9
#measured by spearman’s correlation
fig, ax = plt.subplots(figsize=(15,13))
corr_spearman=data_sorted_.corr(method='spearman')
sns.heatmap(corr_spearman, annot=True, cmap='coolwarm')
plt.savefig(output_folder /'correlation matrix_Spearman.png',dpi=300)












