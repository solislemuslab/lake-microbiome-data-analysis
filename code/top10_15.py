import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import operator
import os

data = pd.read_csv("coverm_431_MAGS_metagenomes_reads_count.csv")
data=data.T
data.columns = data.iloc[0]

# drop the first row
data= data.iloc[1:]

total_row=data.sum(axis=1)
data = data.div(total_row, axis=0)

data_sorted = data.sort_index(ascending=True)

five=data_sorted.loc[['2020-07-24_5m','2020-10-19_5m']].sum()
ten=data_sorted.loc[['2020-08-25_10m','2020-08-05_10m']].sum()
fifteen=data_sorted.loc[['2020-07-24_15m','2020-08-05_15m','2020-08-25_15m','2020-09-11_15m','2020-10-08_15m','2020-10-19_15m']].sum()
twenth2_5=data_sorted.loc[['2020-07-24_23.5m','2020-08-05_23.5m','2020-08-25_23.5m','2020-09-11_23.5m','2020-10-08_23.5m','2020-10-19_23.5m']].sum()

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

july=data_sorted.iloc[0:3].sum()
aug=data_sorted.iloc[3:9].sum()
sep=data_sorted.iloc[9:11].sum()
october=data_sorted.iloc[11:16].sum()


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


# Extract the first 10 and 15 values from 'Genome'
top_10_genomes = df['Genome'].head(10)
top_15_genomes = df['Genome'].head(15)

# Directory path
directory = 'feature_selection'

# Create directory if it doesn't exist
if not os.path.exists(directory):
    os.makedirs(directory)

# Save these DataFrames as CSV files in the directory
top_10_genomes.to_csv(os.path.join(directory, 'top_10_genomes.csv'), index=False)
top_15_genomes.to_csv(os.path.join(directory, 'top_15_genomes.csv'), index=False)