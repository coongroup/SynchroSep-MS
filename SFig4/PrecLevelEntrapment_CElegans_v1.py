# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 08:46:38 2025

@author: nlancaster
"""

# Precursor level FDR Assessment
def diann_filter(report_path):
    import pandas as pd
    report = pd.read_parquet(report_path)
    filter_columns = ['Q.Value']
    for col in filter_columns:
        report = report[report[col]<=0.01]
    
    report = report.drop_duplicates(['Run','Precursor.Id'])
    report = report.reset_index(drop = True)
    
    return report

#read in report
import pandas as pd

prec_report  = pd.read_parquet('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_CElegans_EntrapmentSearch_v2/PrecLevel_EntrapmentSearch_v2/report.parquet')


def fdp(target,entrap):
    fdp_calc = entrap/(target + entrap)*100
    return fdp_calc
#%%
import numpy as np
file_names = []
fdp_values = []



df = prec_report

for file in df['Run'].unique():
    #iterate through each run in search
    file_names.append(file)
    temp_df = df[df['Run']==file].reset_index(drop = True)
    
    #initiate 
    target_count = 0
    entrap_count = 0
           
    for pg in temp_df['Protein.Names']:

        if '_MOUSE' in pg:
            target_count = target_count + 1
        
        elif '_CAEEL' in pg:
            entrap_count = entrap_count +1

    fdp_values.append(fdp(target_count,entrap_count))
result_df = pd.DataFrame({'Run':file_names,'FDP':fdp_values})

#%% Prepare Bar chart of error rates
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize =( 6,4))
ax.spines[['top','right']].set_visible(False)


x = [1,2,3,5,6,7,9,10,11,13,14,15]
ax.bar(x,result_df['FDP'],color = '#2FA8DF')

for i,val in enumerate(result_df['FDP']):
    ax.annotate(str(round(val,2)),(x[i],val+ 0.04),rotation = 90,ha = 'center', va = 'bottom')


ax.set_xticks([1,2,3,5,6,7,9,10,11,13,14,15],[str(x) for x in [1,2,3,1,2,3,1,2,3,1,2,3]])
ax.set_ylim(0,2)
ax.set_ylabel('False Discovery Proportion (%)')
ax.axhline(1,linestyle = '--',color = 'black',linewidth = 1.5)

labels = ['__________\nSingle Inj\nColumn 1','__________\nDual Inj\nColumn 1','__________\nSingle Inj\nColumn 2','__________\nDual Inj\nColumn 2']
positions = [2,6,10,14]
for i,label in enumerate(labels):
    ax.text(positions[i],-0.36,label,ha = 'center')

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/Revision_Round01/PythonFigures/20250717_Rd07_CElegansEntrapment_PrecLevel_v1.svg',dpi = 1000)
