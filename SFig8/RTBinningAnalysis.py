# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 08:59:19 2025

@author: nlancaster
"""

# Binning Missing IDs by retention time
import pandas as pd

report = pd.read_parquet('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualColSearch_FullSearch_NoNorm/ModDIANN_RTGap2_RTShift4_FullSearch_NoNorm/report.parquet')

def diann_filter(report):
    
    filter_columns = ['Q.Value','Lib.Q.Value','Global.Q.Value',
                     'Lib.PG.Q.Value','Global.PG.Q.Value',
                     'Lib.Peptidoform.Q.Value','Global.Peptidoform.Q.Value',
                     'PG.Q.Value']
    for col in filter_columns:
        report = report[report[col]<=0.01]
    
    report = report.reset_index(drop = True)
    
    #in the DIA-NN report, I found some -inf and 0 values in the PG.MaxLFQ column

    return report

report = diann_filter(report)

report = report.drop_duplicates(['Run','Precursor.Id'])





#%%
def prep_report(df,file_dual,file_single):
    df_single = df[df['Run']==file_single]
    df_dual = df[df['Run']==file_dual]
    


            
    df_single['In Dual'] = df_single['Precursor.Id'].isin(df_dual['Precursor.Id'])

    return df_single

col1 = prep_report(report,'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3','[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3')
col2 = prep_report(report,'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3','[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3')
#%%
rt_bins = [*range(0,50,4 )]
rt_labels = [] #for bar chart
col1_missing_fraction = []
col2_missing_fraction = []
col1_total = []
col2_total = []

for i,x in enumerate(rt_bins):
    if i+1 == len(rt_bins):
        break
    
    rt_labels.append((rt_bins[i+1]-rt_bins[i])/2 + rt_bins[i])
    
    temp_df = col1[(col1['RT']>rt_bins[i])&(col1['RT']<=rt_bins[i+1])]

    
    if False in temp_df['In Dual'].value_counts().index:
        missing_count = temp_df['In Dual'].value_counts().loc[False]
    else:
        missing_count = 0
    if True in temp_df['In Dual'].value_counts().index:
        present_count = temp_df['In Dual'].value_counts().loc[True]
    else:
        present_count = 0
    
    if present_count == 0 and missing_count == 0:
        col1_missing_fraction.append(0)
    else:
        col1_missing_fraction.append(missing_count/(missing_count+ present_count))
    col1_total.append(present_count)
    
    temp_df = col2[(col2['RT']>rt_bins[i])&(col2['RT']<=rt_bins[i+1])]
    
    if False in temp_df['In Dual'].value_counts().index:
        missing_count = temp_df['In Dual'].value_counts().loc[False]
    else:
        missing_count = 0
    if True in temp_df['In Dual'].value_counts().index:
        present_count = temp_df['In Dual'].value_counts().loc[True]
    else:
        present_count = 0
    
    if present_count == 0 and missing_count == 0:
        col2_missing_fraction.append(0)
    else:
        col2_missing_fraction.append(missing_count/(missing_count+ present_count))
    col2_total.append(present_count)
    

#%%
import matplotlib.pyplot as plt
import numpy as np
fig, (ax1,ax2) = plt.subplots(2,1,figsize = (4,6),layout = 'constrained')
ax1.spines[['top','right']].set_visible(False)
x = np.arange(len(rt_labels))
width = 0.3

ax1.bar(x-width/2,col1_missing_fraction,width = width,label = 'Column 1',color = 'dimgrey')
ax1.bar(x+width/2,col2_missing_fraction,width = width,label = 'Column 2',color = 'silver')

ax1.set_xticks(x,[str(int(x)) for x in rt_labels])
ax1.legend(frameon = False)
ax1.set_xlabel('RT Bin Center (min)')
ax1.set_ylabel('Fraction of Missing Precursors')


ax2.spines[['top','right']].set_visible(False)
ax2.plot(pd.Series(col1_total).index,col1_total,marker = 's',color = 'dimgrey',label = 'Column 1')
ax2.plot(pd.Series(col2_total).index,col2_total,marker = 's',color = 'silver',label = 'Column 2')
ax2.set_xticks(x,[str(int(x)) for x in rt_labels])
ax2.set_xlabel('RT Bin Center (min)')
ax2.set_ylabel('Number of Precursors')
ax2.legend(frameon = False)
ax2.set_ylim(0,)

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'

fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/Revision_Round01/PythonFigures/20250717_Rd07_SingleDualComparison_RTbinning_PrecursorIdentifications.svg',dpi = 1000)
