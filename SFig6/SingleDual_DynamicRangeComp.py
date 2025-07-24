# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 08:50:03 2025

@author: nlancaster
"""

# Comparison of files 
import pandas as pd
import numpy as np

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
    report['PG.MaxLFQ'] = report['PG.MaxLFQ'].replace([np.inf,0,-np.inf],np.nan)
    return report

report = diann_filter(report)

report = report.drop_duplicates(['Run','Protein.Group'])

#perform a pivot table to go to the wide format
pg_pivot = pd.pivot(report,index=['Protein.Group'] 
                             ,columns = 'Run',values = ['PG.MaxLFQ'])

#%%
import matplotlib.pyplot as plt


def prep_report(df,file_dual,file_single):
    import numpy as np
    df_single = df[df['Run']==file_single]
    df_dual = df[df['Run']==file_dual]
    
    in_dual = []
    for val in df_single['Protein.Group']:
        if val in list(df_dual['Protein.Group']):
            in_dual.append('True')
        else:
            in_dual.append('False')
            
    df_single['In Dual'] = in_dual
    df_single = df_single.sort_values(by = 'PG.MaxLFQ',ascending = False).reset_index(drop = True)
    df_single['PG.MaxLFQ'] = [np.log2(x) for x in df_single['PG.MaxLFQ']]
    return df_single

col1 = prep_report(report,'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3','[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3')
col2 = prep_report(report,'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3','[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3')
#%%


fig,(ax1,ax2) = plt.subplots(1,2,figsize = (7.5,3),layout = 'constrained')

alpha = 0.08
markersize = 15
ax1.scatter(col1.index[col1['In Dual']== 'True'], col1[col1['In Dual']== 'True']['PG.MaxLFQ'],alpha = alpha,edgecolor = 'None',
            label = 'In Dual Injection (n = ' + str(len(col1.index[col1['In Dual']== 'True'])) + ')',color = 'silver',s = markersize,rasterized=True)
ax1.scatter(col1.index[col1['In Dual']== 'False'], col1[col1['In Dual']== 'False']['PG.MaxLFQ'],alpha = alpha,edgecolor = 'None',
            label = 'Not In Dual Injection (n = ' + str(len(col1.index[col1['In Dual']== 'False'])) + ')',color = 'red',s = markersize,rasterized=True)
ax1.set_title('Column 1 - Single Injection')

ax2.scatter(col2.index[col2['In Dual']== 'True'], col2[col2['In Dual']== 'True']['PG.MaxLFQ'],alpha = alpha,edgecolor = 'None',
            label = 'In Dual Injection (n = ' + str(len(col2.index[col2['In Dual']== 'True'])) + ')',color = 'silver',s = markersize,rasterized=True)
ax2.scatter(col2.index[col2['In Dual']== 'False'], col2[col2['In Dual']== 'False']['PG.MaxLFQ'],alpha = alpha,edgecolor = 'None',
            label = 'Not In Dual Injection (n = ' + str(len(col2.index[col2['In Dual']== 'False'])) + ')',color = 'red',s = markersize,rasterized=True)
ax2.set_title('Column 2 - Single Injection')

axes = (ax1,ax2)
for ax in axes:
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlabel('Protein Rank')
    ax.set_ylabel('Log2 PG MaxLFQ')
    ax.legend(frameon = False)
    leg = ax.legend(frameon = False)   
    for lh in leg.legend_handles: 
        lh.set_alpha(0.8)

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'

fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/Revision_Round01/PythonFigures/20250717_Rd07_SingleDual_Intensity_OverlapIDs.svg',dpi = 1000)

