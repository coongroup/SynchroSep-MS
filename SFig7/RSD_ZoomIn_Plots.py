# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 14:48:59 2025

@author: nlancaster
"""
 
# Zoom In on RSD values

import pandas as pd


single_col_comp = ['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC2',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC4',
'[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC2',
'[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3',
'[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC4'
]


dual_col_comp = ['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']






# assemble protein RSDs

def diann_filter(report_path):
    import pandas as pd
    report = pd.read_parquet(report_path)
    filter_columns = ['Q.Value','Lib.Q.Value','Global.Q.Value',
                     'Lib.PG.Q.Value','Global.PG.Q.Value',
                     'Lib.Peptidoform.Q.Value','Global.Peptidoform.Q.Value',
                     'PG.Q.Value']
    for col in filter_columns:
        report = report[report[col]<=0.01]
    
    report = report.reset_index(drop = True)
    report = report.drop_duplicates(subset = ['Run','Protein.Group'])
    return report


def rsd_calc(input_list):
    std_dev = np.std(input_list)
    average = np.average(input_list)
    if average == 0:
        return np.nan
    return std_dev/average * 100
    

def rsd_list(report,files):
    report = report[report['Run'].isin(files)]
    rsd_vals = []
    for row in report.groupby('Protein.Group'):
        if len(row[1]['Protein.Group'])==3:
            rsd_value = rsd_calc(row[1]['PG.MaxLFQ'])
            rsd_vals.append(rsd_value)
    rsd_vals = pd.Series(rsd_vals)
    rsd_vals = rsd_vals.dropna()
    
    return rsd_vals

import numpy as np


single_col_comp_1 = ['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC2',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC4']

single_col_comp_2 = ['[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC2',
'[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3',
'[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC4']

dual_col_comp_1 = ['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']

dual_col_comp_2 = ['[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']




single_comp_1_report = diann_filter('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_TriplicateSearches/Rd07_Single_Col1_TriplicateSearch/report.parquet')
single_comp_2_report = diann_filter('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_TriplicateSearches/Rd07_Single_Col2_TriplicateSearch/report.parquet')
dual_comp_1_report = diann_filter('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_TriplicateSearches/Rd07_Dual_Col1_TriplicateSearch/report.parquet')
dual_comp_2_report = diann_filter('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_TriplicateSearches/Rd07_Dual_Col2_TriplicateSearch/report.parquet')

single_comp_1_rsd = rsd_list(single_comp_1_report,single_col_comp_1)
single_comp_2_rsd = rsd_list(single_comp_2_report,single_col_comp_2)
dual_comp_1_rsd = rsd_list(dual_comp_1_report,dual_col_comp_1)
dual_comp_2_rsd = rsd_list(dual_comp_2_report,dual_col_comp_2)
#%%
import matplotlib.pyplot as plt
import numpy as np
fig, ax3 = plt.subplots(figsize = (3.5,2))



#Quant RSD values
ax3.spines[['top','right']].set_visible(False)
x = np.arange(4)

data = [single_comp_1_rsd,dual_comp_1_rsd,single_comp_2_rsd,dual_comp_2_rsd]


medians = [str(round(np.median(x),1)) for x in data]
edge_color = 'black'
body_color= '#2CA8E0'
violin1 = ax3.violinplot(data,x,bw_method = 0.3,showextrema = False,showmedians = True)
for violin in violin1['bodies']:
    violin.set_facecolor(body_color)
    violin.set_edgecolor(edge_color)
    violin.set_linewidth(0)
    violin.set_alpha(1)
violin1['cmedians'].set_linewidth(1)
violin1['cmedians'].set_color('black')




# ,bodycolor = "#2FA8DF"
ax3.set_xticks(x,['Single\nInjection\nColumn 1','Dual\nInjection\nColumn 1','Single\nInjection\nColumn 2','Dual\nInjection\nColumn 2'])

ax3.set_ylabel('%RSD (n=3)')
ax3.set_ylim(0,25)




plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/Revision_Round01/PythonFigures/20250717_Rd07_ProteinRSD_ZoomIn.svg',dpi = 1000)