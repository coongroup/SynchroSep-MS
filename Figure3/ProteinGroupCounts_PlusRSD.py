# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 14:48:59 2025

@author: nlancaster
"""
 
# Figure 3 Protein Group Identifications
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
    return report

report = diann_filter(report)



#%%
def pg_count(df,file):
    
    df = df[df['Run']==file]
    pg_count = len(df['Protein.Group'].unique())
    return pg_count


dil_points = [10,10,10,50,50,50,100,100,100,200,200,200]

       

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



import numpy as np


files = ['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3',
        '[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3',
        '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3',
        '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']




def cum_count(report, run):
    run_report = report[report['Run']==run]
    
    import numpy as np
    num_points = 1000
    rt_bins = np.linspace(0,50,num_points)
    cum_count = []
    for rt in rt_bins:
        #loop through each point and filter for unique Protein Groups
        cum_count.append(len(run_report[run_report['RT']<=rt]['Protein.Group'].unique()))
    return rt_bins, cum_count

rt_bins, qc_1 = cum_count(report,files[0])
rt_bins, qc_2 = cum_count(report,files[1])
rt_bins, dual_1 = cum_count(report,files[2])
rt_bins, dual_2 = cum_count(report,files[3])

total = list(np.array(dual_1) + np.array(dual_2))
#%%
#count pg detections for single dual comparison
#store results for average, and error bars
single_col_results = [[],[],[]]

dual_col_results = [[],[],[]]



temp_count = []


#loop through files an average every three files and give max/min
for file in single_col_comp:
    count = pg_count(report,file)
    temp_count.append(count)
    
    
    
    if len(temp_count)==3:
        single_col_results[0].append(np.average(temp_count))
        single_col_results[1].append(max(temp_count) - np.average(temp_count))
        single_col_results[2].append(np.average(temp_count)-min(temp_count))
        temp_count = []


for file in dual_col_comp:
    count = pg_count(report,file)
    temp_count.append(count)
    
    
    
    if len(temp_count)==3:
        dual_col_results[0].append(np.average(temp_count))
        dual_col_results[1].append(max(temp_count) - np.average(temp_count))
        dual_col_results[2].append(np.average(temp_count)-min(temp_count))
        temp_count = []


#%% assemble protein RSDs

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
# dual_comp_1_report = diann_filter('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_TriplicateSearches/Rd07_Dual_Col1_Triplicate_LegacyNorm/report.parquet')
# dual_comp_2_report = diann_filter('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_TriplicateSearches/Rd07_Dual_Col2_Triplicate_LegacyNorm/report.parquet')

single_comp_1_rsd = rsd_list(single_comp_1_report,single_col_comp_1)
single_comp_2_rsd = rsd_list(single_comp_2_report,single_col_comp_2)
dual_comp_1_rsd = rsd_list(dual_comp_1_report,dual_col_comp_1)
dual_comp_2_rsd = rsd_list(dual_comp_2_report,dual_col_comp_2)
#%%
import matplotlib.pyplot as plt
import numpy as np
fig, (ax1,ax2,ax3) = plt.subplots(3,1,figsize = (3.5,6))




ax1.spines[['top','right']].set_visible(False)


x = np.arange(2)

width = 0.2
colors = ['dimgrey','cornflowerblue','silver','lightsteelblue']
capsize = 5
ax1.bar(x - width *(1/2),single_col_results[0],width = width,yerr = (single_col_results[1],single_col_results[2]),label = 'Single\nInjection',color = colors[0],capsize = capsize)
ax1.bar(x + width *(1/2),dual_col_results[0],width = width,yerr = (dual_col_results[1],dual_col_results[2]),label = 'Dual\nInjection',color = colors[2],capsize = capsize)
ax1.set_xticks(x,['Column 1','Column 2'])
ax1.tick_params('both',labelsize = 10)
ax1.set_ylabel('Protein Groups',fontsize = 10)
ax1.legend(frameon = False, fontsize = 9,ncol = 1,loc = 'upper center', handletextpad = 0.4,
           columnspacing = 0.5, handlelength = 1, labelspacing = 0.4,alignment = 'center'
)


ax3.spines[['top','right']].set_visible(False)


#cumulative protein group identifications

ax2.spines[['top','right']].set_visible(False)



ax2.plot(rt_bins,qc_1,label = 'Column 1 - Single Inj',color = 'black',linestyle = '--')
ax2.plot(rt_bins,qc_2,label = 'Column 2 - Single Inj',color = 'grey',linestyle = '--')
ax2.plot(rt_bins,total,label = 'Both - Dual Inj',color = 'maroon',linestyle = '--')
ax2.legend(frameon = False,labelspacing = 0.2,handletextpad = 0.3,handlelength = 1.5,
           loc = [0.35,0])
ax2.set_xlim(0,)
ax2.set_ylim(0,)
ax2.set_xlabel('Retention Time (min)')
ax2.set_ylabel('Cumulative Protein Groups')


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


for i,vals in enumerate(data):
    ax3.annotate(str(round(np.median(vals),2))+'%\nn='+str(len(vals)),(x[i],150),ha = 'center',va = 'bottom',fontsize = 8)



# ,bodycolor = "#2FA8DF"
ax3.set_xticks(x,['Single\nInjection\nColumn 1','Dual\nInjection\nColumn 1','Single\nInjection\nColumn 2','Dual\nInjection\nColumn 2'])

ax3.set_ylabel('%RSD (n=3)')
ax3.set_ylim(0,180)
ax3.axhline(20,linestyle = '--',color = 'black')



plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/20250611_Rd07_PG_Ids_v5.svg')