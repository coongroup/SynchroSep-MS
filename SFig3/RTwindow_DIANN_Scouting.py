# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 14:23:39 2025

@author: nlancaster
"""

#DIANN report parsing - window scouting
import pandas as pd


def diann_filter(report):
    
    filter_columns = ['Q.Value','Lib.Q.Value','Global.Q.Value',
                     'Lib.PG.Q.Value','Global.PG.Q.Value',
                     'Lib.Peptidoform.Q.Value','Global.Peptidoform.Q.Value',
                     'PG.Q.Value']
    for col in filter_columns:
        report = report[report[col]<=0.01]
    
    report = report.reset_index(drop = True)
    return report

def lin_fit(xi,x,y):
    from scipy import stats
    #clean out na values
    x = pd.Series(x)
    y = pd.Series(x)
    drop_x = x.isnull()
    x = x[drop_x==False]
    y = y[drop_x==False]
    
    drop_y = y.isnull()
    x = x[drop_y==False]
    y = y[drop_y==False]
    
    lin_reg = stats.linregress(x,y)
    
    m = lin_reg[0]
    b = lin_reg[1]
    r2 = lin_reg[2]**2
    
    yi = m*xi + b
    
    return yi,m,b,r2

#%% Prepare scouting data


window_widths = [1.5,1.75,2,2.25,2.5,2.75,3,3.25,3.5]



error_col1 = []

error_col2 = []

#these two lists will store the total number of precursors identified in the combined search
comb_col1 = []
comb_col2 = []

scouting_reports = ['P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap1p5_RTShift4/report.parquet',
                    'P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap1p75_RTShift4/report.parquet',
                    'P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap2_RTShift4/report.parquet',
                    'P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap2p25_RTShift4/report.parquet',
                    'P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap2p5_RTShift4/report.parquet',
                    'P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap2p75_RTShift4/report.parquet',
                    'P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap3_RTShift4/report.parquet',
                    'P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap3p25_RTShift4/report.parquet',
                    'P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap3p5_RTShift4/report.parquet'
                    ]

for file in scouting_reports:
    report = pd.read_parquet(file)
    report = diann_filter(report)
    prec_pivot = pd.pivot(report,index=['Precursor.Id', 'Stripped.Sequence', 'Precursor.Charge','Precursor.Mz', 'Modified.Sequence'] 
                             ,columns = 'Run',values = ['Precursor.Quantity','RT'])


    #make plots of column assignment accuracy
    column1_results = report[report['Run'].str.contains('[4.0]')]
    column2_results = report[report['Run'].str.contains('[4.0]')]
    qc_col1_rt1 = len(column1_results[column1_results['Run'] == '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3']['Precursor.Id'].unique())
    qc_col1_rt2 = len(column2_results[column2_results['Run'] == '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3']['Precursor.Id'].unique())
    qc_col2_rt1 = len(column1_results[column1_results['Run'] == '[0.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3']['Precursor.Id'].unique())
    qc_col2_rt2 = len(column2_results[column2_results['Run'] == '[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3']['Precursor.Id'].unique())
    comb_col1_rt1 = len(column1_results[column1_results['Run'] == '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']['Precursor.Id'].unique())
    comb_col2_rt2 = len(column2_results[column1_results['Run'] == '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']['Precursor.Id'].unique())
    
    temp_col2_error = qc_col1_rt2/(qc_col1_rt1+qc_col1_rt2) * 100
    temp_col1_error = qc_col2_rt1/(qc_col2_rt1+qc_col2_rt2) * 100
    
    error_col1.append(temp_col1_error)
    error_col2.append(temp_col2_error)
    comb_col1.append(comb_col1_rt1)
    comb_col2.append(comb_col2_rt2)
#%% Bar Chart and Scatter Plot generation
import matplotlib.pyplot as plt
import numpy as np


layout = [
    ['T','T','T','U','U','U'],
    ['T','T','T','U','U','U'],
    ['A','A','A','A','B','B'],
    ['A','A','A','A','C','C']]
fig, axd = plt.subplot_mosaic(layout, figsize=(8,8))

ax1 = axd['A']
ax2 = axd['B']
ax3 = axd['C']
ax4 = axd['T']
ax5 = axd['U']

#Column 1 error rate/precursors
ax4.spines[['top','right']].set_visible(False)
ax4_2 = ax4.twinx()
ax4_2.spines[['top','left']].set_visible(False)
ax4_2.spines['right'].set_color('maroon')
ax4_2.tick_params(axis='y', colors='maroon')

ax4.plot(window_widths,error_col1,marker = 's',color = 'black')
ax4_2.plot(window_widths,comb_col1,marker = 'o',color = 'maroon',linestyle = '--')
ax4.set_title('Column 1 RT Search')

ax4.set_ylabel('Column Assignment Error Rate (%)')
ax4.set_ylim(0,)
ax4_2.set_ylabel('Precursors',color = 'maroon')
ax4_2.set_ylim(0,160000)
ax4.set_xlabel('RT Window Width (min)')
ax4.axvline(2,linestyle = '--',color = 'black',linewidth = 1.5)




#Column 2 Error rate/precursors
ax5.spines[['top','right']].set_visible(False)
ax5_2 = ax5.twinx()
ax5_2.spines[['top','left']].set_visible(False)
ax5_2.spines['right'].set_color('maroon')
ax5_2.tick_params(axis='y', colors='maroon')

ax5.plot(window_widths,error_col2,marker = 's',color = 'black' )
ax5_2.plot(window_widths,comb_col2,marker = 'o',color = 'maroon',linestyle = '--')





ax5.set_ylabel('Column Assignment Error Rate (%)')
ax5.set_ylim(0,)
ax5_2.set_ylabel('Precursors',color = 'maroon')
ax5_2.set_ylim(0,130000)
ax5.set_xlabel('RT Window Width (min)')
ax5.set_title('Column 2 RT Search')
ax5.axvline(2,linestyle = '--',color = 'black',linewidth = 1.5)

#Example of 3.5 min window search showing the analysis picking up on the other column
axes = (ax1,ax2,ax3)
for ax in axes:
    ax.spines[['top','right']].set_visible(False)


labels = ['Column 1\nInjection\nOnly','Column 2\nInjection\nOnly','Combined\nInjection']
x = np.arange(3)
ax1.set_ylabel('Precursors')
ax1.set_xticks(x,labels)

axes = (ax2,ax3)
for ax in axes:
    ax.set_xlabel('Single Injection Ref RT (min)',fontsize = 10)
    ax.set_ylabel('Dual Injection\nRT (min)',fontsize = 10)

width = 0.45

#3.5 min window
report = pd.read_parquet('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_ModDIANN_DualCol_WindowScouting_v1/ModDIANN_RTGap3p5_RTShift4/report.parquet')
report = diann_filter(report)
prec_pivot = pd.pivot(report,index=['Precursor.Id', 'Stripped.Sequence', 'Precursor.Charge','Precursor.Mz', 'Modified.Sequence'] 
                             ,columns = 'Run',values = ['Precursor.Quantity','RT'])

#make plots of column assignment accuracy
column1_results = report[report['Run'].str.contains('[4.0]')]
column2_results = report[report['Run'].str.contains('[4.0]')]
qc_col1_rt1 = len(column1_results[column1_results['Run'] == '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3']['Precursor.Id'].unique())
qc_col1_rt2 = len(column2_results[column2_results['Run'] == '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3']['Precursor.Id'].unique())
qc_col2_rt1 = len(column1_results[column1_results['Run'] == '[0.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3']['Precursor.Id'].unique())
qc_col2_rt2 = len(column2_results[column2_results['Run'] == '[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3']['Precursor.Id'].unique())
comb_col1_rt1 = len(column1_results[column1_results['Run'] == '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']['Precursor.Id'].unique())
comb_col2_rt2 = len(column2_results[column1_results['Run'] == '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']['Precursor.Id'].unique())


col1_data = [prec_pivot['RT']['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3'].tolist(),
            prec_pivot['RT']['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3'].tolist()]

col2_data = [prec_pivot['RT']['[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3'].tolist(),
            prec_pivot['RT']['[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3'].tolist()]

ax1.bar(x-width/2,[qc_col1_rt1,qc_col2_rt1,comb_col1_rt1],width = width, label = 'Column 1',color = '#6A6A6A',align = 'center' )
ax1.bar(x+width/2,[qc_col1_rt2,qc_col2_rt2,comb_col2_rt2],width = width, label = 'Column 2', color = '#AAAAAA',align = 'center' )
ax1.set_xticks(x,labels,rotation = 0)
ax1.set_xticks(x,labels,rotation = 0)
ax1.set_ylabel('Precursors',fontsize = 10)
ax1.tick_params('both',labelsize = 10)
ax1.set_title('RT Window Width = 3.5 min')

ylim_max = 200000
ax1.set_ylim(0,ylim_max)

for i,val in enumerate([qc_col1_rt1,qc_col2_rt1,comb_col1_rt1]):
    ax1.annotate(str(val),(x[i]-width/2,val+40),ha = 'center',fontsize = 10)


#add labels and lines to the bar chart 
ymin = qc_col1_rt2
ymax = qc_col1_rt1
ax1.axvline(0.25,ymin = ymin/ylim_max, ymax= ymax/ylim_max,linestyle = '--',color = 'black',linewidth = 1)
ax1.plot((0.05,0.45),(ymax,ymax),color = 'black',linewidth = 1)
pcnt = str(round(ymin/(ymin + ymax)*100,1))
ax1.annotate(pcnt + "%",(0.25,(ymax-ymin)/2+ymin),va = 'center',ha = 'center',fontsize = 10)

ax1.ticklabel_format(axis = 'y',style = 'scientific',scilimits = (0,1))

for i,val in enumerate([qc_col1_rt2,qc_col2_rt2,comb_col2_rt2]):
    t = ax1.annotate(str(val),(x[i]+width/2,val+40),ha = 'center',fontsize = 10)


#add labels and lines to the bar chart 
ymin = qc_col2_rt1
ymax = qc_col2_rt2
ax1.axvline(0.75,ymin = ymin/ylim_max, ymax= ymax/ylim_max,linestyle = '--',color = 'black',linewidth = 1)
ax1.plot((0.55,0.95),(ymax,ymax),color = 'black',linewidth = 1)
pcnt = str(round(ymin/(ymin + ymax)*100,1))
ax1.annotate(pcnt + "%",(0.75,(ymax-ymin)/2 + ymin),va = 'center',ha = 'center',fontsize = 10)

ax1.legend(title = 'Predicted RTs',frameon = False,fontsize = 10,title_fontsize = 10,
           handletextpad = 0.2,labelspacing=0.2,handlelength = 1, loc ='upper right',
           ncol = 2, columnspacing = 0.5)



#make RT scatter plots

col1_data = [prec_pivot['RT']['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3'].tolist(),
            prec_pivot['RT']['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3'].tolist()]

col2_data = [prec_pivot['RT']['[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3'].tolist(),
            prec_pivot['RT']['[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3'].tolist()]


ax2.plot([5,43],[lin_fit(5,col1_data[0],col1_data[1])[0],lin_fit(43,col1_data[0],col1_data[1])[0]],linestyle = '--',color = 'black')
ax2.annotate('y=' + str(lin_fit(5,col1_data[0],col1_data[1])[1]) + 'x+' + str(lin_fit(5,col1_data[0],col1_data[1])[2]) +
             '\nR\u00B2=' + str(lin_fit(5,col1_data[0],col1_data[1])[3]),(1,0),fontsize = 10,ha = 'right',va = 'bottom',xycoords = 'axes fraction')

ax3.plot([5,43],[lin_fit(5,col2_data[0],col2_data[1])[0],lin_fit(43,col2_data[0],col2_data[1])[0]],linestyle = '--',color = 'black')
ax3.annotate('y=' + str(lin_fit(8,col2_data[0],col2_data[1])[1]) + 'x+' + str(lin_fit(8,col2_data[0],col2_data[1])[2]) +
             '\nR\u00B2=' + str(lin_fit(8,col2_data[0],col2_data[1])[3]),(1,0),fontsize = 10,ha = 'right',va = 'bottom',xycoords = 'axes fraction')



ax2.spines[['top','right']].set_visible(False)
ax2.scatter(col1_data[0],col1_data[1],marker = '.',color = '#2CA8E0', linewidth = 0.5,s = 15,rasterized=True)

ax3.spines[['top','right']].set_visible(False)
ax3.scatter(col2_data[0],col2_data[1],marker = '.',color = '#2CA8E0', linewidth = 0.5, s = 15,rasterized=True)

ax2.tick_params('both',labelsize = 10)
ax3.tick_params('both',labelsize = 10)
ax2.set_title('Column 1',fontsize = 12)
ax3.set_title('Column 2',fontsize = 12)


ax2.set_xlabel('Single Injection Ref RT (min)',fontsize = 10)
ax2.set_ylabel('Dual Injection\nRT (min)',fontsize = 10)

ax3.set_xlabel('Single Injection Ref RT (min)',fontsize = 10)
ax3.set_ylabel('Dual Injection\nRT (min)',fontsize = 10)

ax2.set_xlim([5,43])
ax2.set_ylim([5,43])
ax3.set_xlim([5,43])
ax3.set_ylim([5,43])



plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout(pad = 0.4)
fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/20250610_Rd07_PrecursorRTAssignment_DIANN_MethodScout_v1.svg',dpi = 1000)
