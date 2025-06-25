# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 15:51:28 2025

@author: nlancaster
"""

#Figure 4 - Quant Comparisons - prep dilution series
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
    report['PG.MaxLFQ'] = report['PG.MaxLFQ'].replace([np.inf,0,-np.inf],np.nan)
    return report

report = diann_filter(report)

report = report.drop_duplicates(['Run','Protein.Group'])

#perform a pivot table to go to the wide format
pg_pivot = pd.pivot(report,index=['Protein.Group'] 
                             ,columns = 'Run',values = ['PG.MaxLFQ'])


col1_dil = ['[0.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_1',
'[0.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_2',
'[0.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_3',
'[0.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_1',
'[0.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_2',
'[0.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_3',
'[0.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_1',
'[0.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_2',
'[0.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_3',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']

col2_dil =['[4.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_1',
'[4.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_2',
'[4.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_3',
'[4.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_1',
'[4.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_2',
'[4.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_3',
'[4.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_1',
'[4.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_2',
'[4.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_3',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']

col1_inv_dil = ['[0.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_1',
'[0.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_2',
'[0.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_3',
'[0.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_1',
'[0.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_2',
'[0.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_3',
'[0.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_1',
'[0.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_2',
'[0.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_3',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_1',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_2',
'[0.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_3']

col2_inv_dil = ['[4.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_1',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_2',
'[4.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_3',
'[4.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_1',
'[4.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_2',
'[4.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_3',
'[4.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_1',
'[4.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_2',
'[4.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_3',
'[4.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_1',
'[4.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_2',
'[4.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_3']

inv_dil_points = [10,10,10,50,50,50,100,100,100,200,200,200]
#%% calculate r2 values
import scipy.stats as stats
import numpy as np

conc_dict  = {'[0.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_1':10,
                '[0.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_2':10,
                '[0.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_3':10,
                '[0.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_1':50,
                '[0.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_2':50,
                '[0.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_3':50,
                '[0.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_1':100,
                '[0.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_2':100,
                '[0.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_3':100,
                '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1':200,
                '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2':200,
                '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3':200,
                '[4.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_1':10,
                '[4.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_2':10,
                '[4.0]20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_3':10,
                '[4.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_1':50,
                '[4.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_2':50,
                '[4.0]20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_3':50,
                '[4.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_1':100,
                '[4.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_2':100,
                '[4.0]20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_3':100,
                '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1':200,
                '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2':200,
                '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3':200,
                '[0.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_1':10,
                '[0.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_2':10,
                '[0.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_3':10,
                '[0.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_1':50,
                '[0.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_2':50,
                '[0.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_3':50,
                '[0.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_1':100,
                '[0.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_2':100,
                '[0.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_3':100,
                '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_1':200,
                '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_2':200,
                '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_3':200,
                '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_1':10,
                '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_2':10,
                '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_3':10,
                '[4.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_1':50,
                '[4.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_2':50,
                '[4.0]20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_3':50,
                '[4.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_1':100,
                '[4.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_2':100,
                '[4.0]20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_3':100,
                '[4.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_1':200,
                '[4.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_2':200,
                '[4.0]20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_3':200}

def ignore_na_r2(x,y):
    df = pd.DataFrame({'x':x,'y':y})
    df = df.dropna(axis = 0,subset = ('y'))
    conc_points = [conc_dict[file] for file in df['x']]
    if len(set(conc_points))>=3:
        r2 = stats.pearsonr(conc_points, df['y'])[0]**2
        return r2
    return np.nan #this might be a bit aggressive, but seems like an ok way to treat bad values

index = 0
 

inv_r2_1 = []
dil_r2_1 = []
inv_r2_2 = []
dil_r2_2 = []



for index in range(0,len(pg_pivot.index.get_level_values('Protein.Group'))):
    print(index)

    
    data = pg_pivot['PG.MaxLFQ'][col1_dil].iloc[index]
    pearson_r1 = ignore_na_r2(data.index,data)
    dil_r2_1.append(pearson_r1)

  
    data = pg_pivot['PG.MaxLFQ'][col1_inv_dil].iloc[index]
    pearson_r2 = ignore_na_r2(data.index,data)
    inv_r2_1.append(pearson_r2)

    
    data = pg_pivot['PG.MaxLFQ'][col2_dil].iloc[index]
    pearson_r1 = ignore_na_r2(data.index,data)
    dil_r2_2.append(pearson_r1)

    
    
    data = pg_pivot['PG.MaxLFQ'][col2_inv_dil].iloc[index]
    pearson_r2 = ignore_na_r2(data.index,data)
    inv_r2_2.append(pearson_r2)

dil_r2_1 = pd.Series(dil_r2_1).dropna()
dil_r2_2 = pd.Series(dil_r2_2).dropna()
inv_r2_1 = pd.Series(inv_r2_1).dropna()
inv_r2_2 = pd.Series(inv_r2_2).dropna()
#save r2_calc
dil_r2_1.to_csv('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles/int_med_pg_dil_r2_1.csv',index = False)
dil_r2_2.to_csv('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles/int_med_pg_dil_r2_2.csv',index = False)
inv_r2_1.to_csv('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles/int_med_pg_inv_r2_1.csv',index = False)
inv_r2_2.to_csv('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles/int_med_pg_inv_r2_2.csv',index = False)

#%%
import numpy as np
def pg_count(df,file):
    
    df = df[df['Run']==file]
    pg_count = len(df['Protein.Group'].unique())
    return pg_count

#count prec detections for single dual comparison
#store results for average, and error bars
single_col_results = [[],[],[]]

dual_col_results = [[],[],[]]



temp_count = []


#loop through files an average every three files and give max/min


#count pg identifications for dilution series
#store results for average, and error bars
col1_dir = [[],[],[]]
col1_inv = [[],[],[]]
col2_dir = [[],[],[]]
col2_inv = [[],[],[]]



temp_count = []
#loop through files an average every three files and give max/min
for file in col1_dil:
    temp_count.append(pg_count(report,file))
    
    if len(temp_count)==3:
        col1_dir[0].append(np.average(temp_count))
        col1_dir[1].append(max(temp_count) - np.average(temp_count))
        col1_dir[2].append(np.average(temp_count)-min(temp_count))
        temp_count = []

for file in col2_dil:
    temp_count.append(pg_count(report,file))
    
    if len(temp_count)==3:
        col2_dir[0].append(np.average(temp_count))
        col2_dir[1].append(max(temp_count) - np.average(temp_count))
        col2_dir[2].append(np.average(temp_count)-min(temp_count))
        temp_count = []

for file in col1_inv_dil:
    temp_count.append(pg_count(report,file))
    
    if len(temp_count)==3:
        col1_inv[0].append(np.average(temp_count))
        col1_inv[1].append(max(temp_count) - np.average(temp_count))
        col1_inv[2].append(np.average(temp_count)-min(temp_count))
        temp_count = []

for file in col2_inv_dil:
    temp_count.append(pg_count(report,file))
    
    if len(temp_count)==3:
        col2_inv[0].append(np.average(temp_count))
        col2_inv[1].append(max(temp_count) - np.average(temp_count))
        
        col2_inv[2].append(np.average(temp_count)-min(temp_count))
        temp_count = []
        
#note that I am flipping the order of the column 2 data to plot the data more intuitively
col2_inv[0].reverse()
col2_inv[1].reverse()
col2_inv[2].reverse()
#%%read in R2 plots
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
def read_file(path):
    df = pd.read_csv(path)
    return [float(x) for x in df['0']]

dil_r2_1 = read_file('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles/int_med_pg_dil_r2_1.csv')
dil_r2_2 = read_file('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles/int_med_pg_dil_r2_2.csv')
inv_r2_1 = read_file('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles/int_med_pg_inv_r2_1.csv')
inv_r2_2 = read_file('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles/int_med_pg_inv_r2_2.csv')


# Make plots
import matplotlib.pyplot as plt
import numpy as np
import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap

fig = plt.figure(figsize = (8.5,6))
gs = gridspec.GridSpec(2, 2, height_ratios=[1,1],hspace = 0.45, wspace = 0.4)


#plot PG IDs


#adjust the relative sizes of the top row plots
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
gs_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0, :], width_ratios=[2, 1.85],wspace = 0.2)
ax1.set_position(gs_top[0].get_position(fig))
ax2.set_position(gs_top[1].get_position(fig))
ax1.set_subplotspec(gs_top[0])
ax2.set_subplotspec(gs_top[1])



divider = make_axes_locatable(ax1)
ax1Left = divider.append_axes('left',size = '100%',pad = 0.8)
ax1Right = ax1


ax1Left.spines[['top','right']].set_visible(False)
ax1Right.spines[['top','right']].set_visible(False)


x = np.arange(4)
capsize = 1
width = 0.4
colors = ['dimgrey','silver',]
ax1Left.bar(x - width *(1/2),col1_dir[0],width = width,yerr = (col1_dir[1],col1_dir[2]),label = 'Column 1',color = colors[0],capsize = capsize)
ax1Left.bar(x + width *(1/2),col2_dir[0],width = width,yerr = (col2_dir[1],col2_dir[2]),label = 'Column 2',color = colors[1],capsize = capsize)
ax1Right.bar(x - width *(1/2),col1_inv[0],width = width,yerr = (col1_inv[1],col1_inv[2]),label = 'Column 1',color = colors[0],capsize = capsize)
ax1Right.bar(x + width *(1/2),col2_inv[0],width = width,yerr = (col2_inv[1],col2_inv[2]),label = 'Column 2',color = colors[1],capsize = capsize)



ax1Left.set_title('Direct',fontsize = 10)
ax1Left.set_xticks(x,['10\n10','50\n50','100\n100','200\n200'],linespacing = 1.3 )
ax1Left.set_xlabel('Loading Mass (ng)',fontsize = 10)
ax1Left.tick_params('both',labelsize = 10)
ax1Left.set_ylabel('Protein Groups',fontsize = 10)
ax1Left.legend(frameon = True,fontsize = 10,ncol = 1, loc = [0.02,0.01],handletextpad = 0.2,
           columnspacing = 0.5, handlelength = 1, labelspacing = 0.2,borderpad = 0.1)



ax1Right.set_title('Inverted',fontsize = 10)
ax1Right.set_xticks(x,['10\n200','50\n100','100\n50','200\n10'],linespacing = 1.3)
ax1Right.set_xlabel('Loading Mass (ng)',fontsize = 10)
ax1Right.tick_params('both',labelsize = 10)
ax1Right.set_ylabel('Protein Groups',fontsize = 10)
ax1Right.legend(frameon = True,fontsize = 10,ncol = 1, loc = [0.02,0.01],handletextpad = 0.2,
           columnspacing = 0.5, handlelength = 1, labelspacing = 0.2,borderpad = 0.1)


ax1Left.tick_params(axis = 'both',labelsize = 10)
ax1Right.tick_params(axis = 'both',labelsize = 10)
ax1Left.set_ylim(0,9300)
ax1Right.set_ylim(0,9300)



#make violin plot of R2 values
x = np.arange(4)

width = 0.2

colors = ['dimgrey','cornflowerblue','silver','lightsteelblue']

data = [dil_r2_1,dil_r2_2,inv_r2_1,inv_r2_2]


labels = ['Column 1\nDirect','Column 2\nDirect','Column 1\nInverted','Column 2\nInverted']







medians = [str(round(np.median(x),1)) for x in data]
edge_color = 'black'
body_color= '#2CA8E0'
violin1 = ax2.violinplot(data,bw_method = 0.3,showextrema = False,showmedians = True)
for violin in violin1['bodies']:
    violin.set_facecolor(body_color)
    violin.set_edgecolor(edge_color)
    violin.set_linewidth(0)
    violin.set_alpha(1)
violin1['cmedians'].set_linewidth(1)
violin1['cmedians'].set_color('black')

pos = [1,2,3,4]
for i,vals in enumerate(data):
    ax2.annotate(str(round(np.median(vals),3))+'\nn='+str(len(vals)),(pos[i],1.01),ha = 'center',va = 'bottom',fontsize = 8)

ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0],[0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax2.set_ylim(0,1.1)
ax2.spines[['top','right']].set_visible(False)
plt.setp(ax2.collections, edgecolor="k")
ax2.set_xticks(pos,labels,rotation = 0,fontsize = 10)
ax2.set_ylabel('R\u00B2',fontsize = 10)

ax2.tick_params('both',labelsize = 10)

ax2.axhline(0.8,linestyle = '--',color = 'black',linewidth = 1)




def prep_data(x,y):
    #clean out na values
    x = pd.Series(x)
    y = pd.Series(y)
    int_df  = pd.DataFrame({'x':list(x),'y':list(y)})
    int_df = int_df.dropna()
    int_df = int_df[int_df['x']!=0]
    int_df = int_df[int_df['y']!=0]
    
    x = list(int_df['x'])
    y = list(int_df['y'])
    x = [np.log2(val) for val in x]
    y = [np.log2(val) for val in y]
    return x,y
pg_report = report.drop_duplicates(['Run','Protein.Group'])

pg_pivot = pd.pivot(pg_report,index=['Protein.Group'] 
                             ,columns = 'Run',values = ['PG.MaxLFQ'])



#loop through files to get more data points

file_list = ['20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_1',
'20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_2',
'20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_3',
'20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_1',
'20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_2',
'20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_3',
'20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1',
'20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2',
'20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3',
'20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_1',
'20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_2',
'20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_3']

x_col_comp = []
y_col_comp = []
for file in file_list:
    x_temp, y_temp = prep_data(pg_pivot['PG.MaxLFQ']['[0.0]' + file],
          pg_pivot['PG.MaxLFQ']['[4.0]' + file])
    x_col_comp = x_col_comp + x_temp
    y_col_comp = y_col_comp + y_temp



file_list = ['20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_1',
'20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_2',
'20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_3',
'20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_1',
'20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_2',
'20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_3',
'20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_1',
'20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_2',
'20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_3',
'20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_1',
'20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_2',
'20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_3',
'20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_1',
'20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_2',
'20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_3',
'20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1',
'20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2',
'20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3',
'20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_1',
'20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_2',
'20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_3',
'20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_1',
'20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_2',
'20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_3']




#sum up the results across files - Column 1

x_ser_comp = []
y_ser_comp = []

file_list_direct = ['20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_1',
'20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_2',
'20241216_Rd07_CM_50ngMouse_Zoe_50ngMouse_3',
'20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_1',
'20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_2',
'20241216_Rd07_CM_100ngMouse_Zoe_100ngMouse_3',
'20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1',
'20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_2',
'20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3',
'20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_1',
'20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_2',
'20241216_Rd07_CM_10ngMouse_Zoe_10ngMouse_3']

file_list_inverse = [
'20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_1',
'20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_2',
'20241216_Rd07_CM_50ngMouse_Zoe_100ngMouse_3',
'20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_1',
'20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_2',
'20241216_Rd07_CM_100ngMouse_Zoe_50ngMouse_3',
'20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_1',
'20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_2',
'20241216_Rd07_CM_200ngMouse_Zoe_10ngMouse_3',
'20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_1',
'20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_2',
'20241216_Rd07_CM_10ngMouse_Zoe_200ngMouse_3']

for i,file in enumerate(file_list_direct):
    x_temp, y_temp = prep_data(pg_pivot['PG.MaxLFQ']['[0.0]' + file],
          pg_pivot['PG.MaxLFQ']['[0.0]' + file_list_inverse[i]])
    x_ser_comp = x_ser_comp + x_temp
    y_ser_comp = y_ser_comp + y_temp


def lin_fit(xi,x,y):
    from scipy import stats
    lin_reg = stats.linregress(x,y)
    
    m = lin_reg[0]
    b = lin_reg[1]
    r2 = lin_reg[2]**2
    
    yi = m*xi + b
    
    return yi,m,b,r2





# "Viridis-like" colormap with white background
white_viridis = LinearSegmentedColormap.from_list('white_viridis', [
    (0, '#ffffff'),
    (1e-20, '#440053'),
    (0.2, '#404388'),
    (0.4, '#2a788e'),
    (0.6, '#21a784'),
    (0.8, '#78d151'),
    (1, '#fde624'),
], N=256)

ax3 = fig.add_subplot(2, 2, 3, projection='scatter_density')
density = ax3.scatter_density(x_col_comp, y_col_comp, cmap=white_viridis)
fig.colorbar(density, label='Number of Points')






ax3.plot([5,43],[lin_fit(5,x_col_comp,y_col_comp)[0],lin_fit(43,x_col_comp,y_col_comp)[0]],linestyle = '--',color = 'black',linewidth = 1)
ax3.annotate('y=' + str(round(lin_fit(8,x_col_comp,y_col_comp)[1],3)) + 'x+' + str(round(lin_fit(8,x_col_comp,y_col_comp)[2],3)) +
             '\nR\u00B2=' + str(round(lin_fit(8,x_col_comp,y_col_comp)[3],3)),(1,0),fontsize = 10,ha = 'right',va = 'bottom',xycoords = 'axes fraction')



ax3.spines[['top','right']].set_visible(False)


ax3.tick_params('both',labelsize = 10)

ax3.set_title('Column Comparison',fontsize = 10,fontweight = 'bold')


ax3.set_xlabel('Log2 (Column 1 Quantity)',fontsize = 10)
ax3.set_ylabel('Log2 (Column 2 Quantity)',fontsize = 10)
ax3.set_xlim([5,43])
ax3.set_ylim([5,43])


ax4 = fig.add_subplot(2,2,4,projection = 'scatter_density')


ax4.spines[['top','right']].set_visible(False)
density = ax4.scatter_density(x_ser_comp, y_ser_comp, cmap=white_viridis)
fig.colorbar(density, label='Number of Points')
ax4.tick_params('both',labelsize = 10)
ax4.set_title('Dilution Series Comparison',fontsize = 10,fontweight = 'bold')
ax4.plot([5,43],[lin_fit(5,x_ser_comp, y_ser_comp)[0],lin_fit(43,x_ser_comp, y_ser_comp)[0]],linestyle = '--',color = 'black',linewidth = 1)
ax4.annotate('y=' + str(round(lin_fit(8,x_ser_comp, y_ser_comp)[1],3)) + 'x+' + str(round(lin_fit(8,x_ser_comp, y_ser_comp)[2],3)) +
             '\nR\u00B2=' + str(round(lin_fit(8,x_ser_comp, y_ser_comp)[3],3)),(1,0),fontsize = 10,ha = 'right',va = 'bottom',xycoords = 'axes fraction')
ax4.set_xlim([5,43])
ax4.set_ylim([5,43])
ax4.set_xlabel('Log2 (Direct Series Quantity)',fontsize = 10)
ax4.set_ylabel('Log2 (Inverted Series Quantity)',fontsize = 10)




plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'

fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/20250612_Rd07_PGQuantPerformance_v3.svg',dpi = 1000)


# ADD MULTIPLE FILES FOR THE SERIES COMPARISON


