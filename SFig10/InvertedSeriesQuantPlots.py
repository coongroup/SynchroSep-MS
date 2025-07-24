# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 09:36:06 2025

@author: nlancaster
"""

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


#perform a pivot table to go to the wide format
prec_pivot = pd.pivot(report,index=['Precursor.Id', 'Stripped.Sequence', 'Precursor.Charge','Precursor.Mz', 'Modified.Sequence'] 
                             ,columns = 'Run',values = ['Precursor.Quantity','RT'])

#%%Generate Inverted Series Dilution sieres
import pandas as pd
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

# Make plots
import matplotlib.pyplot as plt

import numpy as np
import mpl_scatter_density # adds projection='scatter_density'
from matplotlib.colors import LinearSegmentedColormap



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


def lin_fit(xi,x,y):
    from scipy import stats
    lin_reg = stats.linregress(x,y)
    
    m = lin_reg[0]
    b = lin_reg[1]
    r2 = lin_reg[2]**2
    
    yi = m*xi + b
    
    return yi,m,b,r2


#loop through files to get more data points

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

x_col_comp = []
y_col_comp = []
for i,file in enumerate(col1_inv_dil):
    x_temp, y_temp = prep_data(pg_pivot['PG.MaxLFQ'][ col1_inv_dil[i]],
          pg_pivot['PG.MaxLFQ'][col2_inv_dil[i] ])
    x_col_comp = x_col_comp + x_temp
    y_col_comp = y_col_comp + y_temp




fig = plt.figure(figsize = (8,3))



ax3 = fig.add_subplot(1, 2, 1, projection='scatter_density')
density = ax3.scatter_density(x_col_comp, y_col_comp, cmap=white_viridis)
fig.colorbar(density, label='Number of Points')


ax3.plot([5,43],[lin_fit(5,x_col_comp,y_col_comp)[0],lin_fit(43,x_col_comp,y_col_comp)[0]],linestyle = '--',color = 'black',linewidth = 1)

ax3.annotate('y=' + str(round(lin_fit(8,x_col_comp,y_col_comp)[1],3)) + 'x+' + str(round(lin_fit(8,x_col_comp,y_col_comp)[2],3)) +
             '\nR\u00B2=' + str(round(lin_fit(8,x_col_comp,y_col_comp)[3],3)),(1,0),fontsize = 10,ha = 'right',va = 'bottom',xycoords = 'axes fraction')




ax3.spines[['top','right']].set_visible(False)


ax3.tick_params('both',labelsize = 10)

ax3.set_title('Inverted Series - PG-level',fontsize = 10,fontweight = 'bold')
ax3.set_xlim([5,43])
ax3.set_ylim([5,43])
ax3.set_ylabel('Log2 (Column 2 Quantity)',fontsize = 10)
ax3.set_xlabel('Log2 (Column 1 Quantity)',fontsize = 10)








x_col_comp = []
y_col_comp = []
for i,file in enumerate(col1_inv_dil):
    x_temp, y_temp = prep_data(prec_pivot['Precursor.Quantity'][ col1_inv_dil[i]],
          prec_pivot['Precursor.Quantity'][col2_inv_dil[i] ])
    x_col_comp = x_col_comp + x_temp
    y_col_comp = y_col_comp + y_temp


ax4 = fig.add_subplot(1, 2, 2, projection='scatter_density')
density = ax4.scatter_density(x_col_comp, y_col_comp, cmap=white_viridis)
fig.colorbar(density, label='Number of Points')






ax4.plot([5,43],[lin_fit(5,x_col_comp,y_col_comp)[0],lin_fit(43,x_col_comp,y_col_comp)[0]],linestyle = '--',color = 'black',linewidth = 1)
ax4.annotate('y=' + str(round(lin_fit(8,x_col_comp,y_col_comp)[1],3)) + 'x+' + str(round(lin_fit(8,x_col_comp,y_col_comp)[2],3)) +
             '\nR\u00B2=' + str(round(lin_fit(8,x_col_comp,y_col_comp)[3],3)),(1,0),fontsize = 10,ha = 'right',va = 'bottom',xycoords = 'axes fraction')



ax4.spines[['top','right']].set_visible(False)


ax4.tick_params('both',labelsize = 10)

ax4.set_title('Inverted Series - Precursor-level',fontsize = 10,fontweight = 'bold')


ax4.set_xlabel('Log2 (Column 1 Peak Area)',fontsize = 10)
ax4.set_ylabel('Log2 (Column 2 Peak Area)',fontsize = 10)
ax4.set_xlim([5,43])
ax4.set_ylim([5,43])



plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
plt.tight_layout()
# fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/Revision_Round01/PythonFigures/20250715_Rd07_PG_Prec_InvSer_ColumnComp_v1.svg',dpi = 1000)
