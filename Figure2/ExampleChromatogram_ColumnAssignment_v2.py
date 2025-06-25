# -*- coding: utf-8 -*-
"""
Created on Mon Jun  9 14:47:52 2025

@author: nlancaster
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 29 11:31:55 2025

@author: nlancaster
"""
# Figure 2 - delta RT modeling
"""
Created on Fri Apr 18 08:02:28 2025

@author: nlancaster

Created on Sun Nov  3 17:11:43 2024

@author: bzhao92

Description: This is to build iRT-RT model for multi-LCs-MS proteomics data processing based on DIA-Umpire

prerequisite

It has been modified by Noah Lancaster a bit to allow returning a delta RT function that depends on RT

Noah also added an expectation value filter as recommended by the FragPipe team

"""



def iRT_RT_model(FragPipe_DIA_files,
                 col_number = 2):
    '''
    

    Parameters
    ----------
    FragPipe_DIA_files : string
        List of absolute path to .pepXML files containing FragPipe DIA files in format, {file-name}_rank{rank-number}.pepXML
    col_number : TYPE, optional
        Number of LC columns used in experiment. The default is 2.

    Returns
    -------
    iRT-RT nonlinear regression result.

    '''

    from pyteomics import pepxml
    import pandas as pd
    import numpy as np
    from scipy.interpolate import interp1d
    import statsmodels.api as sm

    # Read extended iRT peptides as pandas dataframe and then create a dict

    # Read FragPipe DIA outputs using pyteomics: https://pyteomics.readthedocs.io/en/latest/index.html
    fragpipe_combined = pd.DataFrame()
    for file in FragPipe_DIA_files:
        print(file)
        fragpipe_combined = pd.concat([fragpipe_combined,pepxml.DataFrame(file)])

    
    #added filter to the expectation value based on recommendation from Dan Polasky from FragPipe team. This 
    #also anecodotally seemed to improve the agreement with the PRTC experimental differences
    fragpipe_combined = fragpipe_combined[fragpipe_combined['expect']<=1e-5]#default starting point was 1e-5
    
    # Collapse assumed_charge and modified_sequence to form precursors
    fragpipe_combined['precursor'] = fragpipe_combined['modified_peptide'].astype(str) +\
                                 fragpipe_combined['assumed_charge'].astype(str)

    # Parse precursors to Unimod format:

    fragpipe_combined['precursor'] = fragpipe_combined['precursor'].apply(lambda x: x.replace('[160]','(Unimod:4)'))
    fragpipe_combined['precursor'] = fragpipe_combined['precursor'].apply(lambda x: x.replace('[147]','(Unimod:35)'))
    fragpipe_combined['precursor'] = fragpipe_combined['precursor'].apply(lambda x: x.replace('n[43]','(Unimod:1)'))

    # Reduce FragPipe DIA results for multiply identified precursors
    multiple_ids = fragpipe_combined['precursor'].value_counts()    
    multiple_ids = multiple_ids[multiple_ids == col_number].index

    fragpipe_combined =fragpipe_combined[fragpipe_combined['precursor'].isin(multiple_ids)]
    
    
    #filter out rt differences less than 30 seconds as these are likely to be false
    #for our typical experimental setup. Note that 30 seconds was chosen as an initial 
    #default setting and could be tuned if needed.
    rt_diff_result = (
    fragpipe_combined.sort_values(by='retention_time_sec')  # Optional: sort if order matters
    .groupby('precursor')['retention_time_sec']
    .apply(lambda x: x.iloc[1] - x.iloc[0])  # difference
    .reset_index(name='RT_diff')
    )
    rt_diff_result = rt_diff_result[rt_diff_result['RT_diff']>= 30]

    

    #NOTE THAT WITH MY MODIFICATIONS FOR FILTERING THE DUPLICATE FEATURES, THE CODE
    #WILL ONLY WORK FOR DUAL COLUMN ANALYSES - Noah Lancaster, 23Apr2025
    
    # Find col_numberly identified precursors shared in extended calibration set
    # anchors = list(set(rt_diff_result['precursor']).intersection(set(extended_iRT_dict)))
    anchors = list(set(rt_diff_result['precursor']))

    # fragpipe_combined_reduced = fragpipe_combined[fragpipe_combined['precursor'].isin(anchors)]
    fragpipe_combined_reduced = fragpipe_combined
    fragpipe_combined_reduced_sorted = fragpipe_combined_reduced.sort_values(by='retention_time_sec')
       
    
    RTs = [[] for _ in range(col_number)] 

    for x in anchors:
        rts = list((fragpipe_combined_reduced_sorted[fragpipe_combined_reduced_sorted['precursor']==x]['retention_time_sec'])/60)
        for i in range(col_number):
            RTs[i].append(rts[i])
        
    #AFTER HERE NML MODIFIED FOR OTHER PURPOSES
    
    functions = []
    x = []
    y = []
    for i in range(col_number):
        lowess = sm.nonparametric.lowess(RTs[i], RTs[0], frac=.035) 


        lowess_x = list(zip(*lowess))[0]
        lowess_y = list(zip(*lowess))[1]

        
        f1 = interp1d(lowess_x, lowess_y, bounds_error=False)  
        xnew = np.linspace(min(RTs[0]), max(RTs[0]), 100000) #Generating 100000 points for predicting ynew using the model
        ynew = f1(xnew)
        x.append(xnew)
        y.append(ynew)
        functions.append(f1)

    def ydiff_fxn(rt):
        return functions[1](rt) - functions[0](rt)
    
    #only include this return statement for plots of model accuracy
    ydiff = np.array(y[1])-np.array(y[0])

    
    return (x[0],ydiff),ydiff_fxn




# create plots for model accuracy assessment - Rd07
fragpipe_dia_files = ["P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_FragPipeDIA_Searches/20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1_rank1.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_FragPipeDIA_Searches/20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1_rank2.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_FragPipeDIA_Searches/20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1_rank3.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_FragPipeDIA_Searches/20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1_rank4.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_FragPipeDIA_Searches/20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1_rank5.pepXML"]

plot_points, functions = iRT_RT_model(fragpipe_dia_files)
#%% Figure 2 - DIANN report parsing
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


#perform a pivot table to go to the wide format
prec_pivot = pd.pivot(report,index=['Precursor.Id', 'Stripped.Sequence', 'Precursor.Charge','Precursor.Mz', 'Modified.Sequence'] 
                             ,columns = 'Run',values = ['Precursor.Quantity','RT'])


#%% Figure 2 - raw data file parsing

#Use 200ng mouse dataset to show example chromatograms
cm_only = "P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/raw/20241104_Round07_DualColumn/LCMS_CookieMonster_ver02/20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC2.raw"
zoe_only = "P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/raw/20241104_Round07_DualColumn/LCMS_CookieMonster_ver02/20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC2.raw"
combined = 'P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/raw/20241104_Round07_DualColumn/LCMS_CookieMonster_ver02/20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_1.raw'


def ms1_tic(rawfile_path):
    from pymsfilereader import MSFileReader
    raw_file = MSFileReader(rawfile_path)
    
    total_scans = raw_file.NumSpectra
    
    rt = []
    tic  = []
    for x in range(1,total_scans):
        if raw_file.GetMSOrderForScanNum(x)==1 and raw_file.RTFromScanNum(x)<45:
            rt.append(raw_file.RTFromScanNum(x))
            tic.append(raw_file.GetScanHeaderInfoForScanNum(x)['TIC'])
            
    return rt,tic

cm_data = ms1_tic(cm_only)
zoe_data = ms1_tic(zoe_only)
comb_data = ms1_tic(combined)

#get integrated TIC signal
import numpy as np
cm_int = np.format_float_scientific(sum(cm_data[1]),2)
zoe_int = np.format_float_scientific(sum(zoe_data[1]),2)
comb_int = np.format_float_scientific(sum(comb_data[1]),2)

#%% Figure 2 - Plot generation
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable


layout = [
    ['A','A','B'],
    ['A','A','C'],
    ['E','D','D'],
    ['F','D','D']]
fig, axd = plt.subplot_mosaic(layout, figsize=(7.5,7),layout = 'constrained')

ax1 = axd['A']
ax2 = axd['B']
ax3 = axd['C']
ax4 = axd['D']
ax5 = axd['E']
ax6 = axd['F']


axes = (ax1,ax2,ax3,ax4,ax5, ax6)
for ax in axes:
    ax.spines[['top','right']].set_visible(False)
    


#make plots of chromatograms
data= [[cm_data[0],[x+3.1e10 for x in cm_data[1]]],[zoe_data[0],[x + 1.6e10 for x in zoe_data[1]]],comb_data]
data= [cm_data,zoe_data,comb_data]
divider = make_axes_locatable(ax1)
axTop = divider.append_axes('top',size = '100%',pad = 0.2)
axBottom = divider.append_axes('bottom',size = '100%',pad = 0.2)

axTop.set_xticks([])
ax1.set_xticks([])


axTop.plot(data[0][0],data[0][1],color = 'black',linewidth = 1)
ax1.plot(data[1][0],data[1][1],color = 'black',linewidth = 1)
axBottom.plot(data[2][0],data[2][1],color = 'black',linewidth = 1)

axes = (axTop,ax1,axBottom)
for ax in axes:
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlim(0,45)
    ax.set_ylim(0,1.7e10)
    ax.set_ylabel('TIC',fontsize = 10)
    ax.set_yticks([0.0,0.5e10,1e10,1.5e10])

axBottom.set_xlabel('Retention Time (min)')
axBottom.set_xticks([0,10,20,30,40])


labels = ['Column 1: 200ng\nColumn 2: Blank','Column 1: Blank\nColumn 2: 200ng','Column 1: 200ng\nColumn 2: 200ng']


axTop.annotate(labels[0],(0.73,1),va = 'top',ha = 'left',xycoords = 'axes fraction',fontsize = 10)
ax1.annotate(labels[1],(0.73,1),va = 'top',ha = 'left',xycoords = 'axes fraction',fontsize = 10)
axBottom.annotate(labels[2],(0.73,1),va = 'top',ha = 'left',xycoords = 'axes fraction',fontsize = 10)


axTop.annotate(f'Integrated\nTIC: {cm_int}',(0.01,0.98),va = 'top',ha = 'left',xycoords = 'axes fraction',fontsize = 10)
ax1.annotate(f'Integrated\nTIC: {zoe_int}',(0.01,0.98),va = 'top',ha = 'left',xycoords = 'axes fraction',fontsize = 10)
axBottom.annotate(f'Integrated\nTIC: {comb_int}',(0.01,0.98),va = 'top',ha = 'left',xycoords = 'axes fraction',fontsize = 10)



#make delta RT model plots
prtc_rt = pd.read_excel('P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_DualCol_PRTCcalc/20250110_DualColumn_PRTC_RT_results.xlsx')

ax2.spines[['top','right']].set_visible(False)
ax2.plot(plot_points[0],plot_points[1],color = 'black',label = 'Model',linewidth = 1)
ax2.scatter(prtc_rt['CM Average'],prtc_rt['delta RT'],marker = 'o',label = 'Experiment',s = 5,linewidth = 4)
ax2.legend(loc = 'lower left',frameon = False,borderpad = 0,fontsize = 10,handlelength = 1,labelspacing = 0.2)


ax2.set_xlim(0,42)
ax2.set_ylim(0,)
ax2.set_xlabel('Retention Time (min)',fontsize =10 )
ax2.set_ylabel('\u0394RT (min)',fontsize =10)


predicted_rt_diff = [functions(x) for x in prtc_rt['CM Average']]

rt_diff_error = []
for i,x in enumerate(prtc_rt['delta RT']):
    rt_diff_error.append(x-predicted_rt_diff[i])

ax3.spines[['top','right']].set_visible(False)   
ax3.scatter(prtc_rt['CM Average'],rt_diff_error,marker = 'o',s = 5,linewidth = 4)
ax3.set_xlabel('Retention Time (min)',fontsize =10)
ax3.set_ylabel('\u0394RT Error (min)',fontsize =10)
ax3.set_xlim(0,42)
ax3.set_ylim(-0.3,0.5)

ax2.tick_params(axis = 'both',which = 'major',labelsize = 10)
ax3.tick_params(axis = 'both',which = 'major',labelsize = 10)
ax2.set_xticks([0,5,10,15,20,25,30,35,40])
ax3.set_xticks([0,5,10,15,20,25,30,35,40])
ax2.set_yticks([0.0,1.0,2.0,3.0,4.0],[0,1,2,3,4])
ax3.set_yticks([-0.2,0.0,0.2,0.4])


#make plots of column assignment accuracy
column1_results = report[report['Run'].str.contains('[4.0]')]
column2_results = report[report['Run'].str.contains('[4.0]')]
qc_col1_rt1 = len(column1_results[column1_results['Run'] == '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3']['Precursor.Id'].unique())
qc_col1_rt2 = len(column2_results[column2_results['Run'] == '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3']['Precursor.Id'].unique())
qc_col2_rt1 = len(column1_results[column1_results['Run'] == '[0.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3']['Precursor.Id'].unique())
qc_col2_rt2 = len(column2_results[column2_results['Run'] == '[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3']['Precursor.Id'].unique())
comb_col1_rt1 = len(column1_results[column1_results['Run'] == '[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']['Precursor.Id'].unique())
comb_col2_rt2 = len(column2_results[column1_results['Run'] == '[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3']['Precursor.Id'].unique())


labels = ['Column 1\nInjection\nOnly','Column 2\nInjection\nOnly','Combined\nInjection']
x = np.arange(3)

width = 0.45

ax4.bar(x-width/2,[qc_col1_rt1,qc_col2_rt1,comb_col1_rt1],width = width, label = 'Column 1',color = '#6A6A6A',align = 'center' )
ax4.bar(x+width/2,[qc_col1_rt2,qc_col2_rt2,comb_col2_rt2],width = width, label = 'Column 2', color = '#AAAAAA',align = 'center' )
ax4.set_xticks(x,labels,rotation = 0)
ax4.set_xticks(x,labels,rotation = 0)
ax4.set_ylabel('Precursors',fontsize = 10)
ax4.tick_params('both',labelsize = 10)

ylim_max = 200000
ax4.set_ylim(0,ylim_max)

for i,val in enumerate([qc_col1_rt1,qc_col2_rt1,comb_col1_rt1]):
    ax4.annotate(str(val),(x[i]-width/2,val+40),ha = 'center',fontsize = 10)


#add labels and lines to the bar chart 
ymin = qc_col1_rt2
ymax = qc_col1_rt1
ax4.axvline(0.25,ymin = ymin/ylim_max, ymax= ymax/ylim_max,linestyle = '--',color = 'black',linewidth = 1)
ax4.plot((0.05,0.45),(ymax,ymax),color = 'black',linewidth = 1)
pcnt = str(round(ymin/(ymin + ymax)*100,1))
ax4.annotate(pcnt + "%",(0.25,(ymax-ymin)/2+ymin),va = 'center',ha = 'center',fontsize = 10)

ax4.ticklabel_format(axis = 'y',style = 'scientific',scilimits = (0,1))

for i,val in enumerate([qc_col1_rt2,qc_col2_rt2,comb_col2_rt2]):
    t = ax4.annotate(str(val),(x[i]+width/2,val+40),ha = 'center',fontsize = 10)


#add labels and lines to the bar chart 
ymin = qc_col2_rt1
ymax = qc_col2_rt2
ax4.axvline(0.75,ymin = ymin/ylim_max, ymax= ymax/ylim_max,linestyle = '--',color = 'black',linewidth = 1)
ax4.plot((0.55,0.95),(ymax,ymax),color = 'black',linewidth = 1)
pcnt = str(round(ymin/(ymin + ymax)*100,1))
ax4.annotate(pcnt + "%",(0.75,(ymax-ymin)/2 + ymin),va = 'center',ha = 'center',fontsize = 10)

ax4.legend(title = 'Predicted RTs',frameon = False,fontsize = 10,title_fontsize = 10,
           handletextpad = 0.2,labelspacing=0.2,handlelength = 1, loc ='upper right',
           ncol = 2, columnspacing = 0.5)



#make RT scatter plots

col1_data = [prec_pivot['RT']['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_Blank_CM_QC3'].tolist(),
            prec_pivot['RT']['[0.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3'].tolist()]

col2_data = [prec_pivot['RT']['[4.0]20241216_Rd07_CM_Blank_Zoe_200ngMouse_Zoe_QC3'].tolist(),
            prec_pivot['RT']['[4.0]20241216_Rd07_CM_200ngMouse_Zoe_200ngMouse_3'].tolist()]



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

ax5.plot([5,43],[lin_fit(5,col1_data[0],col1_data[1])[0],lin_fit(43,col1_data[0],col1_data[1])[0]],linestyle = '--',color = 'black')
ax5.annotate('y=' + str(lin_fit(5,col1_data[0],col1_data[1])[1]) + 'x+' + str(lin_fit(5,col1_data[0],col1_data[1])[2]) +
             '\nR\u00B2=' + str(lin_fit(5,col1_data[0],col1_data[1])[3]),(1,0),fontsize = 10,ha = 'right',va = 'bottom',xycoords = 'axes fraction')

ax6.plot([5,43],[lin_fit(5,col2_data[0],col2_data[1])[0],lin_fit(43,col2_data[0],col2_data[1])[0]],linestyle = '--',color = 'black')
ax6.annotate('y=' + str(lin_fit(8,col2_data[0],col2_data[1])[1]) + 'x+' + str(lin_fit(8,col2_data[0],col2_data[1])[2]) +
             '\nR\u00B2=' + str(lin_fit(8,col2_data[0],col2_data[1])[3]),(1,0),fontsize = 10,ha = 'right',va = 'bottom',xycoords = 'axes fraction')



ax5.spines[['top','right']].set_visible(False)
ax5.scatter(col1_data[0],col1_data[1],marker = '.',color = '#2CA8E0', linewidth = 0.5,s = 15,rasterized=True)

ax6.spines[['top','right']].set_visible(False)
ax6.scatter(col2_data[0],col2_data[1],marker = '.',color = '#2CA8E0', linewidth = 0.5, s = 15,rasterized=True)

ax5.tick_params('both',labelsize = 10)
ax6.tick_params('both',labelsize = 10)
ax5.set_title('Column 1',fontsize = 12)
ax6.set_title('Column 2',fontsize = 12)


ax5.set_xlabel('Single Injection Ref RT (min)',fontsize = 10)
ax5.set_ylabel('Dual Injection\nRT (min)',fontsize = 10)

ax6.set_xlabel('Single Injection Ref RT (min)',fontsize = 10)
ax6.set_ylabel('Dual Injection\nRT (min)',fontsize = 10)

ax5.set_xlim([5,43])
ax5.set_ylim([5,43])
ax6.set_xlim([5,43])
ax6.set_ylim([5,43])





plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/20250610_Rd07_PrecursorRT_Assignment_DIANN_v5.svg',dpi = 1000)

