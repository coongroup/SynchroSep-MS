# -*- coding: utf-8 -*-
"""
Created on Tue Jun 10 09:57:29 2025

@author: nlancaster
"""

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
    fragpipe_combined = fragpipe_combined[fragpipe_combined['expect']<=1e-5]#default starting point is 1e-5
    
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

    return (x[0],ydiff),ydiff_fxn,(RTs[0],np.array(RTs[1])-np.array(RTs[0]))

# create plots for model accuracy assessment - Rd07
files_0min = ["P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank5.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank1.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank2.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank3.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank4.pepXML"]

files_1min = ["P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank5.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank1.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank2.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank3.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank4.pepXML"]

files_2min = ["P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank5.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank1.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank2.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank3.pepXML",
"P:/Projects/NML_2024_Internal_MultiColumn/Proteomics/data/software_output/Round07_GradOpt_FragPipeDIA/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank4.pepXML"]
#


files_0min = ["D:/Rd07_GradOptFragPipeDIA_Search/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank5.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank1.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank2.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank3.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/NoDelay/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_NoDelay_20241212104519_rank4.pepXML"]

files_1min = ["D:/Rd07_GradOptFragPipeDIA_Search/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank5.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank1.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank2.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank3.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/1min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_1min_rank4.pepXML"]

files_2min = ["D:/Rd07_GradOptFragPipeDIA_Search/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank5.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank1.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank2.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank3.pepXML",
"D:/Rd07_GradOptFragPipeDIA_Search/2min/20241212_DualColRd07_GradTest_CM_250ngHC_Zoe_250ngHC_2min_rank4.pepXML"]
#


plot_points_0min, functions, scatter_0min = iRT_RT_model(files_0min)
plot_points_1min, functions, scatter_1min = iRT_RT_model(files_1min)
plot_points_2min, functions, scatter_2min = iRT_RT_model(files_2min)
#%% save plot points
file_root = 'C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles//'

def save_file(points,name,file_root):
    import pandas as pd
    df = pd.DataFrame({'1':points[0],'2':points[1]})
    df.to_csv(file_root + name)
    
    
    

save_file(plot_points_0min,'deltaRT_Tune_model_points_0min_v2.csv',file_root)
save_file(plot_points_1min,'deltaRT_Tune_model_points_1min_v2.csv',file_root)
save_file(plot_points_2min,'deltaRT_Tune_model_points_2min_v2.csv',file_root)
save_file(scatter_0min,'deltaRT_Tune_scatter_points_0min_v2.csv',file_root)
save_file(scatter_1min,'deltaRT_Tune_scatter_points_1min_v2.csv',file_root)
save_file(scatter_2min,'deltaRT_Tune_scatter_points_2min_v2.csv',file_root)
#%% read in files
file_root = 'C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/IntermediateFiles//'

def read_file(name,file_root):
    import pandas as pd
    df = pd.read_csv(file_root + name)
    return (list(df['1']),list(df['2']))

plot_points_0min = read_file('deltaRT_Tune_model_points_0min_v2.csv',file_root)
plot_points_1min = read_file('deltaRT_Tune_model_points_1min_v2.csv',file_root)
plot_points_2min = read_file('deltaRT_Tune_model_points_2min_v2.csv',file_root)
scatter_0min = read_file('deltaRT_Tune_scatter_points_0min_v2.csv',file_root)
scatter_1min = read_file('deltaRT_Tune_scatter_points_1min_v2.csv',file_root)
scatter_2min = read_file('deltaRT_Tune_scatter_points_2min_v2.csv',file_root)





#%%
import matplotlib.pyplot as plt
fig,ax1 = plt.subplots(figsize = (3,3))

ax1.spines[['top','right']].set_visible(False)

#blank plot for legend entry
ax1.plot([],[],color = 'black',label = 'Model',linewidth = 1.2)

ax1.scatter(scatter_0min[0],scatter_0min[1],color = 'royalblue',label = '0 min Added Delay',marker = '.',s = 10,alpha = 0.4,
            rasterized = True, edgecolor = 'none')
ax1.plot(plot_points_0min[0],plot_points_0min[1],color = 'black',label = '_nolegend_',linewidth = 1.2)
ax1.scatter(scatter_1min[0],scatter_1min[1],color = 'maroon',label = '1 min Added Delay',marker = '.',s = 10,alpha = 0.4,
            rasterized = True,edgecolor = 'none')
ax1.plot(plot_points_1min[0],plot_points_1min[1],color = 'black',label = '_nolegend_',linewidth = 1.2)
ax1.scatter(scatter_2min[0],scatter_2min[1],color = 'grey',label = '2 min Added Delay',marker = '.',s = 10,alpha = 0.4,
            rasterized = True,edgecolor = 'none')
ax1.plot(plot_points_2min[0],plot_points_2min[1],color = 'black',label = '_nolegend_',linewidth = 1.2)

ax1.set_xlabel('Retention Time (min)',fontsize = 10)
ax1.set_ylabel('\u0394RT (min)')

ax1.set_ylim(0,8.5)
ax1.set_xlim(0,40)

ax1.legend(frameon = False, fontsize = 10,labelspacing = 0.2,handletextpad = 0.3,
           markerscale = 2,handlelength = 1)


plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
fig.tight_layout()
fig.savefig('C:/Users/nlancaster/OneDrive - UW-Madison/MultiColumnManuscript/ManuscriptDraft/PythonFigures/20250610_Rd07_GradOpt_deltaRT_Tuning_v4.svg',dpi = 1200)

