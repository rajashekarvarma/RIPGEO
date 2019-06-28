import GEOparse # Python package to upload a geo data
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats import multitest
import seaborn as sns
import matplotlib.pyplot as plt


############### Fetch Agilent data ###############
print('\n\n',"******...Hi Welcome to RIPGEO...******",'\n\n')
GSE_ID = input('Please enter your GSE ID (ex:GSE62893): ')

print('\n',"Provided GSE ID: ",GSE_ID)

print('\n',"Intitating data extraction...",'\n\n')

gse = GEOparse.get_GEO(geo=GSE_ID, destdir="./")
plt_name=[]
# print(gse.gpls)
for pl_name, pl in gse.gpls.items():
    plt_name.append(pl_name)
plt_name=''.join(plt_name)

print("Platform Name:", plt_name)

pivoted_control_samples = gse.pivot_samples('VALUE')
# print(pivoted_control_samples.head())

######## Filter probes that are not expressed worst 25% genes are filtered out

pivoted_control_samples_average = pivoted_control_samples.median(axis=1)
# print("Number of probes before filtering: ", len(pivoted_control_samples_average))
expression_threshold = pivoted_control_samples_average.quantile(0.25)
expressed_probes = pivoted_control_samples_average[pivoted_control_samples_average >= expression_threshold].index.tolist()
# print("Number of probes above threshold: ", len(expressed_probes))

samples = gse.pivot_samples("VALUE").ix[expressed_probes]
# print(samples.head())
# print(gse.gpls[plt_name].table)

######## Annotate matrix table

samples_annotated = samples.reset_index().merge(gse.gpls[plt_name].table[["ID", "GB_ACC"]], left_on='ID_REF', right_on="ID").set_index('ID_REF')

# print(samples_annotated.head())
del samples_annotated["ID"]
# print(samples_annotated.head())
samples_annotated = samples_annotated.dropna(subset=["GB_ACC"])
samples_annotated = samples_annotated[~samples_annotated.GB_ACC.str.contains("///")]
samples_annotated = samples_annotated.groupby("GB_ACC").median()
# print(samples_annotated.index)

print('\n','Column names from the matrix: ',samples_annotated.columns)

######## Extract matrix data to a csv file
exprs = []
gsmNames = []
metadata = {}

for gsm_name, gsm in gse.gsms.items():
    # print(gsm.metadata['type'][0])
    if gsm.metadata['type'][0]=='RNA':
        # Expression data
        if len(gsm.table)>0:
            tmp = gsm.table['VALUE']
            # print(tmp)
            tmp.index = gsm.table['ID_REF']
            gsmNames.append(gsm_name)
            if len(exprs)==0:
                exprs = tmp.to_frame()
            else:
                exprs = pd.concat([exprs,tmp.to_frame()],axis=1)

print('\n','Extracting metadata...','\n')

######## extract metadata to csv file

for gsm_name, gsm in gse.gsms.items():
    if gsm.metadata['type'][0]=='RNA':
                for key,value in gsm.metadata.items():
                # print(key)
                # print(value)
                    if (key=='characteristics_ch1' or key=='characteristics_ch2') and (len([i for i in value if i!=''])>1 or value[0].find(': ')!=-1):
                        # print(value)
                        tmpVal = 0
                        for tmp in value:
                            splitUp = [i.strip() for i in tmp.split(':')]
                            # print(splitUp)
                            if len(splitUp)==2:
                                if not splitUp[0] in metadata:
                                    metadata[splitUp[0]] = {}
                                metadata[splitUp[0]][gsm_name] = splitUp[1]
                            else:
                                if not key in metadata:
                                    metadata[key] = {}
                                metadata[key][gsm_name] = splitUp[0]
                    else:
                        if not key in metadata:
                            metadata[key] = {}
                        if len(value)==1:
                            metadata[key][gsm_name] = ' '.join([j.replace(',',' ') for j in value])

# Write expression data matrix to file
exprs.columns = gsmNames
with open(GSE_ID+'exprs.csv','w') as outFile:
    exprs.to_csv(outFile)

# Write metadata matrix to file
with open(GSE_ID+'metadata.csv','w') as outFile:
    outFile.write('Metadata,'+','.join(gsmNames))
    for key in metadata:
        tmp = [key]
        for gsm_name in gsmNames:
            if gsm_name in metadata[key]:
                tmp.append(metadata[key][gsm_name])
            else:
                tmp.append('NA')
        outFile.write('\n'+','.join(tmp))
print('\n','Data matrix and metadata for',GSE_ID,'have been written to',GSE_ID+'exprs.csv',GSE_ID+'metadata.csv @ cwd','\n')

######## select control and test sample columns

samples_annotated = samples_annotated.astype(float)

control_sample = input('Please enter column numbers of control samples (ex:0,2,4): ')
control_samples = control_sample.split(',')
control_samples = [int(i) for i in control_samples]
Test_sample =  input('Please enter column numbers of test samples (ex:3,5,7): ')
Test_samples = Test_sample.split(',')
Test_samples = [int(i) for i in Test_samples]

print('\n','control samples column names entered:',control_samples,'\n')
print('\n','Test samples column names entered:',Test_samples,'\n')

######## perform t-test for the data

print('\n','Performing independednt T Test on the selected data...'+'\n')
samples_annotated["Ttest"] = stats.ttest_ind(samples_annotated.T.iloc[Test_samples, :],samples_annotated.T.iloc[control_samples, :], equal_var=True, nan_policy="omit")[1]

######## perform anova
# samples_annotated['Anova_one'] = [stats.f_oneway(samples_annotated.T.iloc[Test_samples, x],samples_annotated.T.iloc[control_samples, x])[1] for x in range(samples_annotated.shape[0])]
# samples_annotated['Ttest'].to_csv('pvalues.csv')

######## filter data based FDR (<0.05)

samples_annotated["FDR"] = multitest.multipletests(samples_annotated['Ttest'], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
# print(samples_annotated.head())
samples_annotated = samples_annotated.sort_values(by="FDR")
filtered_samples = samples_annotated[samples_annotated["FDR"] < 0.05]
print('\n','Number of genes remaining after FDR filter of 0.05:',len(filtered_samples)) 
f_samples = pd.DataFrame() 
f_samples['control'] = filtered_samples.T.iloc[Test_samples,:].mean() 
f_samples['test'] = filtered_samples.T.iloc[control_samples,:].mean() 
f_samples['p-value'] = filtered_samples['Ttest']
f_samples['FDR'] = filtered_samples['FDR']
######## calculate log2FC

f_samples['log2FC'] = f_samples['test'].apply(np.log2) - f_samples['control'].apply(np.log2) 

print('\n','Calculating Log2 values of the data...','\n')

######## filter gene list based on log2FC

up_c = float(input('Please enter log2FC cutoff value for up regulation(ex:0.5): '))
dwn_c = float(input('Please enter log2FC cutoff value for down regulation(ex:-0.5): '))

diff_up = f_samples[f_samples['log2FC'] >= up_c]
diff_down = f_samples[f_samples['log2FC'] <= dwn_c]

######## write up and down regulated genes to csv

diff_up.to_csv('Upregulated_genes.csv')
diff_down.to_csv('Downregulated_genes.csv')

######## plot log difference of upregulated and down regulated genes

plot_y = input('Do you want to plot bar plot for log2 fold difference (yes/no): ')

if plot_y == 'yes':
    diff = pd.concat([diff_up,diff_down])
    diff_vals = diff['log2FC'].sort_values()
    counter = np.arange(len(diff_vals.values))

    fig, ax = plt.subplots(figsize = (20,10))
    ax.bar(counter,diff_vals.values, width=0.5)
    ax.set_xticks(counter)
    ax.set_xticklabels(diff_vals.index.values, rotation=90)
    ax.set_title("Gene expression differences of Control vs Test")
    ax.set_ylabel("log2 difference")
    plt.show()
    print('\n','Task completed...Output written successfully to current working directory.')

else:
    print('\n','Task completed...Output written successfully to current working directory.')

