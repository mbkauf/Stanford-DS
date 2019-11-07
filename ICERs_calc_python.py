# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 15:25:11 2019

@author: Tess Ryckman
"""
#NOTES
#this file constains code to calculate ICERs, including removing dominated/extended dominated strategies
#results is a pandas dataframe w/ columns w/ discounted costs & discounted QALYs for each strategy
#n is the number of strategies (should be same as the # of rows in results)
    
import numpy as np
import pandas as pd

os.chdir('[directory name]')
results=pd.read_csv("CEA_Results.csv", delimiter=",", header=0)
n=np.shape(results)[0]

#Remove Dominated and Extended Dominated Strategies#
#also see: http://journals.sagepub.com/doi/abs/10.1177/0272989X15583496?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%3dpubmed#articleShareContainer
results=results.sort_values('DiscCosts')
results=results.reset_index(drop=True)
results['IncrCosts']=results.DiscCosts.diff()
results['IncrQALYs']=results.DiscQALYs.diff()
results_small=results[(results.IncrQALYs>0)|(results['IncrQALYs'].isnull())].copy() #removing dominated options (negative incr. QALYs)
results_small['ICERs']=results_small.IncrCosts/results_small.IncrQALYs
results_small['IncrICER']=results_small.ICERs.diff()
for j in range(0,n):
    for k in range(0,n):
        results_small['IncrCosts']=results_small.DiscCosts.diff()
        results_small['IncrQALYs']=results_small.DiscQALYs.diff()
        results_small=results_small[(results_small.IncrQALYs>0)|(results_small['IncrQALYs'].isnull())].copy() #remove dominated options (negative incr. QALYs) & recalculate
    results_small['ICERs']=results_small['IncrCosts']/results_small['IncrQALYs']
    results_small['IncrICER']=results_small.ICERs.diff() #negative "incremental ICER" indicates extended dominance
    results_small['Drop']=(results_small.IncrICER.shift(-1)<0)
    results_small['DropCount']=(results_small.groupby('Drop').cumcount()+1)*(results_small['Drop']==1)
    results_small=results_small[(results_small.DropCount!=1)].copy() #delete strategy if next strategy has a lower ICER than it (extended dominated)
   
results_small.drop(columns=['IncrICER', 'Drop','DropCount'], inplace=True) #remove unnecessary columns
CEA_results=results_small.reset_index(drop=True) #file with output for non-dominated strategies only
