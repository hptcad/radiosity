import pandas as pd
import numpy as np
import os

def Extract(diameter=500,depth=2500,sticking=0.02,basepath=''):
    if basepath == '':
        basepath = os.path.abspath(os.path.dirname(__file__))
        basepath = os.path.join(basepath,'viennats_data/')
        
    eps = 1e-12
    df = pd.read_csv(basepath+'holescell_'+str(diameter)+'_'+str(depth)+'_sp'+str(int(sticking*100))+'_rates_griddelta.csv',header=None,delimiter=' ')
    df.columns = ['X', 'Y','Z','ABSORBED']
    top = df.loc[np.abs(df['Z']-df.Z.min())<eps] # select only top layer
    bottom = df.loc[np.abs(df.Z-df.Z.max())<eps] # select only bottom layer
    wall = df.loc[(np.abs((df.Z-df.Z.max()))>eps) & ((np.abs(df.Z-df.Z.min()))>eps)] # select only wall
    lines = wall.groupby(['X', 'Y'])  # group by xy coordinates
    
    # extract average/min/max line values
    groups = dict(list(lines))
    groupnames = list(groups.keys())
    ngroups = len(groupnames)
    firstline = groupnames[0]
    n = groups[firstline].shape[0]
    #print(n)
    #print(ngroups)
    avglinez = (groups[firstline]['Z'].values + depth*1e-7) / (depth*1e-7)
    alllines = np.zeros((n,ngroups))
    i = 0
    for name, line in lines:
        data = np.array(line['ABSORBED'].values)
        alllines[:,i] = data
        i = i + 1 
        
    avglinevalues = np.average(alllines,axis=1)      
    maxlinevalues = np.max(alllines,axis=1) 
    minlinevalues = np.min(alllines,axis=1) 
    linez = []
    linevalues = []
    for name, line in lines: # add start/endpoints from top/bottom layers additionally to each line
        #print(name)
        topex = top.loc[(np.abs(top['X']-name[0])<eps ) & (np.abs(top['Y']-name[1])<eps )]
        bottomex = bottom.loc[(np.abs(bottom['X']-name[0])<eps ) & (np.abs(bottom['Y']-name[1])<eps )]
        frames = [topex, line,  bottomex]
        line = pd.concat(frames)
        #line.sort_values(by=['Z'], ascending=[1])
        linez.append((line['Z'].values + depth*1e-7) / (depth*1e-7) )
        linevalues.append(line['ABSORBED'].values)
        #print(line.shape)  
     
    # sort bottom points using radius
    bottomplot = bottom.copy()    
    bottomplot['radius'] = np.sqrt(bottomplot['X']**2+bottomplot['Y']**2)
    bottomplot = bottomplot.sort_values(by=['radius'], ascending=[0])
    # sort top points using radius
    topplot = top.copy()  
    topplot['radius'] = np.sqrt(topplot['X']**2+topplot['Y']**2)
    topplot = topplot.sort_values(by=['radius'], ascending=[0])
    
    radius = (diameter*1e-7*0.5)
    topr = topplot['radius'].values / radius 
    topvalues = topplot['ABSORBED'].values
    bottomr = bottomplot['radius'].values / radius
    bottomvalues = bottomplot['ABSORBED'].values 
    
    return topr, topvalues, bottomr, bottomvalues, linez, linevalues, avglinez, avglinevalues, minlinevalues, maxlinevalues