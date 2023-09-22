#!/usr/bin/env python
# coding: utf-8

# In[1]:


exec(open('imports.py').read())


# * Load data 
# * find the difference between the means of the various variables 
# * plot distributions of the data 

# In[ ]:





# In[2]:


# load the bootstrapped data into collections 
cascade_bard_v1 = collections.defaultdict(list)
guan_waliser = collections.defaultdict(list)
mundhenk = collections.defaultdict(list)
reid250=collections.defaultdict(list)

flds=['cascade_bard_v1','guan_waliser','mundhenk','reid250'] #algorithm names 
vbls=['IVT','MSL','TCWV','PV','Z','T'] #variable names 

#loop through algorithms 
for f in flds:
    for v in vbls:
        eval(f)[v].append(xr.open_dataset(f'Bootstrapped Data/{f}_boot-{v}.nc')) #compute the means of the data 


# In[3]:


import matplotlib.ticker as mtk


# In[41]:


import cmocean.cm as cmo


# In[51]:


#POP : Make a plot that has most of the atmospheric paramters under consideration on it
fig,axes = plt.subplots(figsize=(16,10), ncols =2,nrows=2, subplot_kw={'projection':ccrs.PlateCarree()})
title=['TECA_BARD v1','Guan & Waliser','Mundhenk v3','Reid 250']

for ax,f,tt in zip(axes.flat,flds,title):
    
    eval(f)['IVT'][0]['IVT'].mean('tstep').plot(ax=ax, cmap=cmo.delta,cbar_kwargs={'orientation':'horizontal','label':'IVT ($kgm^{-2}s^{-2}$)'},vmin=-300,vmax=300)
    
    cb=(eval(f)['MSL'][0]['MSL']/100).mean('tstep').plot.contour(ax=ax,cmap='grey',linestyles='solid', linewidth=0.5)
    fmt = mtk.StrMethodFormatter("{x:1.0f} $hPa$")
    cb.clabel(fmt=fmt,inline = True,fontsize=8,rightside_up=True)
    
    cbb=(eval(f)['Z'][0]['Z']/1000).mean(['tstep','level']).plot.contour(ax=ax,cmap='spring',linestyles='dotted', linewidth=0.5)
    fmt = mtk.StrMethodFormatter("{x:2.1f} $km$")
    cbb.clabel(fmt=fmt,inline = True,fontsize=8,rightside_up=True)
    
    cbb=eval(f)['TCWV'][0]['TCWV'].mean(['tstep']).plot.contour(ax=ax,cmap='#DFFF00',linestyles='dotted',linewidth=0.1,vmin=4,levels=np.arange(4,11,3))
    fmt = mtk.StrMethodFormatter("{x:1.0f} $kg/m^2$")
    cbb.clabel(fmt=fmt,inline = True,fontsize=8,colors='#DFFF00',rightside_up=True)
    
    #cbb=(eval(f)['PV'][0]['PV']/10**-6).mean(['tstep','level']).plot.contour(ax=ax,cmap='#7C0A02',linestyles='solid',levels=np.arange(0.01,2.2,0.2),linewidth=0.1) #divide by 10^-6 to get PV units 
    #fmt = mtk.StrMethodFormatter("{x:2.2f} PVU")
    #cbb.clabel(fmt=fmt,inline = True,fontsize=8,colors='#7C0A02',rightside_up=True)
    
    cbb=eval(f)['T'][0]['T'].mean(['tstep','level']).plot.contour(ax=ax,cmap='bwr',linestyles='solid',linewidth=0.1) #divide by 10^-6 to get PV units 
    fmt = mtk.StrMethodFormatter("{x:3.0f} K")
    cbb.clabel(fmt=fmt,inline = True,fontsize=8,colors='red',rightside_up=True)
    
    ax.set_title(f'{tt}')

    cart_plot(ax,rlon=-155)   
plt.savefig('Images/Compund_plot-[IVT,MSL,Z,TCWV,T]2.png',facecolor='white',dpi=350)

'''
# # WORK ON BOOTSTRAPPED DATA

# In[2]:


# load the bootstrapped data into collections 
cascade_bard_v1 = collections.defaultdict(list)
guan_waliser = collections.defaultdict(list)
mundhenk = collections.defaultdict(list)
reid250=collections.defaultdict(list)

flds=['cascade_bard_v1','guan_waliser','mundhenk','reid250'] #algorithm names 
vbls=['IVT','MSL','TCWV','PV','Z','T'] #variable names 

#loop through algorithms 
for f in flds:
    for v in vbls:
        eval(f)[v].append(xr.open_dataset(f'Bootstrapped Data/{f}_boot-{v}.nc')) #compute the means of the data 


# In[3]:


#create a collection of the differences between the means for all algorithms 
diffs_IVT = collections.defaultdict(list) ; diffs_PV = collections.defaultdict(list) 
diffs_MSL = collections.defaultdict(list) ; diffs_Z = collections.defaultdict(list)
diffs_TCWV = collections.defaultdict(list) ; diffs_T = collections.defaultdict(list)
flds=['cascade_bard_v1','guan_waliser','mundhenk','reid250'] #algorithm names 
diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250']
#loop through the algorithms
for fx,f in enumerate(flds[1:]):
    #loop through the variables 
    print(f'{flds[0]}, {f}')
    for v in vbls:
        #compute the differences between cascasde and all of the others for every variable
        #if fx<len(flds)-1:
        diff = eval(flds[0])[v][0]-eval(f)[v][0]  #compute the difference between the variables for the algorithms 
        eval(f'diffs_{v}')[diff_labels[fx]].append(diff) #store the results as a list 

diff_labels = ['RD250-MDK','RD250-GW']
flds = np.flip(flds)

for fx,f in enumerate(flds[1:-1]):
    #loop through the variables 
    print(f'{flds[0]}, {f}')
    for v in vbls:
        #compute the differences between cascasde and all of the others for every variable
        #if fx<len(flds)-2:
        diff = eval(flds[0])[v][0]-eval(f)[v][0]  #compute the difference between the variables for the algorithms 
        eval(f'diffs_{v}')[diff_labels[fx]].append(diff) #store the results as a list 

diff_labels=['MDK-GW']
flds = np.flip(flds)

for v in vbls:
    #compute the differences between cascasde and all of the others for every variable

    diff = eval(flds[2])[v][0]-eval(flds[1])[v][0]  #compute the difference between the variables for the algorithms 
    eval(f'diffs_{v}')[diff_labels[0]].append(diff) #store the results as a list 


# In[4]:


convert = lambda x: chr(ord('`')+(x+1)) if x<26 else chr(ord('`')+((x)//25))+chr(ord('`')+(x+1)%26)


# In[ ]:





# In[5]:


fig,axes=plt.subplots(ncols=3,nrows=2, figsize=(15,8), subplot_kw={'projection':ccrs.PlateCarree()})
diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']
units = ['kgm$^{-1}$s$^{-1}$' ,'hPa' ,'kgm$^{-2}$' ,'m$^2$s$^{-1}$Kkg$^{-1}$' ,'m' ,'K']

ct=0
for ax,v in zip(axes.flat,diff_labels):
    diffs_IVT[v][0]['IVT'].mean('tstep').plot(ax=ax,cbar_kwargs={'orientation':'horizontal','label':f'IVT ({units[0]})'})
    ax.set_title(f'({convert(ct)}) {v}')
    ct+=1
    cart_plot(ax)


# ## COMPUTE THE STATISICAL SIGNIFICANCE OF THE DIFFERNCE 

# In[5]:


#  mean differences 
diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']
vbls = ['IVT','MSL','TCWV','PV','Z','T']
mean_diffs_IVT = collections.defaultdict(list)
mean_diffs_MSL = collections.defaultdict(list)
mean_diffs_TCWV = collections.defaultdict(list)
mean_diffs_PV = collections.defaultdict(list)
mean_diffs_Z = collections.defaultdict(list)
mean_diffs_T = collections.defaultdict(list)

for vb in (vbls):
    for dl in (diff_labels):
        eval(f'mean_diffs_{vb}')[dl].append(eval(f'diffs_{vb}')[dl][0].mean('tstep'))


# In[12]:


mean_diffs_IVT


# In[6]:


#start the process 
vbls = ['IVT','MSL','TCWV','PV','Z','T']
diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']
sig_IVT = collections.defaultdict(list)
sig_MSL = collections.defaultdict(list)
sig_TCWV = collections.defaultdict(list)
sig_PV = collections.defaultdict(list)
sig_Z = collections.defaultdict(list)
sig_T = collections.defaultdict(list)
conf_lev = 0.99 #confidence level to compute pvalue 


#loop through the variables 
for vb in vbls:
    print(f'Working on {vb}')
    #loop through the differences 
    for dl in diff_labels:
        print(f'Started with {vb}-{dl} ...')
        #keep that data as variable dd 
        mean_data =  eval(f'mean_diffs_{vb}')[dl][0][vb] # mean data for difference in algorithms and for a variable
        significance = []
        #get only lon and lat 
        for lon in mean_data.longitude.values: #loop through lon 
            for lat in mean_data.latitude.values: #loop through lat 
                
                #if data at a lon,lat is > 0 
                if mean_data.sel(longitude=lon,latitude=lat).values > 0:
                    ds = eval(f'diffs_{vb}')[dl][0].sel(longitude=lon,latitude=lat)[vb]#select the corresponding data 
                    num_opposite_sign = np.sum(ds<=0) #find the sum of the bootstrapped data for the point 
                    
                else: 
                    ds = eval(f'diffs_{vb}')[dl][0].sel(longitude=lon,latitude=lat)[vb] #select the corresponding data 
                    num_opposite_sign = np.sum(ds>=0) #find the sum of the bootstrapped data for the point 
                    
                pval = num_opposite_sign/len(ds) #compute pvalue 
                
                #compute confidence interval 
                confidence_level = conf_lev
                
                significance.append(pval)
                #if pval <= (1-confidence_level):
                    #significance.append([lon,lat]) #set 1 for significantly different  
                #elif pval>=(1-confidence_level):
                    #significance.append(0)  #set 0 for not significant different 
        print(f'Done with P-value computation')
        #reshape data 
        
        print('Reshaping significance level binary data ...')
        
        significance = np.array(significance).reshape(len(mean_data.longitude.values),len(mean_data.latitude.values))
        lon,lat = np.meshgrid(list(diffs_IVT['TECA-GW'][0].longitude.values),list(diffs_IVT['TECA-GW'][0].latitude.values))
        ft=xr.DataArray(data=significance,dims=dict(lon=(['longitude','latitude'],lon),lat=(['longitude','latitude'],lat)))

        FT=xr.Dataset({'sig':ft})

        FT.coords['lon']=diffs_IVT['TECA-GW'][0].longitude.values
        FT.coords['lat']=diffs_IVT['TECA-GW'][0].latitude.values
        print(f' Appending significance value to sig_{vb}[{dl}]')
        eval(f'sig_{vb}')[dl].append(FT)


# ### Make plots for the significance as dots and differences as spatial plots

# In[9]:



diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']
units = ['kgm$^{-1}$s$^{-1}$' ,'hPa' ,'kgm$^{-2}$' ,'m$^2$s$^{-1}$Kkg$^{-1}$' ,'m' ,'K']
vbls = ['IVT','MSL','TCWV','PV','Z','T']

text_str = "Significance \n.. 90% \n// 95% \n\\\ 99% "
props = dict(boxstyle='round', facecolor='white', alpha=0.5)

for vx,vb in enumerate(vbls):
    ct=0
    fig,axes=plt.subplots(ncols=3,nrows=2, figsize=(15,8), subplot_kw={'projection':ccrs.PlateCarree()})
    fig.subplots_adjust(bottom=0.15, top = 0.95, left=0.05,hspace=0.05)

    for ax,dl in zip(axes.flat,diff_labels):
        if vx>=3: 
            eval(f'diffs_{vb}')[dl][0][vb].mean('tstep').plot(cmap='bwr',ax=ax,cbar_kwargs={'orientation':'horizontal','label':f' 500hPa {vb} ({units[vx]})'})
        else:
            eval(f'diffs_{vb}')[dl][0][vb].mean('tstep').plot(cmap='bwr',ax=ax,cbar_kwargs={'orientation':'horizontal','label':f'{vb} ({units[vx]})'})

        #xlon =  np.stack(np.array(eval(f'sig_{vb}')[dl][0])[:,0]) ; ylat = np.stack(np.array(eval(f'sig_{vb}')[dl][0])[:,1]) 
        nd=eval(f'sig_{vb}')[dl][0]['sig'].T
        nd.where(nd.values<=0.1).plot.contourf(ax=ax,colors='none',hatches= ['...',None],add_colorbar=False) #90%
        nd.where(nd.values<=0.05).plot.contourf(ax=ax,colors='none',hatches=['//',None], linewidth=0.05,add_colorbar=False) #95%
        nd.where(nd.values<=0.01).plot.contourf(ax=ax,colors='none',hatches=['\\',None],linewidth=0.05,add_colorbar=False) #99%
        
        #plt.text()
        #ax.scatter(xlon,ylat,s=0.001,marker='*',alpha=1,linewidth=0.2,edgecolor ="black")#,linewidth=0.1,alpha=0,edgecolors='grey',hatch='////')
        #plt.scatter()
        
        ax.set_title(f'({convert(ct)}) {dl}')
        ct+=1
        cart_plot(ax)
        fig.text(0.92,0.78,s=text_str,fontweight='light',fontsize=14,bbox=props)
        #fig.text(0.92,0.85,s='// 95% significance',fontweight='light',fontsize=14)
        #fig.text(0.92,0.80,s='\\\ 99% significance',fontweight='light',fontsize=14)
    plt.savefig(f'anomaly-images/combined_significance_{vb}.png', facecolor='white')


# # WORK ON DIFFERENCE IN CONSENSUS AND NON-CONSENSUS

# In[6]:


# load the data into collections 
cascade_bard_v1 = collections.defaultdict(list)
guan_waliser = collections.defaultdict(list)
mundhenk = collections.defaultdict(list)
reid250=collections.defaultdict(list)

diffs_IVT = collections.defaultdict(list) ; diffs_PV = collections.defaultdict(list) 
diffs_MSL = collections.defaultdict(list) ; diffs_Z = collections.defaultdict(list)
diffs_TCWV = collections.defaultdict(list) ; diffs_T = collections.defaultdict(list)

flds=['cascade_bard_v1','guan_waliser','mundhenk','reid250'] #algorithm names 
vbls=['IVT','MSL','TCWV','PV','Z','T'] #variable names 

#loop through algorithms 
for f in flds:
    for v in vbls:
        con_bs = xr.open_dataset(f'Bootstrapped Data/guan_waliser_consensus-boot-{v}.nc') ; con_bs.load()#the boostrapped consensus data is the same for all algorithms
        ds1 = xr.open_dataset(f'Bootstrapped Data/{f}_boot-{v}.nc'); ds1.load()
        store_diff = ds1-con_bs
        eval(f'diffs_{v}')[f].append(store_diff) #compute the means of the data 


# In[7]:


diffs_IVT


# ## SIGNIFICANCE CALCULATION FOR CONSENSUS

# In[ ]:





# In[10]:


#  mean differences 
#diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']

vbls = ['IVT','MSL','TCWV','PV','Z','T']
flds=['cascade_bard_v1','guan_waliser','mundhenk','reid250']

mean_diffs_IVT = collections.defaultdict(list)
mean_diffs_MSL = collections.defaultdict(list)
mean_diffs_TCWV = collections.defaultdict(list)
mean_diffs_PV = collections.defaultdict(list)
mean_diffs_Z = collections.defaultdict(list)
mean_diffs_T = collections.defaultdict(list)

for vb in (vbls):
    for dl in (flds):
        eval(f'mean_diffs_{vb}')[dl].append(eval(f'diffs_{vb}')[dl][0].mean('tstep'))


# In[12]:


#start the process 
vbls = ['IVT','MSL','TCWV','PV','Z','T']
flds=['cascade_bard_v1','guan_waliser','mundhenk','reid250']

diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']
sig_IVT = collections.defaultdict(list)
sig_MSL = collections.defaultdict(list)
sig_TCWV = collections.defaultdict(list)
sig_PV = collections.defaultdict(list)
sig_Z = collections.defaultdict(list)
sig_T = collections.defaultdict(list)
conf_lev = 0.99 #confidence level to compute pvalue 


#loop through the variables 
for vb in vbls:
    print(f'Working on {vb}')
    #loop through the differences 
    for dl in flds:
        print(f'Started with {vb}-{dl} ...')
        #keep that data as variable dd 
        mean_data =  eval(f'mean_diffs_{vb}')[dl][0][vb] # mean data for difference in algorithms and for a variable
        significance = []
        #get only lon and lat 
        for lon in mean_data.longitude.values: #loop through lon 
            for lat in mean_data.latitude.values: #loop through lat 
                
                #if data at a lon,lat is > 0 
                if mean_data.sel(longitude=lon,latitude=lat).values > 0:
                    ds = eval(f'diffs_{vb}')[dl][0].sel(longitude=lon,latitude=lat)[vb]#select the corresponding data 
                    num_opposite_sign = np.sum(ds<=0) #find the sum of the bootstrapped data for the point 
                    
                else: 
                    ds = eval(f'diffs_{vb}')[dl][0].sel(longitude=lon,latitude=lat)[vb] #select the corresponding data 
                    num_opposite_sign = np.sum(ds>=0) #find the sum of the bootstrapped data for the point 
                    
                pval = num_opposite_sign/len(ds) #compute pvalue 
                
                #compute confidence interval 
                confidence_level = conf_lev
                
                significance.append(pval)
                #if pval <= (1-confidence_level):
                    #significance.append([lon,lat]) #set 1 for significantly different  
                #elif pval>=(1-confidence_level):
                    #significance.append(0)  #set 0 for not significant different 
        print(f'Done with P-value computation')
        #reshape data 
        
        print('Reshaping significance level binary data ...')
        
        significance = np.array(significance).reshape(len(mean_data.longitude.values),len(mean_data.latitude.values))
        lon,lat = np.meshgrid(list(diffs_IVT['cascade_bard_v1'][0].longitude.values),list(diffs_IVT['cascade_bard_v1'][0].latitude.values))
        ft=xr.DataArray(data=significance,dims=dict(lon=(['longitude','latitude'],lon),lat=(['longitude','latitude'],lat)))

        FT=xr.Dataset({'sig':ft})

        FT.coords['lon']=diffs_IVT['cascade_bard_v1'][0].longitude.values
        FT.coords['lat']=diffs_IVT['cascade_bard_v1'][0].latitude.values
        print(f' Appending significance value to sig_{vb}[{dl}]')
        eval(f'sig_{vb}')[dl].append(FT)


# ### Make plots for the significance as dots and differences as spatial plots

# In[14]:


convert = lambda x: chr(ord('`')+(x+1)) if x<26 else chr(ord('`')+((x)//25))+chr(ord('`')+(x+1)%26)


# In[19]:



diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']
units = ['kgm$^{-1}$s$^{-1}$' ,'hPa' ,'kgm$^{-2}$' ,'m$^2$s$^{-1}$Kkg$^{-1}$' ,'m' ,'K']
vbls = ['IVT','MSL','TCWV','PV','Z','T']

text_str = "Significance \n.. 90% \n// 95% \n\\\ 99% "
props = dict(boxstyle='round', facecolor='white', alpha=0.5)

for vx,vb in enumerate(vbls):
    ct=0
    fig,axes=plt.subplots(ncols=2,nrows=2, figsize=(15,8), subplot_kw={'projection':ccrs.PlateCarree()})
    fig.subplots_adjust(bottom=0.15, top = 0.95, left=0.05,hspace=0.3)

    for ax,dl in zip(axes.flat,flds):
        if vx>=3: 
            eval(f'diffs_{vb}')[dl][0][vb].mean('tstep').plot(cmap='bwr',ax=ax,cbar_kwargs={'orientation':'horizontal','label':f' 500hPa {vb} ({units[vx]})'})
        else:
            eval(f'diffs_{vb}')[dl][0][vb].mean('tstep').plot(cmap='bwr',ax=ax,cbar_kwargs={'orientation':'horizontal','label':f'{vb} ({units[vx]})'})

        #xlon =  np.stack(np.array(eval(f'sig_{vb}')[dl][0])[:,0]) ; ylat = np.stack(np.array(eval(f'sig_{vb}')[dl][0])[:,1]) 
        nd=eval(f'sig_{vb}')[dl][0]['sig'].T
        nd.where(nd.values<=0.1).plot.contourf(ax=ax,colors='none',hatches= ['...',None],add_colorbar=False) #90%
        nd.where(nd.values<=0.05).plot.contourf(ax=ax,colors='none',hatches=['//',None], linewidth=0.05,add_colorbar=False) #95%
        nd.where(nd.values<=0.01).plot.contourf(ax=ax,colors='none',hatches=['\\',None],linewidth=0.05,add_colorbar=False) #99%
        
        #plt.text()
        #ax.scatter(xlon,ylat,s=0.001,marker='*',alpha=1,linewidth=0.2,edgecolor ="black")#,linewidth=0.1,alpha=0,edgecolors='grey',hatch='////')
        #plt.scatter()
        
        ax.set_title(f'({convert(ct)}) {dl}')
        ct+=1
        cart_plot(ax)
        fig.text(0.92,0.8,s=text_str,fontweight='light',fontsize=14,bbox=props)
        #fig.text(0.92,0.85,s='// 95% significance',fontweight='light',fontsize=14)
        #fig.text(0.92,0.80,s='\\\ 99% significance',fontweight='light',fontsize=14)
    plt.savefig(f'anomaly-images/consensus-combined_significance_{vb}.png', facecolor='white')


# # DIFF BETWEEN CASCADE AND OTHER ALGORITHMS 

# In[3]:


# load the data into collections 
cascade_bard_v1 = collections.defaultdict(list)
guan_waliser = collections.defaultdict(list)
mundhenk = collections.defaultdict(list)
reid250=collections.defaultdict(list)

diffs_IVT = collections.defaultdict(list) ; diffs_PV = collections.defaultdict(list) 
diffs_MSL = collections.defaultdict(list) ; diffs_Z = collections.defaultdict(list)
diffs_TCWV = collections.defaultdict(list) ; diffs_T = collections.defaultdict(list)

flds=['guan_waliser','mundhenk','reid250'] #algorithm names 
vbls=['IVT','MSL','TCWV','PV','Z','T'] #variable names 

#loop through algorithms 
for f in flds:
    for v in vbls:
        con_bs = xr.open_dataset(f'Bootstrapped Data/cascade_bard_v1_boot-{v}.nc') ; con_bs.load()#the boostrapped consensus data is the same for all algorithms
        ds1 = xr.open_dataset(f'Bootstrapped Data/{f}_boot-{v}.nc'); ds1.load()
        store_diff = ds1-con_bs
        eval(f'diffs_{v}')[f].append(store_diff) #compute the means of the data 


# In[4]:


#  mean differences 
#diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']

vbls = ['IVT','MSL','TCWV','PV','Z','T']
flds=['guan_waliser','mundhenk','reid250']

mean_diffs_IVT = collections.defaultdict(list)
mean_diffs_MSL = collections.defaultdict(list)
mean_diffs_TCWV = collections.defaultdict(list)
mean_diffs_PV = collections.defaultdict(list)
mean_diffs_Z = collections.defaultdict(list)
mean_diffs_T = collections.defaultdict(list)

for vb in (vbls):
    for dl in (flds):
        eval(f'mean_diffs_{vb}')[dl].append(eval(f'diffs_{vb}')[dl][0].mean('tstep'))


# In[7]:


#start the process 
vbls = ['IVT','MSL','TCWV','PV','Z','T']
flds=['guan_waliser','mundhenk','reid250']

diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']
sig_IVT = collections.defaultdict(list)
sig_MSL = collections.defaultdict(list)
sig_TCWV = collections.defaultdict(list)
sig_PV = collections.defaultdict(list)
sig_Z = collections.defaultdict(list)
sig_T = collections.defaultdict(list)
conf_lev = 0.99 #confidence level to compute pvalue 


#loop through the variables 
for vb in vbls:
    print(f'Working on {vb}')
    #loop through the differences 
    for dl in flds:
        print(f'Started with {vb}-{dl} ...')
        #keep that data as variable dd 
        mean_data =  eval(f'mean_diffs_{vb}')[dl][0][vb] # mean data for difference in algorithms and for a variable
        significance = []
        #get only lon and lat 
        for lon in mean_data.longitude.values: #loop through lon 
            for lat in mean_data.latitude.values: #loop through lat 
                
                #if data at a lon,lat is > 0 
                if mean_data.sel(longitude=lon,latitude=lat).values > 0:
                    ds = eval(f'diffs_{vb}')[dl][0].sel(longitude=lon,latitude=lat)[vb]#select the corresponding data 
                    num_opposite_sign = np.sum(ds<=0) #find the sum of the bootstrapped data for the point 
                    
                else: 
                    ds = eval(f'diffs_{vb}')[dl][0].sel(longitude=lon,latitude=lat)[vb] #select the corresponding data 
                    num_opposite_sign = np.sum(ds>=0) #find the sum of the bootstrapped data for the point 
                    
                pval = num_opposite_sign/len(ds) #compute pvalue 
                
                #compute confidence interval 
                confidence_level = conf_lev
                
                significance.append(pval)
                #if pval <= (1-confidence_level):
                    #significance.append([lon,lat]) #set 1 for significantly different  
                #elif pval>=(1-confidence_level):
                    #significance.append(0)  #set 0 for not significant different 
        print(f'Done with P-value computation')
        #reshape data 
        
        print('Reshaping significance level binary data ...')
        
        significance = np.array(significance).reshape(len(mean_data.longitude.values),len(mean_data.latitude.values))
        lon,lat = np.meshgrid(list(diffs_IVT['guan_waliser'][0].longitude.values),list(diffs_IVT['guan_waliser'][0].latitude.values))
        ft=xr.DataArray(data=significance,dims=dict(lon=(['longitude','latitude'],lon),lat=(['longitude','latitude'],lat)))

        FT=xr.Dataset({'sig':ft})

        FT.coords['lon']=diffs_IVT['guan_waliser'][0].longitude.values
        FT.coords['lat']=diffs_IVT['guan_waliser'][0].latitude.values
        print(f' Appending significance value to sig_{vb}[{dl}]')
        eval(f'sig_{vb}')[dl].append(FT)


# In[8]:


convert = lambda x: chr(ord('`')+(x+1)) if x<26 else chr(ord('`')+((x)//25))+chr(ord('`')+(x+1)%26)


# In[12]:



diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']
units = ['kgm$^{-1}$s$^{-1}$' ,'hPa' ,'kgm$^{-2}$' ,'m$^2$s$^{-1}$Kkg$^{-1}$' ,'m' ,'K']
vbls = ['IVT','MSL','TCWV','PV','Z','T']

text_str = "Significance \n.. 90% \n// 95% \n\\\ 99% "
props = dict(boxstyle='round', facecolor='white', alpha=0.5)

for vx,vb in enumerate(vbls):
    ct=0
    fig,axes=plt.subplots(ncols=3, figsize=(18,4), subplot_kw={'projection':ccrs.PlateCarree()})
    fig.subplots_adjust(bottom=0.15, top = 0.95, left=0.05,hspace=0.3)

    for ax,dl in zip(axes.flat,flds):
        if vx>=3: 
            eval(f'diffs_{vb}')[dl][0][vb].mean('tstep').plot(cmap='bwr',ax=ax,cbar_kwargs={'orientation':'horizontal','label':f' 500hPa {vb} ({units[vx]})'})
        else:
            eval(f'diffs_{vb}')[dl][0][vb].mean('tstep').plot(cmap='bwr',ax=ax,cbar_kwargs={'orientation':'horizontal','label':f'{vb} ({units[vx]})'})

        #xlon =  np.stack(np.array(eval(f'sig_{vb}')[dl][0])[:,0]) ; ylat = np.stack(np.array(eval(f'sig_{vb}')[dl][0])[:,1]) 
        nd=eval(f'sig_{vb}')[dl][0]['sig'].T
        nd.where(nd.values<=0.1).plot.contourf(ax=ax,colors='none',hatches= ['...',None],add_colorbar=False) #90%
        nd.where(nd.values<=0.05).plot.contourf(ax=ax,colors='none',hatches=['//',None], linewidth=0.05,add_colorbar=False) #95%
        nd.where(nd.values<=0.01).plot.contourf(ax=ax,colors='none',hatches=['\\',None],linewidth=0.05,add_colorbar=False) #99%
        
        #plt.text()
        #ax.scatter(xlon,ylat,s=0.001,marker='*',alpha=1,linewidth=0.2,edgecolor ="black")#,linewidth=0.1,alpha=0,edgecolors='grey',hatch='////')
        #plt.scatter()
        
        ax.set_title(f'({convert(ct)}) {dl}')
        ct+=1
        cart_plot(ax)
        fig.text(0.92,0.7,s=text_str,fontweight='light',fontsize=14,bbox=props)
        #fig.text(0.92,0.85,s='// 95% significance',fontweight='light',fontsize=14)
        #fig.text(0.92,0.80,s='\\\ 99% significance',fontweight='light',fontsize=14)
    plt.savefig(f'anomaly-images/teca-combined_significance_{vb}.png', facecolor='white')


# In[ ]:





# In[ ]:





# In[22]:


xr.open_dataset('Bootstrapped Data/guan_waliser_consensus-boot-IVT.nc')


# In[ ]:





# In[ ]:


# calculate the difference between the two points
mean_diff = ... # this is a scalar

# get the bootstrap samples of the differences
bootstrap_diff = ... # this is a numpy array

# calculate the probability of the null hypothesis (the mean difference is 0)
if mean_diff > 0:
  num_opposite_sign = np.sum(bootstrap_diff <= 0)
else:
  num_opposite_sign = np.sum(bootstrap_diff >= 0)

p_val = num_opposite_sign / len(bootstrap_diff)

# test at the 95% confidence level
confidence_level = 0.95
if pval <= (1 - confidence_level):
  # reject the null hypothesis that they are from the same distribution
  print(f"the probability of the null hypothesis being true is less than {100*(1 - confidence_level)}%")


# In[ ]:





# In[22]:


import seaborn as sb


# In[23]:


sb.kdeplot(diffs_IVT['TECA-GW'][0]['IVT'].sel(longitude=-122,latitude=37).values,label='TECA-GW')
sb.kdeplot(diffs_IVT['TECA-MDK'][0]['IVT'].sel(longitude=-122,latitude=37).values,label='TECA-MDK')
sb.kdeplot(diffs_IVT['TECA-RD250'][0]['IVT'].sel(longitude=-122,latitude=37).values,label='TECA-RD250')
#sb.kdeplot(diffs_IVT['MDK-GW'][0]['IVT'].sel(longitude=-123,latitude=37).values)
sb.kdeplot(diffs_IVT['MDK-GW'][0]['IVT'].sel(longitude=-122,latitude=37).values,label='MDK-GW')

sb.kdeplot(diffs_IVT['RD250-GW'][0]['IVT'].sel(longitude=-122,latitude=37).values,label='RD250-GW')
sb.kdeplot(diffs_IVT['RD250-MDK'][0]['IVT'].sel(longitude=-122,latitude=37).values,label='RD250-MDK')
plt.legend()


# ### PDF DISTRIBUTION PLOT FOR A POINT

# In[24]:


diff_labels = ['TECA-GW','TECA-MDK','TECA-RD250','RD250-MDK','RD250-GW','MDK-GW']
colors = ['blue','black','red','green','purple','orange']
vbls = ['IVT', 'MSL', 'TCWV', 'PV', 'Z', 'T']
units = ['kgm$^{-1}$s$^{-1}$' ,'hPa' ,'kgm$^{-2}$' ,'m$^2$s$^{-1}$Kkg$^{-1}$' ,'m' ,'K']

lon=-122 ; lat = 37
fig,axes= plt.subplots(figsize=(18,10),ncols=3,nrows=2)
fig.subplots_adjust(wspace=0.25,bottom=0.1, top = 0.95, left=0.1,hspace=0.25)


ct=0
for vb,ax in zip(vbls,axes.flat):
    for lx,lb in enumerate(diff_labels):
        if ct <3:
            sb.kdeplot( eval( f'diffs_{vb}' )[lb][0][vb].sel( longitude=lon , latitude=lat ).values, color=colors[lx],ax = ax)
        elif ct>=3 : 
            sb.kdeplot( eval( f'diffs_{vb}' )[lb][0][vb].sel( longitude=lon , latitude=lat ).mean('level').values, color=colors[lx],ax = ax)
        ax.set_title(vb,fontsize=18)
        ax.set_xlabel(f' {vb} ({units[ct]})',fontsize=14)
        ax.set_ylabel('Density',fontsize=14)
    ct+=1
    if ct == 5:
        fig.legend(labels=diff_labels,ncol=2,loc='lower right')
plt.savefig('anomaly-images/PDF-Bootstrapped.png',facecolor='white')


# # PLOT FOR CONSENSUS TIMES AND NCTS ARDTs 

# In[99]:


import cmocean as cmo


# In[100]:


main_path = '/global/cscratch1/sd/kquagra/ARs Work/larger_region_AR_capture/' #main path
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
#ARDTs = single_data(main_path,flds,'landfall_ARs__0.9-0.4.nc')
ARDTs = ARDTs_load(main_path,flds,'newlandfall_ARs__0.8-1.0') #load main AR data c


# In[101]:


#SELECT DAYS THAT APPEAR IN ALL THE ALGORITHMS 
consensus_days = reduce(np.intersect1d,(ARDTs['guan_waliser'][0].time.values,ARDTs['cascade_bard_v1'][0].time.values,ARDTs['reid250'][0].time.values,ARDTs['mundhenk'][0].time.values))


# In[102]:


cons_ARDTs =   { } # specific consensus days when ARs are captured 
non_cons_ARDTs ={ }# select the specific non_consensus days for AR captured 

for fx,f in enumerate(flds):
    f = f[1:-1]
    cons_ARDTs[f]=[] ;non_cons_ARDTs[f]=[] #empty list to hold data 
    cons_ARDTs[f].append(ARDTs[f][0].sel(time=np.isin(ARDTs[f][0].time.values,consensus_days))) #compute the consensus days 
    non_cons_ARDTs[f].append(ARDTs[f][0].sel(time = np.isin(ARDTs[f][0].time.values,consensus_days,invert=True))) #invert to get non-consensus days 


# In[110]:


fig,axes= plt.subplots(ncols=2,nrows=2,figsize=(15,10),subplot_kw={'projection':ccrs.PlateCarree()})
title=['TECA_BARD v1','Guan & Waliser','Reid 250','Mundhenk v3']

ct=0
for ax,f in zip(axes.flat,non_cons_ARDTs.keys()):
    non_cons_ARDTs[f][0]['ar_binary_tag'].sum('time').plot.contourf(cmap=cmo.cm.matter,vmax=1500,levels=np.arange(0,2000,50),ax=ax,cbar_kwargs={'orientation':'horizontal','label':'AR count'})
    pb=non_cons_ARDTs[f][0]['ar_binary_tag'].sum('time').plot.contour(cmap='black',vmin=80,ax=ax, add_colorbar=False)#,cbar_kwargs={'orientation':'horizontal','label':'AR binary tag'})
    pb.clabel(inline=True)
    frac = (len(non_cons_ARDTs[f][0]['ar_binary_tag'].values)/(len(cons_ARDTs[f][0]['ar_binary_tag'].values)+len(non_cons_ARDTs[f][0]['ar_binary_tag'].values)))*100
    ax.set_title(f'({convert(ct)}) {title[ct]} \n[AR fraction = {round(frac,2)} %]')
    cart_plot(ax)
    ct+=1
plt.savefig('Images/Final_Images/Non-consensus.png',facecolor='white')


# In[109]:


fig,axes= plt.subplots(ncols=2,nrows=2,figsize=(15,10),subplot_kw={'projection':ccrs.PlateCarree()})
title=['TECA_BARD v1','Guan & Waliser','Reid 250','Mundhenk v3']

ct=0
for ax,f in zip(axes.flat,cons_ARDTs.keys()):
    cons_ARDTs[f][0]['ar_binary_tag'].sum('time').plot.contourf(cmap=cmo.cm.matter,vmax=500,levels=np.arange(0,500,20),ax=ax,cbar_kwargs={'orientation':'horizontal','label':'AR count'})
    pb=cons_ARDTs[f][0]['ar_binary_tag'].sum('time').plot.contour(cmap='black',vmin=200,levels=np.arange(200,300,40),ax=ax, add_colorbar=False)#,cbar_kwargs={'orientation':'horizontal','label':'AR binary tag'})
    pb.clabel(inline=True)
    frac = (len(cons_ARDTs[f][0]['ar_binary_tag'].values)/(len(cons_ARDTs[f][0]['ar_binary_tag'].values)+len(non_cons_ARDTs[f][0]['ar_binary_tag'].values)))*100
    ax.set_title(f'({convert(ct)}) {title[ct]} \n[AR fraction = {round(frac,2)} %]')
    cart_plot(ax)
    ct+=1
plt.savefig('Images/Final_Images/Consensus.png',facecolor='white')


# In[106]:


(len(cons_ARDTs[f][0]['ar_binary_tag'].values))/


# # PLOTS FOR VARIABLES CONSENSUS AND NON- CONSENSUS

# ### PLOT FOR CONSENSUS DAYS

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# ## COMPUTE CONSENSUS AND NON-CONSENSUS DAYS FOR IVT

# In[2]:


main_path = '/global/cscratch1/sd/kquagra/ARs Work/larger_region_AR_capture/' #main path
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
#ARDTs = single_data(main_path,flds,'landfall_ARs__0.9-0.4.nc')
ARDTs = ARDTs_load(main_path,flds,'newlandfall_ARs__0.8-1.0') #load main AR data c
#SELECT DAYS THAT APPEAR IN ALL THE ALGORITHMS 
consensus_days = reduce(np.intersect1d,(ARDTs['guan_waliser'][0].time.values,ARDTs['cascade_bard_v1'][0].time.values,ARDTs['reid250'][0].time.values,ARDTs['mundhenk'][0].time.values))


# In[3]:


cons_ARDTs =   { } # specific consensus days when ARs are captured 
non_cons_ARDTs ={ }# select the specific non_consensus days for AR captured 

for fx,f in enumerate(flds):
    f = f[1:-1]
    cons_ARDTs[f]=[] ;non_cons_ARDTs[f]=[] #empty list to hold data 
    cons_ARDTs[f].append(ARDTs[f][0].sel(time=np.isin(ARDTs[f][0].time.values,consensus_days))) #compute the consensus days 
    non_cons_ARDTs[f].append(ARDTs[f][0].sel(time = np.isin(ARDTs[f][0].time.values,consensus_days,invert=True))) #invert to get non-consensus days 


# In[4]:


#Load IVT DATA
ivt_data = xr.open_mfdataset(glob.glob('../../era5_IVT/*.nc',recursive=True))


# In[20]:


non_cascade_bard_v1_ivt = collections.defaultdict(list)
non_guan_waliser_ivt = collections.defaultdict(list)
non_mundhenk_ivt = collections.defaultdict(list)
non_reid250_ivt = collections.defaultdict(list)
flds = ['/cascade_bard_v1/','/guan_waliser/','/mundhenk/','/reid250/'] #algorithms 
vb='IVT'

for fx,f in enumerate(flds):
    f = f[1:-1]
    
    newd = ivt_data.sel(time=cons_ARDTs[f][0].time.values)
    eval(f'non_{f}_ivt')['IVT'].append(newd)
    
    newd = newd.chunk({'time':1})
    SAVE_DIR='../Final_variable_data/IVT_consensus_point/'
    
    newd.to_netcdf(f'{SAVE_DIR}{f}True_concensus_{vb}.nc')
    


# In[14]:


newd


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[28]:


from xrrandom import stats
from xrrandom import bootstrap_samples


# In[2]:


exec(open('imports.py').read())

#data paths 
clim_data_path = ['/global/cscratch1/sd/kquagra/era5_IVT/',
                  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/msl_climatologies/',
                  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/tcwv_climatologies/',
                  '/global/cscratch1/sd/kquagra/new-ERA5_PV/',
                  '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/z_climatologies/',
                 '/global/cscratch1/sd/kquagra/ARs Work/Non-coinciding ARDTs/Climatologies/temp_climatologies/',
                 ]
#variable names 
vbls = ['ivt','msl','tcwv','pv','z','t']

#files paths 
files = [glob.glob(f'{clim_data_path[vx]}*.nc',recursive=True) for vx,v in enumerate(vbls)]
#files.sort()


#iterate through dataset and obtain the climatology as a dictionary
climatology=collections.defaultdict(list)
vbls = ['ivt','msl','tcwv','pv','z','t']
for px,pts in enumerate(clim_data_path):
    dd =[]
    
    #check if the file name is for msl,tcwv,z or temp since they have a different arrangement 
    if 0<px <=2 or 3<px<=5: 
        for yy in range(1980,2018): #loop through the files year by year 
            dim='month'  #set the concatenation dimension 
            tm=pd.date_range(f'{yy}-01',f'{yy+1}-01',freq='1M') #create the time axis 
            tm = [np.datetime64(x) for x in tm]  #convert to datetime64
            ds = xr.open_mfdataset(glob.glob(f'{pts}*{yy}*.nc',recursive=True),parallel=True)  #open the dataset

            ds['month'] =tm  #reassign the time 
            ds=ds.sel(month=ds.month.dt.month.isin([1,2,12]))  #select only DJF
            ds.coords['longitude'] = (ds.coords['longitude']  + 180) % 360 - 180
            ds=ds.sortby(ds.longitude)
            ds = ds.sel(longitude=slice(-170,-60),latitude=slice(55,10))
            dd.append(ds.mean('month')) #append to the dd list 
        new_dd = xr.concat(dd,dim=dim)  #concatenate the list of yearly datasets into one file
    elif px==0 or px ==3:
        for yy in range(1980,2018):
            dim='time'
            ds = xr.open_mfdataset(glob.glob(f'{pts}*{yy}*.nc',recursive=True),parallel=True)
            ds = ds.sel(time=ds.time.dt.month.isin([1,2,12]))
            ds.coords['longitude'] = (ds.coords['longitude']  + 180) % 360 - 180
            ds=ds.sortby(ds.longitude)
            ds = ds.sel(longitude=slice(-170,-60),latitude=slice(55,10))
            
            dd.append(ds.mean('time'))
        new_dd = xr.concat(dd,dim=dim)
    climatology[vbls[px]].append(new_dd.mean(dim)) #save the variable DJF climatology 
    
    
    
    
non_con_paths = []
flds = ['cascade_bard_v1','guan_waliser','reid250','mundhenk']
vbls=['IVT','MSL','TCWV','PV','GEOPTH','TEMP']
cascade_bard_v1 = collections.defaultdict(list)
guan_waliser = collections.defaultdict(list)
reid250 = collections.defaultdict(list)
mundhenk = collections.defaultdict(list)

for cx,c in enumerate(list(climatology.keys())):
    for fx,f in enumerate(flds):
        clim = climatology[c][0]
        #clim.coords['longitude'] = (clim.coords['longitude']  + 180) % 360 - 180
        #clim=clim.sortby(clim.longitude)
        #clim = clim.load()
        if cx<=2:
            ds = xr.open_dataset(f'{vbls[cx]}_non_consensus_point/{f}False_concensus_{vbls[cx]}.nc')#; ds.load()
        if cx>2: 
            ds = xr.open_dataset(f'{vbls[cx]}_non_consensus_point/{f}False_concensus_{vbls[cx]}.nc').sel(level=500)#; ds.load()
        anom = ds- clim
        eval(f)[c].append(anom)


# In[24]:


cascade_bard_v1_boot = collections.defaultdict(list)
guan_waliser_boot = collections.defaultdict(list)
mundhenk_boot = collections.defaultdict(list)
reid250_boot = collections.defaultdict(list)
vb=['IVT','MSL','TCWV','PV','Z','T']

print('Beginning Bootstrap computations')
for fx,f in enumerate(flds):
    
    for vx,v in enumerate(vb):
        data = []
        for num_times in range(1000):
        
            dtime = eval(f)[v.lower()][0].time.values
            dtim = np.random.choice(dtime,size=len(dtime))
            df = eval(f)[v.lower()][0][v]
            df1=df.sel(time=df.time.dt.time.isin(list(dtim)))#.mean('time')
            df1.coords['tstep'] = num_times
            data.append(df1)
        concated = xr.concat(data,dim='tstep')
        #eval(f'{f}_boot')[v].append(concated)
        concated.to_netcdf(f'Bootstrapped Data/{f}_boot-{v}.nc')
print('Done with all procedures')


# In[26]:


df.sel(time=dtim)


# In[29]:


nd = bootstrap_samples(cascade_bard_v1['ivt'][0]['IVT'],'time',100)


# In[44]:


nd=nd.chunk({'sample':1,'time':1})


# In[48]:


nd[0].mean('time')


# In[6]:


nd = nd.chunk({'sample':1})


# In[ ]:


nd[0][0].plot()


# In[4]:


cascade_bard_v1


# In[ ]:
'''



