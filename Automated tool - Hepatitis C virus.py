#!/usr/bin/env python
# coding: utf-8

# Automated data mining of PubChem to fetch required data perform autoQSAR

# In[2]:


import pandas as pd
import numpy as np
import urllib
import os 
import requests
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import seaborn as seabornInstance 
from sklearn.model_selection import train_test_split 
from sklearn.linear_model import LinearRegression
from sklearn import metrics
get_ipython().magic(u'matplotlib inline')
import time


# In[3]:


link1 = "https://pubchem.ncbi.nlm.nih.gov/sdq/cgi2rcgi.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22bioactivity%22,%22where%22:{%22ands%22:[{%22protacxn%22:%22notnull%22},{%22cid%22:%22notnull%22},{%22repacxn%22:%22ABN45876%22}]},%22order%22:[%22acvalue,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22PROTACXN_ABN45876_bioactivity_protein%22}"


# In[4]:


for i in range(0,2):
    try:
        os.remove('downloaded1.csv')
        #print("Deleted old File")
        break
    except Exception as e:
        #print("No file ")
        break
    else:
        break


# In[125]:


for i in range(2):
    try:
        data = pd.read_csv(link1)
        break
    #except IncompleteRead as I:
     #   print("Server Overloading , Proceeding")
      #  break
    except Exception as a:
        print(str(a)+" is the error , Trying {} time".format(i))
        continue
    else:
        break
else:
    print("something Wrong , Try running Again [refer error code for more]")


# In[7]:


data = pd.DataFrame(data)


# In[229]:


new_data = data[['cid','acvalue']]
new_data = new_data.dropna()
new_data


# In[9]:


except_val=0
cid_value = new_data['cid'].to_list()
PIC50_value = (-np.log10(new_data['acvalue']*10**-6)).to_list()

PIC50_value = pd.DataFrame(PIC50_value,columns = ["y"])
y_data = PIC50_value


# In[10]:


link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/new_link_me/property/MolecularWeight,HeavyAtomCount,XLOGP,Complexity,HBondAcceptorCount,MonoisotopicMass,RotatableBondCount,TPSA/CSV"
link_fixed = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/new_link_me/property/MolecularWeight,HeavyAtomCount,XLOGP,Complexity,HBondAcceptorCount,MonoisotopicMass,RotatableBondCount,TPSA/CSV"
sub_link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/replaceme/cids/TXT"
sub_link_fixed = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/cid/replaceme/cids/TXT"
#the link with the required properties.
counter = 1


# In[11]:


final_data_frame = pd.DataFrame()
x_data = pd.DataFrame()
cid_main_counter=0 
for x in cid_value:
    for iter_x in range(1000):
        try:    
            link = link_fixed
            link = link.replace("new_link_me",str(x))
            data1 = pd.read_csv(link)
            data1 = pd.DataFrame(data1)
            x_data = x_data.append(data1)
            data1 = 0
        except Exception as e:
            print("Exception Encountered as {} .! , Trying again !Iteration : {}".format(str(e),str(iter_x)))
            time.sleep(5)
            continue
        else:
            break
    
    else:
        print("Something Wrong with the Trials ! Restart The Algorithm!")


# In[12]:


x_data_saved = x_data
x_data = x_data_saved
n_cid=[]


# In[13]:


#20 sec Wait_time
for i in cid_value:
    max_tries = 10
    
    for iter_ in range(max_tries):
        try:
            #print("inside Try")
            #time.sleep(5)


            link = link_fixed
            link = link.replace("new_link_me",str(i))
            f_data = pd.read_csv(link) 
            
            sub_link = sub_link_fixed
            sub_link = sub_link.replace("replaceme",str(i))
            res = urllib.request.urlopen(sub_link)
            data_sub = res.read()
            data_sub = str(data_sub)
            data_sub = data_sub.replace("\\n",",") ; data_sub = data_sub.replace('b'," ")
            data_sub = data_sub.replace("'",""); data_sub = data_sub.replace(" ","")
            n_count_loop = 0 
            #time.sleep(3)
            #print("Going for the loop")
            for j in data_sub.split(","):
                n_count_loop+=1
                if(len(j)>1):
                    if(n_count_loop<=30):
                        n_cid.append(j)
                    else:
                        break
                else:
                    pass
             
        except:
            if(except_val>=15):
                break
            else:
                print("Re-Trying")
                except_val+=1
                #time.sleep(3)
                continue
        else:
            break
    '''else:
        print("Either The Network is Down(So no point in continuing ) , Or Some Uknown Error Spotted!. Refer Error Code , Continuing with Fetched Data")
        time.sleep(5)
        break'''
    
    
    '''final_data_frame = final_data_frame.append(f_data)
    cid_main_counter+=1
    final_data_frame.to_csv("fdf1.csv")
         
   

    print("Main cid no: {}".format(cid_main_counter))
    #time.sleep(5)'''
    


# In[14]:


k_count = 0
t_count = 0
phase_count=0
f_data = pd.DataFrame()


# In[15]:


for k in n_cid:
    k_count+=1
    t_count+=1
    phase_count+=1
    for sub_iter in range(10):
            try:
                if(t_count>=50):
                    #print("taking Rest")
                    #time.sleep(15)
                    t_count=0
                    
                #print("Inside Try"+str(k_count))
            
                link = link_fixed
                link = link.replace("new_link_me",str(k))
                #print("replacement done"+str(k_count))
                f_data_df = pd.read_csv(link)
                f_data = f_data.append(f_data_df)
                print("Passed Without Exception"+str(k_count))
                break

            except Exception as e:
                print(str(e) + "Encountered , Please Wait :+ "+str(k_count))
                #time.sleep(20)
                continue

    final_data_frame = final_data_frame.append(f_data)
    cid_main_counter+=1
    final_data_frame.to_csv("fdf1.csv")
    print("Data Fetching Continued"+ str(k_count))
    
    


# In[272]:


final_data_frame.drop_duplicates(inplace=True)
print(final_data_frame.shape)
print("The Final Training DataSet: X_data = {} , Y_data = {}".format(x_data.shape,y_data.shape))
print("The Final Data For Preditcion is : {}".format(final_data_frame.shape))


# In[273]:


try:
    final_data_frame.reset_index(inplace=True)
    final_data_frame.drop('index',axis=1,inplace=True)
except:
    pass


# In[274]:


#DON'T RE-RUN AT ANY COST
final_data_frame_reset = final_data_frame
x_data_reset = x_data
y_data_reset = y_data


# In[275]:


#fghjk


#  

#  

# Machine Learning based AutoQSAR to predict drug leads 

# In[276]:


x_data = pd.DataFrame(x_data)


# In[277]:


y_data = pd.DataFrame(y_data)


# In[278]:


print("X :- "+str(x_data.shape)+" Y :- "+str(y_data.shape) )


# In[279]:


train_file = pd.DataFrame()
x_data = x_data.astype("float64")
y_data = y_data.astype("float64")


# In[280]:


x_data.head()


# In[281]:


y_data.head()


# In[282]:


x_data.head()


# In[283]:


new = y_data['y'].to_list()


# In[284]:


train_file  = x_data
train_file['y'] = new


# In[285]:


for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')
train_file


# In[286]:


x = list(train_file.columns)
x = x[:-1]

x_ = train_file[x].values
y_ = train_file['y'].values


# In[287]:


train_file.describe()
x_ = x_.astype("float64")
y_ = y_.astype("float64")


# In[288]:


X_train, X_test, y_train, y_test = train_test_split(x_, y_, test_size=0.1, random_state=0)


#  

# In[289]:


from sklearn.preprocessing import PolynomialFeatures
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing


# In[290]:


#poly regression - degree 1
scaler = preprocessing.StandardScaler()
degree=1
polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
polyreg_scaled.fit(X_train,y_train)


# In[291]:



y_pred = polyreg_scaled.predict(X_test)
main_r2 = r2_score(y_test, y_pred, multioutput='uniform_average')
max_ = main_r2 


# In[292]:


from itertools import combinations
comb_list = [[]]


# In[293]:


#creating all possible combinations of descriptors
def sub(arr,r):
    global comb_list
    for i in r:
        comb = list(combinations(arr,i))
        comb_list.append(comb)
    return comb_list
newone = 0
newone = sub(x , [2,3,4,5,6])
del newone[0]


# In[294]:


#Linear regression based QSAR model 
coef_dict  = 0
loop_index = 0
coef_dict = [[]]
r2_score_new =[]
max_r2 = []
index_r2 =[]
for i in range(0,len(newone)):
    
    index_r=0
    for combi_ in newone[loop_index]:
        
        #print(combi_)
        features = list(combi_)
        features_=train_file[features].values
        output_=train_file['y'].values
        X_train, X_test, y_train, y_test = train_test_split(features_, output_, test_size=0.1, random_state=0)
        scaler = preprocessing.StandardScaler()
        degree=1
        polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
        polyreg_scaled.fit(X_train,y_train)
        y_pred = polyreg_scaled.predict(X_test)
        r2_score_  = r2_score(y_test, y_pred, multioutput='uniform_average')
        r2_score_new.append(r2_score_)
    loop_index+=1
        
    
    max_r2.append(max(r2_score_new))
    index_r = r2_score_new.index(max(r2_score_new))
    index_r2.append(index_r)
    r2_score_new =[]


# In[295]:


#Extracting the best model with highest Rsquare
sec_index = max_r2.index(max(max_r2))
fir_index = index_r2[sec_index]


# In[296]:


if main_r2 > max(max_r2):
    r2_features = x
else:
    
    features  = newone[sec_index][fir_index]
    maxi_r2 = max(max_r2)


# In[297]:


model_dict = {}
model_dict[1]=[maxi_r2,features]


# In[298]:


for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')


# In[299]:


#Building the non-linear(degree 2) regression based QSAR model
coef_dict  = 0
loop_index = 0
coef_dict = [[]]
r2_score_new =[]
max_r2 = []
index_r2 =[]
for i in range(0,len(newone)):
    
    index_r=0
    for combi_ in newone[loop_index]:
        
        #print(combi_)
        features = list(combi_)
        features_=train_file[features].values
        output_=train_file['y'].values
        X_train, X_test, y_train, y_test = train_test_split(features_, output_, test_size=0.1, random_state=0)
        scaler = preprocessing.StandardScaler()
        degree=2
        polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
        polyreg_scaled.fit(X_train,y_train)
        y_pred = polyreg_scaled.predict(X_test)
        r2_score_  = r2_score(y_test, y_pred, multioutput='uniform_average')
        r2_score_new.append(r2_score_)
    loop_index+=1
        
    
    max_r2.append(max(r2_score_new))
    index_r = r2_score_new.index(max(r2_score_new))
    index_r2.append(index_r)
    r2_score_new =[]


# In[300]:


#Extracting the best model with highest Rsquare
sec_index = max_r2.index(max(max_r2))
fir_index = index_r2[sec_index]


# In[301]:


if main_r2 > max(max_r2):
    r2_features = x
else:
    
    features  = newone[sec_index][fir_index]
    maxi_r2 = max(max_r2)


# In[302]:


model_dict[2]=[maxi_r2,features]


# In[303]:


for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')


# In[304]:


model_dict


# In[305]:


#Building the non-linear(degree 2) regression based QSAR model
coef_dict  = 0
loop_index = 0
coef_dict = [[]]
r2_score_new =[]
max_r2 = []
index_r2 =[]
for i in range(0,len(newone)):
    
    index_r=0
    for combi_ in newone[loop_index]:
        
        #print(combi_)
        features = list(combi_)
        features_=train_file[features].values
        output_=train_file['y'].values
        X_train, X_test, y_train, y_test = train_test_split(features_, output_, test_size=0.1, random_state=0)
        scaler = preprocessing.StandardScaler()
        degree=3
        polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
        polyreg_scaled.fit(X_train,y_train)
        y_pred = polyreg_scaled.predict(X_test)
        r2_score_  = r2_score(y_test, y_pred, multioutput='uniform_average')
        r2_score_new.append(r2_score_)
    loop_index+=1
        
    
    max_r2.append(max(r2_score_new))
    index_r = r2_score_new.index(max(r2_score_new))
    index_r2.append(index_r)
    r2_score_new =[]


# In[306]:


#Extracting the best model with highest Rsquare
sec_index = max_r2.index(max(max_r2))
fir_index = index_r2[sec_index]


# In[307]:


if main_r2 > max(max_r2):
    r2_features = x
else:
    
    features  = newone[sec_index][fir_index]
    maxi_r2 = max(max_r2)


# In[308]:


model_dict[3]=[maxi_r2,features]


# In[309]:


for x in train_file.isnull().any():
    if x == True:
        train_file = train_file.fillna(method='ffill')


# In[310]:


model_dict


# In[311]:


top_cid_mix_df = pd.DataFrame()


#   

#  

#  

# Drug lead predictions from the linear and non-linear based QSAR models corresponding to highest Rsquare  

# In[312]:


model_dict[1]

features = list(model_dict[1][1])

x_trained =train_file[features].values
y_trained = train_file['y'].values
x_trained.shape,y_trained.shape


# In[313]:


X_train, X_test, y_train, y_test = train_test_split(x_trained, y_trained, test_size=0.1, random_state=0)


# In[314]:


#poly regression - degree 1
scaler = preprocessing.StandardScaler()
degree=1
polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
polyreg_scaled.fit(X_train,y_train)


# In[315]:


y_pred = polyreg_scaled.predict(X_test)
r2_score(y_test, y_pred, multioutput='uniform_average')

final_data_frameF1 = final_data_frame[final_data_frame['MolecularWeight'] <= 500]
final_data_frameF2 = final_data_frameF1[final_data_frameF1['XLogP'] <=5.6]


# In[316]:


pred_data_cid = final_data_frameF2["CID"].to_list()


# In[317]:


#features.append('CID')
new_features = features.append("CID")
features


# In[318]:


final_data_frame1 = final_data_frameF2[features]
for x in final_data_frame1.isnull().any():
    if x == True:
        final_data_frame1 = final_data_frame1.fillna(method='ffill')
final_data_frame1.drop(columns='CID',inplace=True)
final_data_frame1.shape


# In[319]:


final_pred_1 = polyreg_scaled.predict(final_data_frame1)


# In[320]:


final_pred_1 = list(final_pred_1)
final_data_frame1['CID'] = pred_data_cid
final_data_frame1['y_'] = final_pred_1


# In[321]:


saved_final_data_frame1 = final_data_frame1
sorted_final_df = final_data_frame1.sort_values('y_',ascending =0)
final_larg= sorted_final_df.head(50)


# In[322]:


final_larg = final_larg[['CID','y_']]
top_cid = final_larg["CID"].to_list()


# In[323]:


print("Top 50 Drug Leads identified with PubChem CID's (from degree 1 based regression) : ")
itter_count = 0
for itter in top_cid:
    itter_count+=1
    
    print(str(itter_count)+" : "+str(itter))


# In[324]:


top_cid_mix_df['Degree 1'] = top_cid


# In[325]:


final_cid_1 = pd.DataFrame()
final_cid_1['CID'] = top_cid
final_cid_1.to_csv("final_cid_degree_1.csv")


#  

#  

#  

# DEGREE :2

# In[326]:


model_dict[2]

features = list(model_dict[2][1])
features

x_trained =train_file[features].values
y_trained = train_file['y'].values
x_trained.shape,y_trained.shape


# In[327]:


X_train, X_test, y_train, y_test = train_test_split(x_trained, y_trained, test_size=0.1, random_state=0)


# In[328]:


#poly regression - degree 2
scaler = preprocessing.StandardScaler()
degree=2
polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
polyreg_scaled.fit(X_train,y_train)


# In[329]:


y_pred = polyreg_scaled.predict(X_test)
r2_score(y_test, y_pred, multioutput='uniform_average')
final_data_frameF1 = final_data_frame[final_data_frame['MolecularWeight'] <= 500]
final_data_frameF2 = final_data_frameF1[final_data_frameF1['XLogP'] <=5.6]


# In[330]:


pred_data_cid = final_data_frameF2["CID"].to_list()


# In[331]:


#features.append('CID')
new_features = features.append("CID")
features


# In[332]:


final_data_frame1 = final_data_frameF2[features]
for x in final_data_frame1.isnull().any():
    if x == True:
        final_data_frame1 = final_data_frame1.fillna(method='ffill')
final_data_frame1.drop(columns='CID',inplace=True)
final_data_frame1.shape


# In[333]:


final_pred_1 = polyreg_scaled.predict(final_data_frame1)
final_pred_1 = list(final_pred_1)
final_data_frame1['CID'] = pred_data_cid
final_data_frame1['y_'] = final_pred_1


# In[334]:


saved_final_data_frame1 = final_data_frame1
sorted_final_df = final_data_frame1.sort_values('y_',ascending =0)
final_larg= sorted_final_df.head(50)


# In[335]:


final_larg = final_larg[['CID','y_']]
top_cid = final_larg["CID"].to_list()


# In[336]:


print("Top 50 Drug Leads identified with PubChem CID's (from degree 2 based regression) :")
itter_count = 0
for itter in top_cid:
    itter_count+=1
    
    print(str(itter_count)+" : "+str(itter))


# In[337]:


top_cid_mix_df['Degree 2'] = top_cid


# In[338]:


final_cid_2 = pd.DataFrame()
final_cid_2['CID'] = top_cid
final_cid_2.to_csv("final_cid_degree_2.csv")


#  

#  

#  

# DEGREE :3

# In[339]:


model_dict[3]
features = list(model_dict[3][1])
features
x_trained =train_file[features].values
y_trained = train_file['y'].values
x_trained.shape,y_trained.shape


# In[340]:


X_train, X_test, y_train, y_test = train_test_split(x_trained, y_trained, test_size=0.1, random_state=0)


# In[341]:


#poly regression - degree 3
scaler = preprocessing.StandardScaler()
degree=3
polyreg_scaled=make_pipeline(PolynomialFeatures(degree),scaler,LinearRegression())
polyreg_scaled.fit(X_train,y_train)


# In[342]:


y_pred = polyreg_scaled.predict(X_test)
r2_score(y_test, y_pred, multioutput='uniform_average')
final_data_frameF1 = final_data_frame[final_data_frame['MolecularWeight'] <= 500]
final_data_frameF2 = final_data_frameF1[final_data_frameF1['XLogP'] <=5.6]


# In[343]:


pred_data_cid = final_data_frameF2["CID"].to_list()


# In[344]:


#features.append('CID')
new_features = features.append("CID")
features


# In[345]:


final_data_frame1 = final_data_frameF2[features]
for x in final_data_frame1.isnull().any():
    if x == True:
        final_data_frame1 = final_data_frame1.fillna(method='ffill')
final_data_frame1.drop(columns='CID',inplace=True)
final_data_frame1.shape


# In[346]:


final_pred_1 = polyreg_scaled.predict(final_data_frame1)
final_pred_1 = list(final_pred_1)
final_data_frame1['CID'] = pred_data_cid
final_data_frame1['y_'] = final_pred_1


# In[347]:


saved_final_data_frame1 = final_data_frame1
sorted_final_df = final_data_frame1.sort_values('y_',ascending =0)
final_larg= sorted_final_df.head(50)


# In[348]:


final_larg = final_larg[['CID','y_']]
top_cid = final_larg["CID"].to_list()


# In[349]:


print("Top 50 Drug Leads identified with PubChem CID's (from degree 2 based regression) :")
itter_count = 0
for itter in top_cid:
    itter_count+=1
    
    print(str(itter_count)+" : "+str(itter))


# In[350]:


top_cid_mix_df['Degree 3'] = top_cid


# In[351]:


final_cid_2 = pd.DataFrame()
final_cid_2['CID'] = top_cid
final_cid_2.to_csv("final_cid_degree_2.csv")


# In[352]:


top_cid_mix_df


# In[353]:


top_cid_mix_df.to_csv("TOP_CID_123.csv")



top_cid = top_cid_mix_df['Degree 3']

top_cid


# Automated In Silico Modelling and Analysis

# In[228]:


import urllib.request
from os.path import expanduser
import os
home = expanduser("~")
assert (os.path.isdir(home+'/MGLTools-1.5.6')), "AutoDockTools not found!"
#set preparation pathways
prepare_protein_path = '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -A bonds_hydrogens -U nphs_lps_waters_nonstdres -r'
prepare_ligand_path = '~/MGLTools-1.5.6/bin/pythonsh ~/MGLTools-1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py -A bonds_hydrogens -U nphs_lps -l'



l=[]
for itter in top_cid:
    l.append(itter)
for i in l:


# programmatic download of structure of drug lead compounds


    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/SDF'.format(str(i))
    urllib.request.urlretrieve(url, '{}'.format(str(i)+'.sdf'))
    
    
#   programmatic ligand preparation for autodock 
    print("Preparing ligand {} for docking".format(i))
    os.system('obabel {} -O {} --gen3d --best'.format(str(i)+'.sdf', str(i)+'.mol2'))
    os.system(prepare_ligand_path + '{}'.format(str(i)+'.mol2'))
    
    

    
    print("Running docking procedure for Hepatitis C viral drug target, 1A1V ")
    os.system(
    'vina --config coordinatesandsetups1A1V.txt --receptor 1A1V.pdbqt --ligand {} --out {} --log {}'.format(
    str(i)+'.pdbqt', str(i)+'.1A1V.pdbqt', str(i)+'.1A1V.log')
)
    print("Docking to 1A1V protein completed. Files saved as {} as docked complex and {} as a logfile".format(
        str(i)+'.1A1V.pdbqt',str(i)+'.1A1V.log'))
    

            
    
    print(url)
    
    


# In[ ]:


"""Now postprocessing the results"""
import os, pandas as pd, requests
from biopandas.pdb import PandasPdb

for file in os.listdir():
    if file.endswith('.1A1V.pdbqt'): 
        os.system("obabel {} -O {} -l 1".format(file, file.strip('.1A1V.pdbqt')+'.best.pdb')) 
        prot_df = PandasPdb().read_pdb('1A1V.pdb') 
        flex_df = PandasPdb().read_pdb(file.strip('.1A1V.pdbqt')+'.best.pdb') 
        flex_df.df['ATOM']['chain_id'].replace(to_replace='A', value = "X", inplace = True)  
        prot_df.df['ATOM'] = prot_df.df['ATOM'].append(flex_df.df['ATOM']) 
        prot_df.to_pdb(path='1A1V-{}.cplx.pdb'.format(file.strip('.1A1V.pdbqt')),records=['ATOM'],gz=False,
                     append_newline=True) 
    
                
        


# In[ ]:


"""Uploads pdb and downloads protein-ligand interaction images"""
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.options import Options
import time


# In[ ]:


options=Options()
# options.headless = True
plip = webdriver.Firefox(options=options)
plip.get("https://projects.biotec.tu-dresden.de/plip-web/plip")
print("Connection successful")

for cplx in os.listdir():
    if cplx.endswith('.cplx.pdb'):
        print("Processing file {}".format(cplx))
    
        select_pdb_input = plip.find_element_by_xpath("//*[@id='select-pdb-by-file']").click() 

        browse = plip.find_element_by_xpath(
            '/html/body/div[1]/div[2]/div/form/div[1]/div[1]/div[3]/input'
        ).send_keys(                                                       
            os.getcwd()+'/{}'.format(cplx)
        )

        send_file = plip.find_element_by_xpath("//*[@id='submit']").click() 
        time.sleep(10) 
        try:
            try:
                open_interactions_1 = plip.find_element_by_xpath('/html/body/div/div[2]/div/div[1]/h2[2]').click()
                open_interactions_2 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/div[2]/h3').click()
                open_interactions_3 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/div[2]/div/h4').click()
                pngs = plip.find_elements_by_xpath("//a[contains(@href,'.png')]")
                pymolsessions = plip.find_elements_by_xpath("//a[contains(@href,'.pse')]")
                
            except:
                open_interactions_1 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/h2[1]').click()
                open_interactions_2 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/div[1]/h3').click()
                open_interactions_3 = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/div[1]/div[1]/div/h4').click()
                pngs = plip.find_elements_by_xpath("//a[contains(@href,'.png')]")
                pymolsessions = plip.find_elements_by_xpath("//a[contains(@href,'.pse')]")

            for image in pngs:
                print(image.get_attribute("href"))
                output_image = requests.get(image.get_attribute("href"))
                open(
                    os.getcwd()+'/{}'.format(cplx+'.png'), 'wb'
                ).write(output_image.content)
                print("Image saved as {}".format(cplx+'.png'))

            for pysession in pymolsessions:
                print(pysession.get_attribute("href"))
                pse = requests.get(pysession.get_attribute("href"))
                open(
                    os.getcwd()+'/{}'.format(cplx+'.pse'), 'wb'
                ).write(pse.content)
                print("Pymol sessions saved as {}".format(cplx+'.pse'))
                  
            restart_plip = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/p[3]/a').click()
            time.sleep(5)
        except:
            print("No interactions found for {} or damaged structure".format(cplx))
            try:
                restart_plip = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div/p[3]/a').click()
                time.sleep(5)
            except:
                restart_plip = plip.find_element_by_xpath('/html/body/div[1]/div[2]/div[2]/p/a').click() 


#  

# PROJECT ENDS
