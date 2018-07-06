
# coding: utf-8

# In[6]:

import pandas as pd
df = pd.read_csv("labels.txt", sep='\t', index_col=0, header=0)
print(df.shape)


# In[7]:

print(df.head(5))


# In[8]:

df.sum(axis=0)


# In[9]:

melted_data = pd.melt(df, value_vars=['ENCSR000AKB', 'ENCSR000BJE','ENCSR000DZL'], 
                      var_name='Task', value_name='count')
print(melted_data.groupby(by=['Task', 'count'])['count'].count())


# In[ ]:




