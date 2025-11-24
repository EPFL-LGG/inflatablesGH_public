#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys; sys.path.append('..'); sys.path.append('../../');  sys.path.append('../../gmsh')


# In[2]:


import inflation
import numpy as np


# In[3]:


bsis = inflation.BendingStiffnessIntegralSensitivity()

bsis.update(1, 2, 1, -np.ones(5) * 9, np.ones(5))
