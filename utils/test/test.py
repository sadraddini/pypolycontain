#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May  7 17:11:06 2019

@author: sadra
"""

import numpy as np

from pypolycontain.utils.utils import unique_rows

test={
      "utils.unique":True
      }

if test["utils.unique"]:
    E=np.random.randint(-2,2,(50,2))
    e=np.random.randint(-2,2,(50,1))
    E_1,e_1=unique_rows(E,e)
#    for row in range(E_1.shape[0]):
#        E_1[row,:]