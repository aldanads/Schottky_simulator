# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 16:31:14 2022

@author: ALDANADS
"""

import shelve

filename = 'variables'

my_shelf = shelve.open(filename)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()