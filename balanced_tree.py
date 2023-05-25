# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 10:03:02 2023

@author: ALDANADS
"""

class Node:
    def __init__(self, data):
        self.data = data
        self.left = None
        self.right = None

def build_tree(arr):
    if not arr:
        return None
    if len(arr) == 1:
        return Node(arr[0])
    mid = len(arr)//2
    root = Node(None)
    root.left = build_tree(arr[:mid])
    root.right = build_tree(arr[mid:])
    return root

def update_data(root):
    if root is None:
        return
    
    # We start at the leafs
    update_data(root.left)
    update_data(root.right)
    
    # Superior nodes are the sum of their childrens - Leaf are tuples, but the
    # nodes are floats
    # Tuple[0] = Transition rate
    # Tuple[1] = Type of event
    # Tuple[2] = Particle selected
    if root.left is not None and root.right is not None:
        if type(root.left.data) != tuple: 
            aux_l = root.left.data
        else:
            aux_l = root.left.data[0]
            
        if type(root.right.data) != tuple:
            aux_r = root.right.data
        else:
            aux_r = root.right.data[0]

        root.data = aux_l + aux_r

    elif root.left is not None:
        if type(root.left.data) != tuple: 
            root.data = root.left.data
        else: 
            root.data = root.left.data[0] 
            
    elif root.right is not None:
        if type(root.right.data) != tuple:
            root.data = root.right.data
        else: 
            root.data = root.right.data[0] 
    return root.data

def search_value(root,target):

    # Base case --> We go all the way down and the node is null
    # ATTENTION! -> Take care, some of them are tuples and other floats
    if type(root.data) != tuple:
        aux = root.data
    else:
        aux = root.data[0]
        
    # If it is None, is the leaf
    if (root.left is not None) and (type(root.left.data) != tuple):
        aux_l = root.left.data
    elif (root.left is not None):
        aux_l = root.left.data[0]
    
    if root is None:
        return False
    
    # We find the value among the leafs - Only interested in leafs, so we set the two extra conditions:
    # root.left and root.right is None.
    # The target should be smaller than the leaf we are checking to select that node
    elif (root.left is None) and (root.right is None) and target <= aux:
        return root.data
    else:
        #print(root.data)
        if target <= aux_l:
            return search_value(root.left, target)
        else:
            # Remember: nodes are the sum of their children. 
            # Because of that, if the target is greater than the left side, 
            # we operate target - root.left.data
            return search_value(root.right,target-aux_l)
        
    
# example usage:
# arr = [1, 2, 3, 4, 5, 6, 7]
# arr.sort()
# root = build_tree(arr)
# total = update_data(root)
# print(total)
