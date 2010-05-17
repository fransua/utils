#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/05/16 14:48:48


import os,re
from ete2 import Tree






def consensus(trees,weights=[]):
    '''
    returns weighted consensus tree
    50% majority rule
    '''
    if weights == []: weights = [1] * len (trees)
    dic = {}
    for (tree, weight) in zip (trees, weights):
        tree = Tree(tree)
        for node in tree.traverse():
            if node.is_root(): continue
            cluster = ','.join(sorted(node.get_leaf_names()))
            if dic.has_key(cluster):
                dic[cluster] += weight
            else:
                dic[cluster] =  weight

    sorted_nodes = map(lambda x: x[1], sorted (map (lambda x: (len (x.split(',')), x)\
                                                    , dic.keys()), reverse = True))
    cons_tree = Tree()
    node = cons_tree.add_child(name='NoName')
    node.add_feature('childrens',set (sorted_nodes[0].split(',')))
    for key in sorted_nodes[1:]:
        childrens = set (key.split(','))
        weight = dic[key]
        if weight < float (max (dic.values())) / 2: continue
        for node in cons_tree.traverse(strategy='postorder'):
            if node.is_root(): continue
            if len (node.childrens & childrens) is not 0:
                if len (node.get_children()) is 2:
                    continue
                if len (childrens) is 1:
                    new_node = node.add_child(name=key)
                    new_node.add_feature('childrens',childrens)
                else:
                    new_node = node.add_child(name='NoName')
                    new_node.add_feature('childrens',childrens)
                    new_node.support = weight
                break
        else:
            new_node = cons_tree.add_child(name=key)
            new_node.add_feature('childrens',childrens)

    return cons_tree.write(format=9)
    


def main():
    '''
    '''
    trees = []
    tree_file = open('/home/francisco/toolbox/utils/pmodeltest/bin/intree')
    for line in tree_file:
        trees.append(line)








if __name__ == "__main__":
    exit(main())
