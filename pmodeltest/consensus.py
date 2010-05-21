#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/05/16 14:48:48


import os,re
from ete2 import Tree


def consensus(trees,weights=[],lim=0):
    '''
    returns weighted consensus tree
    50% majority rule
    TODO: fix glup
    '''
    if weights == []: weights = [1] * len (trees)
    dic = {}
    outgroup_name = Tree(trees[0]).get_leaf_names()[1]
    for (tree, weight) in zip (trees, weights):
        tree = Tree(tree)
        outgroup = tree.search_nodes(name=outgroup_name)[0]
        tree.set_outgroup(outgroup)
        dad = outgroup.get_sisters()[0]
        for node in dad.traverse():
            if node.is_root(): continue
            cluster  = ','.join (sorted (node.get_leaf_names()))
            if dic.has_key(cluster):
                dic[cluster] += weight
            else:
                dic[cluster] =  weight

    sorted_nodes = map(lambda x: [x[2], x[1]], sorted (\
        map (lambda x: (len (x.split(',')), x, dic[x]), \
             dic.keys()), reverse = True))
    if lim == 0:
        lim = sorted (sorted_nodes, reverse=True)[:len (tree)*2 - 3][-1][0]
    sorted_nodes = filter (lambda x: x[0] >= lim, sorted_nodes)
    sorted_nodes = map (lambda x: x[1], sorted_nodes)
    if len (sorted_nodes) > len (tree)*2 - 3:
        print 'GLUP!!'
    cons_tree = Tree()
    cons_tree.add_child(name=outgroup_name)
    node = cons_tree.add_child(name='NoName')
    node.add_feature('childrens', set (tree.get_leaf_names())-set([outgroup_name]))
    sorted_nodes.pop(0)
    while len (sorted_nodes) > 0:
        for name in sorted_nodes:
            if not name in sorted_nodes: continue
            for node in cons_tree.traverse(strategy='postorder'):
                if node.is_root(): continue
                if node.name is not 'NoName': continue
                if len (node.childrens & set(name.split(','))) == 0:
                    continue
                weight = dic[name]
                children = set(name.split(','))
                if len (children) == 1:
                    node.add_child(name=name)
                else:
                    n = node.add_child(name='NoName')
                    n.add_feature('childrens', children)
                    n.support = weight
                break
            sorted_nodes.pop(sorted_nodes.index(name))
            sister = node.childrens - children
            name = ','.join (sorted ( list (sister)))
            weight = dic[name]
            if len (sister) == 1:
                node.add_child(name=name)
            else:
                n = node.add_child(name='NoName')
                n.add_feature('childrens', sister)
                n.support = weight
            sorted_nodes.pop(sorted_nodes.index(name))
            break                
    return cons_tree.write(format=9)
    


def main():
    '''
    '''
    trees = []
    tree_file = open('/home/francisco/toolbox/utils_dev/examples/100_tree.cons')
    for line in tree_file:
        trees.append(line)

    weights = []







if __name__ == "__main__":
    exit(main())
