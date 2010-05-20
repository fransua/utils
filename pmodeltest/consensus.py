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
    outgroup_name = Tree(trees[0]).get_leaf_names()[1]
    for (tree, weight) in zip (trees, weights):
        tree = Tree(tree)
        outgroup = tree.search_nodes(name=outgroup_name)[0]
        tree.set_outgroup(outgroup)
        dad = outgroup.get_sisters()[0]
        for node in dad.traverse():
            if node.is_root(): continue
            cluster   = ','.join (sorted (node.get_leaf_names()))
            dad_clust = ','.join (sorted (node.up.get_leaf_names()))
            if dic.has_key(dad_clust):
                if dic[dad_clust].has_key(cluster):
                    dic[dad_clust][cluster] += weight
                else:
                    dic[dad_clust][cluster] =  weight
            else:
                dic[dad_clust] = {cluster: weight}

    sorted_nodes = map(lambda x: x[1], sorted (\
        map (lambda x: (len (x.split(',')), x), \
             dic.keys()), reverse = True))

    # sorted( map (lambda x: [sorted(dic[x].values(),reverse=True),x],dic.keys()))

    cons_tree = Tree()
    cons_tree.add_child(name=outgroup_name)
    node = cons_tree.add_child(name='NoName')
    node.add_feature('childrens', set (tree.get_leaf_names())-set([outgroup_name]))

    while 'NoName' in cons_tree.get_leaf_names():
        for dad in dic.keys():
            if not dad in sorted_nodes: continue
            try:
                node = cons_tree.search_nodes(childrens=set(dad.split(',')))[0]
            except IndexError:
                continue
            better_child = [0, '']
            for child in dic[dad].keys():
                if dic[dad][child] > better_child[0]:
                    better_child = [dic[dad][child], child]
            weight, name = better_child
            brother = set(name.split(','))
            if len (brother) == 1:
                n = node.add_child(name=name)
                n.add_feature('childrens', brother)
            else:
                n = node.add_child(name='NoName')
                n.add_feature('childrens', brother)
            sister = set(dad.split(',')) - brother
            name = ','.join (sorted (list (sister)))
            weight = dic[dad][name]
            if len (sister) == 1:
                n = node.add_child(name=name)
                n.add_feature('childrens', sister)
            else:
                n = node.add_child(name='NoName')
                n.add_feature('childrens', sister)
            if not (   len (sister) == 1 and len(brother) > 1
                    or len (sister) > 1 and len(brother) == 1):
                n.up.support = weight
            sorted_nodes.pop(sorted_nodes.index(dad))
            



    sorted_nodes = map(lambda x: x[1], sorted (map (lambda x: (len (x.split(',')), x)\
                                                    , dic.keys()), reverse = True))
    cons_tree = Tree()
    node = cons_tree.add_child(name='NoName')
    node.add_feature('childrens',set (sorted_nodes[0].split(',')))
    sorted_nodes.pop(0)
    while sorted_nodes != []:
        for key in sorted_nodes:
            print key
            childrens = set (key.split(','))
            weight = dic[key]
            #if weight < float (max (dic.values())) / 2: continue
            for node in cons_tree.traverse(strategy='postorder'):
                if node.is_root(): continue
                if len (node.childrens & childrens) == len (childrens):
                    if not node.get_children() == []:
                        for child in node.get_children():
                            if len (child.childrens & childrens) > 0:
                                if child.support > weight:
                                    sorted_nodes.pop(sorted_nodes.index(key))
                                    break
                                else:
                                    print cons_tree
                                    child.detach()
                    if not key in sorted_nodes: break
                    if len (node.get_children()) == 2:
                        continue
                    if len (childrens) == 1:
                        new_node = node.add_child(name=key)
                        new_node.add_feature('childrens',childrens)
                        new_node.support = weight
                        sorted_nodes.pop(sorted_nodes.index(key))
                    else:
                        sorted_nodes.pop(sorted_nodes.index(key))
                        new_node = node.add_child(name='NoName')
                        new_node.add_feature('childrens',childrens)
                        new_node.support = weight
                        print cons_tree
                    break
            else:
                print 'hu?'
                new_node = cons_tree.add_child(name=key)
                new_node.add_feature('childrens',childrens)
                

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
