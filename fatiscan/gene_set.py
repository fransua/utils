#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/08/17 16:46:02

# easy_install fisher
from fisher import pvalue
# in my extra_stats package
from extra_stats.fdr import bh_qvalues
from bisect import bisect_left

class Gene_set:
    '''
    Fatiscan with upper case, it is an object.
    '''
    def __init__(self, infile, annot, partitions=30, use_order=True):
        '''
        init function, what is done when object is called.
        '''

        # get gene list and corresponding values
        self.infile = infile
        self.use_order = use_order
        self.genes, self.values, self.order = self._parse_infile()

        # use oder of list instead of values.. or not
        # mak a dictionnary of genes to get values quickly
        self.values = dict (zip (self.genes, self.values))

        self.annot = self._parse_annot(annot)
        # sort genes in annot by their values...
        # useful to know which gene in list1 or list2
        self.annot = self._order_genes_in_annot()
        
    def _parse_infile(self):
        '''
        parse in file in format:
        geneID_1 <tab> value1
        geneID_2 <tab> value2
        ...
        genes should be ordered by value
        returns genes, values and order of values
        '''
        genes, values = zip (*sorted ([i.strip().split('\t') \
                              for i in open(self.infile)], \
                             key=lambda x: float(x[1])))
        values = map (float, values)
        if self.use_order:
            order = map (values.index, values)
        else:
            order = values[:]
        return genes, values, order

    def _parse_annot(self, annot):
        '''
        parse annotation file in format:
        annotationA <tab> geneID_1
        annotationA <tab> geneID_2
        ...
        '''
        dico = {}
        for gene, annot in [i.strip().split('\t') for i in open (annot)]:
            # only store genes that we have in our list
            if not self.values.has_key(gene):
                continue
            if dico.has_key(annot):
                dico[annot].append (gene)
            else:
                dico[annot] = [gene]
        return dico

    def _order_genes_in_annot(self):
        '''
        order genes in annot dict by their values
        '''
        dico = {}
        for annot in self.annot.iterkeys():
            dico[annot] = set (sorted (self.annot[annot], \
                                       key=lambda x: self.values[x]))
        return dico

    def run_gsea (self, partitions=30):
        '''
        run gsea needs python fisher, and fdr from extra stats
        '''
        rank      = float (max (self.order))/partitions
        pvalues   = []
        keys      = []
        total_len = len (self.genes)
        # intialize dict
        dico      = {}
        for annot in self.annot.keys():
            dico[annot] = [{}] *partitions
        # define cutoff value for each partition
        dico['thresh'] = [bisect_left (self.order, rank * (p + 1)) \
                         for p in xrange(partitions)]
        # start fishers
        for part in xrange(partitions):
            genes1 = set (self.genes[:self.order.index (dico['thresh'][part])])
            len_genes1 = len (genes1)
            len_genes2 = total_len - len_genes1
            for annot, annot_genes in self.annot.iteritems():
                p1 = len (annot_genes & genes1)
                p2 = len (annot_genes) - p1
                n1 = len_genes1  - p1
                n2 = len_genes2  - p2
                pv = pvalue(p1, n1, p2, n2).two_tail
                dico[annot][part]['p1' ] = p1
                dico[annot][part]['n1' ] = n1
                dico[annot][part]['p2' ] = p2
                dico[annot][part]['n2' ] = n2
                dico[annot][part]['pv' ] = pv
                pvalues.append(pv)
                keys.append((annot, part))
        # compute adjustment of pvalues
        adj_pvalues = bh_qvalues(pvalues)
        for annot, part in keys:
            dico[annot][part]['apv'] = adj_pvalues.pop(0)
        return dico

    def write_gsea (self):
        '''
        write to file, or pickle
        '''
        pass

