#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/05/14 08:58:48
#
# This script 

from optparse import OptionParser
from subprocess import Popen, PIPE
from re import match
from numpy import exp

__version__ = "0.0"
__title__   = "pmodeltest v%s" % __version__

models = [ ['-m', '000000'],
           ['-m', '010010'],
           ['-m', '010020'],
           ['-m', '012210'],
           ['-m', '010212'],
           ['-m', '012012'],
           ['-m', '012230'],
           ['-m', '010232'],
           ['-m', '012032'],
           ['-m', '012314'],
           ['-m', '012345'] ]

freqs = { 'ef': ['-f', '0.25 0.25 0.25 0.25'],
          'uf': ['-f', 'e'                  ] }

invts = { ''  : [],
          '+I': ['-v', 'e'  ] }

gamma = { ''  : ['-a', '1.0'],
          '+G': ['-c', '4', '-a', 'e'] }

modelnames = { '000000' + 'ef': ['JC'     , 0],
               '010010' + 'ef': ['K80'    , 1],
               '010020' + 'ef': ['TrNef'  , 2],
               '012210' + 'ef': ['TPM1'   , 2],
               '010212' + 'ef': ['TPM2'   , 2],
               '012012' + 'ef': ['TPM3'   , 2],
               '012230' + 'ef': ['TIM1ef' , 3],
               '010232' + 'ef': ['TIM2ef' , 3],
               '012032' + 'ef': ['TIM3ef' , 3],
               '012314' + 'ef': ['TVMef'  , 4],
               '012345' + 'ef': ['SYM'    , 5],
               '000000' + 'uf': ['F81'    , 3],
               '010010' + 'uf': ['HKY'    , 4],
               '010020' + 'uf': ['TrN'    , 5],
               '012210' + 'uf': ['TPM1uf' , 5],
               '010212' + 'uf': ['TPM2uf' , 5],
               '012012' + 'uf': ['TPM3uf' , 5],
               '012230' + 'uf': ['TIM1'   , 6],
               '010232' + 'uf': ['TIM2'   , 6],
               '012032' + 'uf': ['TIM3'   , 6],
               '012314' + 'uf': ['TVM'    , 7],
               '012345' + 'uf': ['GTR'    , 8]
               }

def main():
    '''
    main function when called by command line.
    infile must be in phylip format.
    '''
    opts = get_options()

    #ali_apth = opts.algt

    results = {}
    for model in models:
        for freq in freqs.keys():
            for inv in invts.keys():
                for gam in gamma.keys():
                    model_name  = modelnames[model[1] + freq][0]
                    model_param = modelnames[model[1] + freq][1]
                    print 'Model ' + \
                          model_name + inv + gam
                    #print 'Command line = ' +\
                    #     ' '.join(model + freqs[freq] + invts[inv] +gamma[gam])
                    (out, err) = Popen(['bin/phyml',
                                        '-i', opts.algt,
                                        '-d', 'nt',
                                        '-n', '1',
                                        '-b', '0'] + \
                                       model + freqs[freq] + \
                                       invts[inv] + gamma[gam],
                                       stdout=PIPE).communicate()
                    (numspe, lnl) = parse_stats(opts.algt + '_phyml_stats.txt')
                    tree          = get_tree   (opts.algt + '_phyml_tree.txt') 
                    numparam = model_param + \
                               (inv != '') + (gam != '') + numspe*2-3 + 1
                    aic = 2*numparam-2*lnl
                    print 'K = '+str (numparam)+', lnL = '+str(lnl) + \
                          '\nAIC = ' + str (aic)
                    print '-----------------------------------'
                    if err is not None:
                        exit ('problem running phyml: '+out)
                    results[model_name + inv + gam] =  { 'AIC' : aic,
                                                         'lnL' : lnl,
                                                         'K'   : numparam,
                                                         'tree': tree}

    # lala
    ord_aic = sorted (map (lambda x: [results[x]['AIC'], x], results.keys()))
    ord_aic = map (lambda x: x[1], ord_aic)
    min_aic = results[ord_aic[0]]['AIC']
    
    for model in ord_aic:
        results[model]['deltar'] =  results[model]['AIC'] - min_aic
        results[model]['weight'] = exp (-0.5 * results[model]['deltar'])

    sumweight = sum (map (lambda x: results[x]['weight'], results.keys()))

    cumweight = 0
    good_models = []
    for model in ord_aic:
        results[model]['weight'] = results[model]['weight']/sumweight
        cumweight += results[model]['weight']
        results[model]['cumweight'] = results[model]['weight'] + cumweight
        if results[model]['cumweight'] < 0.9999:
            good_models.append([model, int(1000*results[model]['weight']+0.5)])

    print '\n'.join (map (lambda x: '%-12s'%(x) + \
                          '%-10s' % (str (results[x]['AIC'])       ) + '\t' +\
                          '%-9s' % (str (results[x]['deltar'])    ) + '\t' +\
                          '%-17s' % (str (results[x]['weight'])    ) + '\t' +\
                          '%-17s' % (str (results[x]['cumweight']) ) \
                          , ord_aic))


    print '\n\n Models keeped to build consensus: \n' + \
          ', '.join (map(lambda x: x[0], good_models))

    tree_file = open(opts.path + 'intree', 'w')
    for model, weight in good_models:
        for i in range (weight):
            tree_file.write(results[model]['tree'])
    tree_file.close()

    Popen(['yes | bin/consense'], stdout=PIPE)

    final_tree   = get_tree(path + 'outtree')
    better_model = ord_aic[0]
    # FINI!!!! YUJUUUUUUUUUUUUUUUUU

# number of parameters = X (nb of branches) + 1 (topology) + Y (model)

def parse_stats(path):
    '''
    parse stats file of phyml, to extract the likelyhood value
    '''
    for line in open(path):
        if line.startswith('. Log-likelihood:'):
            lnl = float (line.strip().split()[-1])
        elif line.startswith('. Number of taxa:'):
            numspe = int (line.strip().split()[-1])
    return (numspe, lnl)

def get_tree(path):
    '''
    juste return the tree
    '''
    for line in open(path):
        if line.startswith('('):
            return line
    


#def get_alignment(seqs, typ='codons'):
#    '''
#    returns alignment from file
#    TODO: find better way than zip (*algt) to reverse it
#    '''
#    keyseqs = sorted(seqs.keys())
#    div = 3 if typ=='codons' else 1
#    seqlist   = map(lambda x: divide(seqs[x][typ], size=div), keyseqs)
#    align = [keyseqs] + zip ( *seqlist)
#    return align


def get_options():
    '''
    parse option from call
    '''
    parser = OptionParser(
        version=__title__,
        usage="%prog [options] file [options [file ...]]",
        description="""\
Reads sequeneces from file fasta format, and align acording to translation.
"""
        )
    models = '"HKY85,JC69,K80,F81,F84,TN93,GTR"'
    parser.add_option('-i', dest='algt', metavar="PATH", \
                      help='path to input file in fasta format')
    parser.add_option('-o', dest='outfile', metavar="PATH", \
                      help='path to output file in fasta format')
    parser.add_option('-t', '--trimseqs', action='store_true', \
                      dest='trimseq', default=False, \
                      help='[%default] remove bad sequences (uses trimAl).')
    parser.add_option('-m', metavar='LIST', \
                      dest='models', default=models, \
                      help=\
                      '[%default] DNA models.')
    opts = parser.parse_args()[0]
    if not opts.algt:
        exit(parser.print_help())
    return opts









if __name__ == "__main__":
    exit(main())
