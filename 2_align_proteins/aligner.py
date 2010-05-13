#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/05/10 14:17:46

from optparse import OptionParser
from re import match, split
from sys import stderr
from subprocess import Popen, PIPE

__version__ = "0.2"
__title__   = "aligner v%s" % __version__

def main():
    '''
    main function when called by command line.
    '''
    opts = get_options()

    seqs     = {}
    for seq in read_fasta(opts.fastafile):
        seq['trseq'] = translate(seq['seq'], stop=opts.remove_stop)
        seqs[seq['name']] = seq

    prot_path   = opts.outfile + '_prot'
    aali_path   = opts.outfile + '_aa_ali'
    ali_path    = opts.outfile + '_ali'
    trimsq_path = opts.outfile + '_trimseq'
    trimcl_path = opts.outfile + '_trimcol'
    score_path  = opts.outfile + '_score'
    map_path    = opts.outfile + '_map'
    todel       = [prot_path]
    
    write_fasta(seqs, prot_path, clean=True, typ='trseq')

    if opts.only_translate:
        exit()

    ###########
    # ALIGN
    proc = Popen([opts.muscle_bin,
                  '-stable',
                  '-quiet',
                  '-noanchors',
                  '-maxiters' , '999',
                  '-maxhours' , '24 ',
                  '-maxtrees' , '100',
                  '-in'       , prot_path,
                  '-out'      , aali_path,
                  '-scorefile', score_path # must be last!!! because option...
                  ][:None if opts.score else -2], stdout=PIPE)
    if proc.communicate()[1] is not None:
        print >> stderr, proc.communicate()[0]
        exit('\nERROR: runninge muscle')

    ###########
    # TRIM SEQS
    if opts.trimseq:
        todel.append(trimsq_path)
        proc = Popen([opts.trimal_bin,
                      '-in'        , aali_path,
                      '-out'       , trimsq_path,
                      '-resoverlap', opts.resovlp,
                      '-seqoverlap', opts.seqovlp,
                      '-cons'      , '100'
                      ], stdout=PIPE)
        if proc.communicate()[1] is not None:
            print >> stderr, proc.communicate()[0]
            exit('\nERROR: runninge muscle')

        for seq in read_fasta(trimsq_path):
            seqs[seq['name']]['ali'] = seq['seq']
        
        trimmed = filter (lambda x: not seqs[x].has_key('ali'), seqs)
        if not opts.quiet:
            print >> stderr, 'WARNING: trimmed sequences: \n\t' + \
                  '\n\t'.join(trimmed)

        for s in seqs.keys():
            if s in trimmed:
                del(seqs[s])
        aali_path = trimsq_path
    else:
        for seq in read_fasta(aali_path):
            seqs[seq['name']]['ali'] = seq['seq']

    ###########
    # CODON MAP
    seqs = map2codons(seqs)

    ###########
    # TRIM COLS
    if opts.trimcol:
        todel.append(trimcl_path)
        proc = Popen([opts.trimal_bin,
                      '-in' , aali_path,
                      '-out', trimcl_path, 
                      '-automated1',
                      '-colnumbering'
                      ], stdout=PIPE)

        (keeplist, err) = proc.communicate()
        if err is not None:
            exit('ERROR: trimming columns.')

        keeplist = str (keeplist).strip().split(', ')

        #algt = get_alignment(seqs)
        for s in seqs:
            codons = divide (seqs[s]['codons'], rm_cod=False)
            for c in range (len (codons)):
                if not str(c) in keeplist and codons[c] != '---':
                    codons[c] = 'NNN'
            seqs[s]['codons'] = ''.join(codons)

    ###########
    # SEQ MAP
    if opts.printmap:
        _printmap(seqs, map_path, opts.pymap)
    write_fasta(seqs, ali_path, clean=opts.clean, typ='codons')

    Popen(['rm', '-f'] + todel, stdout=PIPE)


def _printmap(seqs, map_path, pymap=True):
    '''
    for a given alignment read from file, and the corresponding
    nucleotides sequences object (before alignment), create a map file
    '''
    out = open(map_path,'w')
    if not pymap:
        out.write('Aligned Sequences map:\n'+\
                  '\t\t(position in sequence/position in alignment)\n\n')
    for s in seqs.keys():
        codons = divide(seqs[s]['codons'], rm_cod=False)
        seqcount = map(lambda x: \
                       str (len (codons[:x])-codons[:x].count('---')) \
                       *(codons[x] != '---'), range(len(codons)))
        if pymap:
            out.write(seqs[s]['name']+'\t'+' '.join(seqcount)+'\n')
            continue
        gapmap = zip (seqcount, map(str, range(len (codons))))
        prestring  =  map('/'.join, gapmap)
        poststring = map(lambda x: ('%20s' % (prestring[x])*(x%5==0)), \
                         range (len (prestring)))
        out.write(' '*3 + ''.join(poststring)+'\n')
        out.write('  ' + ''.join(map(lambda x: ('%20s' % ('|')*(x%5==0)), \
                                      range (len (prestring))))+'\n')
        out.write('%-20s' % (seqs[s]['name'])+' '.join(codons)+'\n\n')
    out.close()

def map2codons(seqs):
    '''
    map amino-acid alignment to nt
    '''
    for s in seqs:
        codons     = divide(seqs[s]['seq'], rm_cod=False)
        ali_codons = ''
        for aa in seqs[s]['ali']:
            if aa == '-':
                ali_codons += '---'
            else:
                ali_codons += codons.pop(0)
        seqs[s]['codons'] = ali_codons
    return seqs

def read_fasta(infile):
    '''
    read file in fasta format and yield each sequence
    '''
    nam   = None
    descr = None
    seq   = ''
    for line in open(infile):
        line = line.strip()
        if line.startswith('>'):
            if nam is not None:
                yield { 'name'  : nam,
                        'descr' : descr,
                        'seq'   : seq
                    }
            items = split('[ |\t]', line, maxsplit=1)
            nam   = items[0].lstrip('>')
            descr = items[1] if len (items) == 2 else None
            seq = ''
            continue
        seq += line
    yield { 'name'  : nam,
            'descr' : descr,
            'seq'   : seq
            }

def write_fasta(seqs, outfile, clean=False, typ='seq'):
    '''
    just write fasta file from seq dic
    '''
    out = open(outfile, 'w')
    for s in seqs:
        seq = seqs[s]
        nam = seq['name'] if clean else seq['name'] + ' ' + seq['descr']
        out.write('>%s\n' % (nam))
        out.write('\n'.join(divide(seq[typ], 60, rm_cod=False))+'\n')
    out.close()

def divide(seq, size=3, rm_cod = True):
    '''
    divide a string into a list of string of given size.
    rest will be lose.
    '''
    codons = []
    for i in xrange(0, len(seq), size):
        if not match('^[aAtTgGcC]*$', str (seq[i:i+size])) and rm_cod:
            continue
        else :
            codons.append(seq[i:i+size])
    return codons

def translate(sequence, stop=False):
    '''
    little function to translate DNA to protein...
    from: http://python.genedrift.org/
    '''
    #dictionary with the genetic code
    gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'0', 'TAG':'0',
    'TGC':'C', 'TGT':'C', 'TGA':'0', 'TGG':'W',
    '---':'-'
    }
    ambig = {'Y':['A', 'G'], 'R':['C', 'T'], 'M':['G', 'T'], 'K':['A', 'C'], \
             'S':['G', 'C'],'W':['A', 'T'], 'V':['C', 'G', 'T'], \
             'H':['A', 'G', 'T'], 'D':['A', 'C', 'T'], 'B':['A', 'C', 'G'], \
             'N':['A', 'C', 'G', 'T']}
    proteinseq = ''
    #loop to read DNA sequence in codons, 3 nucleotides at a time
    for n in range(0, len(sequence), 3):
        #checking to see if the dictionary has the key
        try:
            proteinseq += gencode[sequence[n:n+3]]
        except KeyError:
            if len (sequence[n:n+3])< 3 :
                break
            newcod = []
            for nt in sequence[n:n+3]:
                if ambig.has_key(nt):
                    newcod.append(ambig[nt])
                else :
                    newcod.append(list (nt))
            aa = ''
            for nt1 in newcod[0]:
                for nt2 in newcod[1]:
                    for nt3 in newcod[2]:
                        try:
                            if aa == '':
                                aa  = gencode[nt1+nt2+nt3]
                            elif gencode[nt1+nt2+nt3] != aa:
                                aa = 'X'
                                break
                        except KeyError:
                            aa = 'X'
                            break
            proteinseq += aa
    #return protein sequence
    if stop:
        if proteinseq.endswith('0'):
            return proteinseq[:-1]
        else:
            return proteinseq
    else:
        return proteinseq

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
    parser.add_option('-i', dest='fastafile', metavar="PATH", \
                      help='path to input file in fasta format')
    parser.add_option('-o', dest='outfile', metavar="PATH", \
                      help='path to output file in fasta format')
    parser.add_option('-t', '--trimseqs', action='store_true', \
                      dest='trimseq', default=False, \
                      help='[%default] remove bad sequences (uses trimAl).')
    parser.add_option('-m', '--maskcols', action='store_true', \
                      dest='trimcol', default=False, \
                      help=\
                      '[%default] mask (with "N") bad columns (uses trimAl).')
    parser.add_option('-M', '--printmap', action='store_true', \
                      dest='printmap', default=False, \
                      help=\
                      '''[%default] save a map of alignement not human
                      friendly by default, see "--humanmap" option''')
    parser.add_option('--musclepath', dest='muscle_bin', \
                      metavar="PATH", help=\
                      '[%default] path to muscle binary.', \
                      default='/usr/bin/muscle')
    parser.add_option('--trimalpath', dest='trimal_bin', \
                      metavar="PATH", help=
                      '[%default] path to trimal binary.', \
                      default='/usr/local/bin/trimal')
    parser.add_option('--resoverlap', dest='resovlp', \
                      metavar="FLOAT", default='0.7', help=\
                      '''[%default] Minimum overlap of a positions with
                      other positions in the column to be considered a
                      "good position". (see trimAl User Guide).''')
    parser.add_option('--seqoverlap', dest='seqovlp', \
                      metavar="PERCENT", default='70', help=\
                      '''[%default] Minimum percentage of "good 
                      positions" that a sequence must have in order to
                      be conserved. (see trimAl User Guide).''')
    parser.add_option('--translate', action='store_true', \
                      dest='only_translate', default=False, \
                      help=\
                      '[%default] do not align just translate fasta.')
    parser.add_option('-q', '--quiet', action='store_false', \
                      dest='quiet', default=True, \
                      help='[%default] shut!')
    parser.add_option('-c', '--cleannames', action='store_false', \
                      dest='clean', default=True, \
                      help='[%default] removes sequence decription.')
    parser.add_option('--musclescore', action='store_true', \
                      dest='score', default=False, \
                      help='[%default] generate muscle score file.')
    parser.add_option('-r', '--remove_stop', action='store_false', \
                      dest='remove_stop', default=True, \
                      help=\
                      '[%default] remove stop codons from alignment.')
    parser.add_option('--humanmap', action='store_false', \
                      dest='pymap', default=True, \
                      help=\
                      '[False] print human readable map.')
    opts = parser.parse_args()[0]
    if not opts.outfile or not opts.fastafile:
        exit(parser.print_help())
    return opts


if __name__ == "__main__":
    exit(main())
