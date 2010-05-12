#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/05/10 14:17:46

from optparse import OptionParser
from os import system
from re import match, split
from utils import translate
from sys import stderr

__version__ = "0.1_beta"
__title__   = "aligner v%s" % __version__

def main():

    '''
    main function when called by command line.
    '''
    parser = OptionParser(
        version=__title__,
        usage="%prog [options] file [options [file ...]]",
        description="""\
Reads sequeneces from file fasta format, align them.
"""
        )
    parser.add_option('-i', dest='fastafile', metavar="PATH", \
                      help='path to input file in fasta format')
    parser.add_option('-o', dest='outfile', metavar="PATH", \
                      help='path to output file in fasta format')
    parser.add_option('--musclepath', dest='muscle_bin', \
                      metavar="PATH", help=\
                      'path to muscle binary. [default: %default]', \
                      default='/usr/bin/muscle')
    parser.add_option('--trimalpath', dest='trimal_bin', \
                      metavar="PATH", help=\
                      'path to trimal binary. [default: %default]', \
                      default='/usr/local/bin/trimal')
    parser.add_option('--resoverlap', dest='resovlp', \
                      metavar="PATH", default='0.7', help=\
                      '''Minimum overlap of a positions with other positions
                      in the column to be considered a "good position".
                      (see trimAl User Guide). [default: %default]''')
    parser.add_option('--seqoverlap', dest='seqovlp', \
                      metavar="PATH", default='70', help=\
                      '''Minimum percentage of "good positions" that a
                      sequence must have in order to be conserved.
                      (see trimAl User Guide). [default: %default]''')
    parser.add_option('-t', '--translate', action='store_true', \
                      dest='only_translate', default=False, \
                      help='do not align just translate fasta.')
    parser.add_option('-r', '--remove_stop', action='store_false', \
                      dest='remove_stop', default=True, \
                      help='do not align just translate fasta.')
    opts = parser.parse_args()[0]
    if not opts.outfile or not opts.fastafile:
        exit(parser.print_help())

    seqs     = {}
    for seq in read_fasta(opts.fastafile):
        seq['trseq'] = translate(seq['seq'], stop=opts.remove_stop)
        seqs[seq['name']] = seq

    prot_path   = opts.outfile + '_prot'
    aali_path   = opts.outfile + '_aa_ali'
    ali_path    = opts.outfile + '_ali'
    trim_path   = opts.outfile + '_trim'
    score_path  = opts.outfile + '_score'

    write_fasta(seqs, prot_path, clean=True, typ='trseq')

    if opts.only_translate:
        exit()

    ###########
    # ALIGN
    system(opts.muscle_bin + ' -stable -quiet -maxiters 999 -maxhours 24 ' + \
              '-maxtrees 100 -noanchors -in %s -out %s -scorefile %s '\
              % (prot_path, aali_path, score_path)
              )

    ###########
    # TRIM
    system(opts.trimal_bin + \
              ' -in %s -out %s -resoverlap %s -seqoverlap %s -cons 100' \
              % (aali_path, trim_path, opts.resovlp, opts.seqovlp))

    for seq in read_fasta(trim_path):
        seqs[seq['name']]['ali'] = seq['seq']

    trimmed = filter (lambda x: not seqs[x].has_key('ali'), seqs)
    print >> stderr, 'WARNING: trimmed sequences: \n\t' + '\n\t'.join(trimmed)

    for s in seqs.keys():
        if s in trimmed:
            del(seqs[s])

    ###########
    # MAP
    seqs  = map2codons(seqs)

    write_fasta(seqs, ali_path, clean=True, typ='codons')

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

if __name__ == "__main__":
    exit(main())
