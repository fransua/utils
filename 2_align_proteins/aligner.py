#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/05/10 14:17:46

from optparse import OptionParser
from os import system
from re import match, split
from translate import translate
from sys import stderr

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
    trim_path   = opts.outfile + '_trim'
    score_path  = opts.outfile + '_score'
    map_path    = opts.outfile + '_map'

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
    if opts.trim:
        system(opts.trimal_bin + \
               ' -in %s -out %s -resoverlap %s -seqoverlap %s -cons 100' \
               % (aali_path, trim_path, opts.resovlp, opts.seqovlp))

        for seq in read_fasta(trim_path):
            seqs[seq['name']]['ali'] = seq['seq']
        
        trimmed = filter (lambda x: not seqs[x].has_key('ali'), seqs)
        if not opts.quiet:
            print >> stderr, 'WARNING: trimmed sequences: \n\t' + \
                  '\n\t'.join(trimmed)

        for s in seqs.keys():
            if s in trimmed:
                del(seqs[s])
    else:
        for seq in read_fasta(aali_path):
            seqs[seq['name']]['ali'] = seq['seq']

    ###########
    # MAP
    seqs  = map2codons(seqs)

    if opts.printmap:
        _printmap(seqs, map_path, opts.pymap)

    write_fasta(seqs, ali_path, clean=True, typ='codons')

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

def get_alignment(seqs):
    '''
    returns alignment from file
    '''
    keyseqs   = sorted(seqs.keys())
    seqlist   = map(lambda x: divide(seqs[x]['codons'], rm_cod=False), keyseqs)
    align = zip( *seqlist )
    return align

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
    parser.add_option('-t', '--trim', action='store_true', \
                      dest='trim', default=False, \
                      help='[%default] trim alignment with trimAl.')
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
