#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2010/05/10 14:17:46

from optparse import OptionParser
from re import match
from re import compile as compil
from sys import stderr
from subprocess import Popen, PIPE

__version__ = "0.3"
__title__   = "aligner v%s" % __version__

def main():
    '''
    main function when called by command line.
    '''
    opts = get_options()

    log = '\n\n'
    gencode = _set_code(opts.code)
    seqs     = {}
    for seq in read_fasta(opts.fastafile):
        seq['trseq'] = translate(seq['seq'], gencode, stop=opts.remove_stop)
        seqs[seq['name']] = seq

    log += '   ' + str (len (seqs)) + ' sequences\n\n'
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
    if not opts.input_ali:
        proc = Popen([opts.muscle_bin,
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
       
        log += '   Muscle command line: \n' + \
               ' '.join([opts.muscle_bin, '-quiet', '-noanchors', '-maxiters' , \
                         '999', '-maxhours', '24 ', '-maxtrees', '100', '-in', \
                         prot_path, '-out', aali_path, '-scorefile', \
                         score_path][:None if opts.score else -2]) + '\n\n'

    else:
        proc = Popen(['cp', prot_path, aali_path], stdout=PIPE)
        if proc.communicate()[1] is not None:
            print >> stderr, proc.communicate()[0]
            exit('\nERROR: when skipping muscle.')


    ###########
    # TRIM SEQS
    if opts.trimseq != False:
        todel.append(trimsq_path)
        proc = Popen([opts.trimal_bin,
                      '-in'        , aali_path,
                      '-out'       , trimsq_path,
                      '-resoverlap', str (opts.trimseq[1]),
                      '-seqoverlap', str (opts.trimseq[2]),
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
        if len (trimmed) > 0:
            log += '->trimmed sequences: \n\t' + \
                   '\n\t'.join(trimmed) + '\n'
        else: log += '->no trimmed sequences\n'

        for s in seqs.keys():
            if s in trimmed:
                del(seqs[s])
        aali_path = trimsq_path
        log += '   Trimal (sequences) command line: \n' + \
               ' '.join([opts.trimal_bin, '-in', aali_path, '-out',
                         trimsq_path, '-resoverlap', str (opts.trimseq[1]), \
                         '-seqoverlap', str (opts.trimseq[2]), '-cons', '100']) \
                         + '\n\n'
    else:
        for seq in read_fasta(aali_path):
            seqs[seq['name']]['ali'] = seq['seq']


    ###########
    # CODON MAP
    seqs = map2codons(seqs, opts.input_ali)

    ###########
    # TRIM COLS
    if opts.trimcol != 'None':
        if opts.trimcol == 'specific':
            todel.append(trimcl_path)
            proc = Popen([opts.trimal_bin,
                          '-in' , aali_path,
                          '-out', trimcl_path, 
                          '-gt' , str (opts.gaptreshold),
                          '-st' , str (opts.similarity),
                          '-colnumbering'
                          ], stdout=PIPE)
            (keeplist, err) = proc.communicate()
            if err is not None:
                exit('ERROR: trimming columns.')
            log += '   Trimal (columns) command line: \n' + \
                   ' '.join([opts.trimal_bin,
                          '-in', aali_path,
                          '-out', trimcl_path, 
                          '-gt', str (opts.gaptreshold),
                          '-st', str (opts.similarity),
                          '-colnumbering'
                          ]) + '\n'

        else:
            todel.append(trimcl_path)
            proc = Popen([opts.trimal_bin,
                          '-in' , aali_path,
                          '-out', trimcl_path, 
                          '-' + opts.trimcol,
                          '-colnumbering'
                          ], stdout=PIPE)
            (keeplist, err) = proc.communicate()
            if err is not None:
                exit('ERROR: trimming columns.')
            log += '   Trimal (columns) command line: \n' + \
                   ' '.join([opts.trimal_bin, '-in' , aali_path, '-out',
                             trimcl_path, '-' + opts.trimcol, \
                             '-colnumbering']) + '\n'

        keeplist = str (keeplist).strip().split(', ')

        algt = get_alignment(seqs)
        nnn = compil('[A-Z]{3}')
        if opts.nogap: 
            for (col, num) in zip (algt, range (len (algt))):
                if not str(num) in keeplist:
                    algt[num] = map (lambda x: nnn.sub('', x), col)
                    algt[num] = map (lambda x: compil('---').sub('', x), algt[num])
        else:
            for (col, num) in zip (algt, range (len (algt))):
                if not str(num) in keeplist:
                    algt[num] = map (lambda x: nnn.sub('NNN', x), col)
        for (key, seq) in zip (sorted (seqs.keys()), zip (*algt)):
            seqs[key]['codons'] = ''.join(seq)

    ###########
    # SEQ MAP
    if opts.printmap:
        _printmap(seqs, map_path, opts.pymap)
    write_fasta(seqs, ali_path, clean=opts.clean, typ='codons')

    Popen(['rm', '-f'] + todel, stdout=PIPE)

    if opts.print_log:
        print log

def _printmap(seqs, map_path, pymap=True):
    '''
    for a given alignment read from file, and the corresponding
    nucleotides sequences object (before alignment), create a map file
    '''
    out = open(map_path,'w')
    if not pymap:
        out.write('Aligned Sequences map:\n'+\
                  '\t\t(position in sequence/position in alignment)\n\n')

    cols = zip (*map(lambda x: divide(seqs[x]['codons'], rm_cod=False), seqs))
    nogapcount = []
    for i in cols:
        if '---' not in i:
            nogapcount.append (max (nogapcount)+1)
        else:
            nogapcount.append (0)
    nogapcount = map (lambda x: str(x-1), nogapcount)
    nogapcount = map (lambda x: ('-1' != x)*x+('-1'==x)*'*', nogapcount)
    if pymap:
        out.write ('no gap\t' + ' '.join (nogapcount) + '\n')
    for s in seqs.keys():
        codons = divide(seqs[s]['codons'], rm_cod=False)
        seqcount = map(lambda x: \
                       str (len (codons[:x])-codons[:x].count('---')) \
                       *(codons[x] != '---')+'*'*(codons[x] == '---'), \
                       range(len(codons)))
        if pymap:
            out.write(seqs[s]['name']+'\t'+' '.join(seqcount)+'\n')
            continue
        gapmap = zip (nogapcount, seqcount, map(str, range(len (codons))))
        prestring  =  map('/'.join, gapmap)
        poststring = map(lambda x: ('%-20s' % (prestring[x])*(x%5==0)), \
                         range (len (prestring)))
        out.write('%-20s' % (seqs[s]['name'])+'\n')
        out.write(' '*21 + ''.join(poststring)+'\n')
        out.write('  ' + ''.join(map(lambda x: ('%20s' % ('|')*(x%5==0)), \
                                      range (len (prestring))))+'\n')
        out.write(' '*20+' '.join(codons)+'\n')
    out.close()


def get_alignment(seqs, typ='codons'):
    '''
    returns alignment from file
    TODO: find better way than zip (*algt) to reverse it
    '''
    keyseqs   = sorted(seqs.keys())
    seqlist   = map(lambda x: divide(seqs[x][typ], rm_cod=False), keyseqs)
    align = zip( *seqlist)
    return align

def map2codons(seqs, input_ali):
    '''
    map amino-acid alignment to nt
    '''
    for s in seqs:
        if input_ali:
            seqs[s]['codons'] = seqs[s]['seq']
            continue
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
    blank_re = compil('[ \t]')
    for line in open(infile):
        line = line.strip()
        if line.startswith('>'):
            if nam is not None:
                if seq == '':
                    print >> stderr, 'ERROR: no sequence for ', str(nam)
                    exit()
                yield { 'name'  : nam,
                        'descr' : descr,
                        'seq'   : seq
                    }
            items = blank_re.split(line, maxsplit=1)
            nam   = items[0].lstrip('>')
            descr = items[1] if len (items) == 2 else None
            seq = ''
            continue
        seq += blank_re.sub('', line)
    if seq == '' and nam is not None:
        print >> stderr, 'ERROR: no sequence for ', str(nam)
        exit()
    elif seq == '':
        print >> stderr, 'ERROR: presence of repeated names'
        exit()
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

def _set_code(code):
    '''
    gencode, choose between:
        * [std] Standard
        * [vmt] Vertebrate Mitochondrial
        * [ymt] Yeast Mitochondrial
        * [mmt] Mold Mitochondrial, Protozoan Mitochondrial, Coelenterate Mitochondrial, Mycoplasma and Spiroplasma
        * [imt] Invertebrate Mitochondrial
        * [cnc] Ciliate Nuclear, Dasycladacean Nuclear, Hexamita Nuclear
        * [emi] Echinoderm Mitochondrial and Flatworm Mitochondrial
        * [enu] Euplotid Nuclear
        * [bpp] Bacterial and Plant Plastid
        * [ayn] Alternative Yeast Nuclear
        * [ami] Ascidian Mitochondrial
        * [afm] Alternative Flatworm Mitochondrial
        * [bma] Blepharisma Macronuclear
        * [cmi] Chlorophycean Mitochondrial
        * [tmi] Trematode Mitochondrial
        * [som] Scenedesmus obliquus Mitochondrial
        * [thm] Thraustochytrium Mitochondrial
    '''
    gencode = {
        'std' :{
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T',
            'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'CTA':'L', 'CTC':'L',
            'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R',
            'CGG':'R', 'CGT':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D',
            'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F',
            'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', '---':'-'},
        'vmt' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
            'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D', 'GAC':'D',
            'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
            'TAA':'*', 'TAG':'*', 'AGA':'*', 'AGG':'*', '---':'-'},
        'ymt' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'T', 'CTC':'T', 'CTA':'T', 'CTG':'T',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'mmt' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'imt' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'S', 'AGG':'S', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'cnc' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAA':'Q', 'TAG':'Q',
            'TGT':'C', 'TGC':'C', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L',
            'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H',
            'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R',
            'CGG':'R', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K',
            'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G',
            'GGC':'G', 'GGA':'G', 'GGG':'G', 'TGA':'*', '---':'-'},
        'emi' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'N', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'S', 'AGG':'S', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'enu' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'C', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'bpp' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P',
            'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H', 'CAA':'Q',
            'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'ATT':'I',
            'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T',
            'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S',
            'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V', 'GTA':'V',
            'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D',
            'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G',
            'GGG':'G', 'TAA':'*', 'TAG':'*', 'TGA':'*', '---':'-'},
        'ayn' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'S', 'CCT':'P',
            'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H', 'CAA':'Q',
            'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'ATT':'I',
            'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T',
            'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S',
            'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V', 'GTA':'V',
            'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D',
            'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G',
            'GGG':'G', 'TAA':'*', 'TAG':'*', 'TGA':'*', '---':'-'},
        'ami' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'G', 'AGG':'G', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'afm' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAA':'Y', 'TGT':'C',
            'TGC':'C', 'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L',
            'CTG':'L', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H',
            'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R',
            'CGG':'R', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T',
            'ACC':'T', 'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'N',
            'AAG':'K', 'AGT':'S', 'AGC':'S', 'AGA':'S', 'AGG':'S', 'GTT':'V',
            'GTC':'V', 'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A',
            'GCG':'A', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G',
            'GGC':'G', 'GGA':'G', 'GGG':'G', 'TAG':'*', '---':'-'},
        'bma' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAG':'Q', 'TGT':'C',
            'TGC':'C', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'cmi' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAG':'L', 'TGT':'C',
            'TGC':'C', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'tmi' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCA':'S', 'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C',
            'TGA':'W', 'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
            'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H',
            'CAA':'Q', 'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
            'ATT':'I', 'ATC':'I', 'ATA':'M', 'ATG':'M', 'ACT':'T', 'ACC':'T',
            'ACA':'T', 'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'N', 'AAG':'K',
            'AGT':'S', 'AGC':'S', 'AGA':'S', 'AGG':'S', 'GTT':'V', 'GTC':'V',
            'GTA':'V', 'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
            'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G',
            'GGA':'G', 'GGG':'G', 'TAA':'*', 'TAG':'*', '---':'-'},
        'som' :{
            'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'TCT':'S', 'TCC':'S',
            'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TAG':'L', 'TGT':'C', 'TGC':'C',
            'TGG':'W', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P',
            'CCC':'P', 'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H', 'CAA':'Q',
            'CAG':'Q', 'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'ATT':'I',
            'ATC':'I', 'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T',
            'ACG':'T', 'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S',
            'AGC':'S', 'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V', 'GTA':'V',
            'GTG':'V', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D',
            'GAC':'D', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G',
            'GGG':'G', 'TCA':'*', 'TAA':'*', 'TGA':'*', '---':'-'},
        'thm' :{
            'TTT':'F', 'TTC':'F', 'TTG':'L', 'TCT':'S', 'TCC':'S', 'TCA':'S',
            'TCG':'S', 'TAT':'Y', 'TAC':'Y', 'TGT':'C', 'TGC':'C', 'TGG':'W',
            'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'CCT':'P', 'CCC':'P',
            'CCA':'P', 'CCG':'P', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
            'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'ATT':'I', 'ATC':'I',
            'ATA':'I', 'ATG':'M', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
            'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K', 'AGT':'S', 'AGC':'S',
            'AGA':'R', 'AGG':'R', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
            'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'GAT':'D', 'GAC':'D',
            'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
            'TTA':'*', 'TAA':'*', 'TAG':'*', 'TGA':'*', '---':'-'},
    }
    return gencode [code]
    

def translate(sequence, gencode, stop=False):
    '''
    little function to translate DNA to protein...
    from: http://python.genedrift.org/ completed by biopython
    TODO: do not look at it too much :S
    '''
    #dictionary with the genetic code
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
        if proteinseq.endswith('*'):
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
    parser.add_option('-t', '--trimseqs', action='store_const', const=[True, 0.7, 70], \
                      dest='trimseq', default=False,\
                      help=
                      '''[%default] remove bad sequences (uses trimAl).
                      By default residue overlap is set to 0.7, and
                      sequence overlap to 70. Use --resoverlap and
                      --seqoverlap to change it.
                      ''')
    parser.add_option('--resoverlap', dest='trimseq[1]', action="store", \
                      metavar="FLOAT", default=0.7, type='float', \
                      help=\
                      '''[%default] Minimum overlap of a positions with
                      other positions in the column to be considered a
                      "good position". (see trimAl User Guide).''')
    parser.add_option('--seqoverlap', dest='trimseq[2]', action="store", \
                      metavar="PERCENT", default=70, type='int', help=\
                      '''[%default] Minimum percentage of "good 
                      positions" that a sequence must have in order to
                      be conserved. (see trimAl User Guide).''')
    parser.add_option('--nogap', action='store_true', \
                      dest='nogap', default=False, \
                      help=\
                      '''[%default] removes all gaps from alignement.
                      (uses trimAl).''')
    parser.add_option('-M', '--printmap', action='store_true', \
                      dest='printmap', default=False, \
                      help=\
                      '''[%default] save a map of alignement not human
                      friendly by default, see "--humanmap" option''')
    parser.add_option('--humanmap', action='store_false', \
                      dest='pymap', default=True, \
                      help=\
                      '[False] print human readable map. Only with -M option.')
    parser.add_option('--maskcol', metavar='OPTION', dest='trimcol', \
                      default='None', \
                      choices = ['None', 'automated1', 'softmasking', \
                                 'gapyout', 'strict', 'strictplus', \
                                 'specific'], \
                      help=\
                      '''[%default] mask (with "N") bad columns (uses
                      trimAl). Masking options are: None, automated1,
                      softmasking, gapyout, strict, strictplus or specific.
                      ''')
    parser.add_option('--gt', dest='gaptreshold', action="store", \
                      metavar="FLOAT", default=0, type='float', \
                      help=\
                      '''[%default] 1 - (fraction of sequences with a gap
                      allowed). Only use with specific maskcol option.
                      (see trimAl User Guide).''')
    parser.add_option('--st', dest='similarity', action="store", \
                      metavar="FLOAT", default=0, type='float', \
                      help=\
                      '''[%default] Minimum average similarity allowed
                      (see trimAl User Guide).''')
    parser.add_option('--alignment', dest='input_ali', \
                      action="store_true", help=\
                      '''[%default] Infile is already aligned. This option
                      will unactive muscle alignment process.''', \
                      default=False)
    parser.add_option('--musclepath', dest='muscle_bin', \
                      metavar="PATH", help=\
                      '[%default] path to muscle binary.', \
                      default='/usr/bin/muscle')
    parser.add_option('--trimalpath', dest='trimal_bin', \
                      metavar="PATH", help=
                      '[%default] path to trimal binary.', \
                      default='/usr/local/bin/trimal')
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
    parser.add_option('--log', action='store_true', \
                      dest='print_log', default=False, \
                      help=\
                      '[%default] Print aligner Log.')
    parser.add_option('--gencode', metavar='OPTION', dest='code', \
                      default='std', \
                      choices = ['std', 'vmt', 'ymt', 'mmt', 'imt', 'cnc',\
                                 'emi', 'enu', 'bpp', 'ayn', 'ami', 'afm',\
                                 'bma', 'cmi', 'tmi', 'som', 'thm'], \
                      help=\
                      '''[%default] Choose genetic code between:
                        std -> Standard                           
                        cnc -> Ciliate Nuclear, Dasycladacean Nuclear,
                        ---    Hexamita Nuclear                           
                        bpp -> Bacterial and Plant Plastid                           
                        ayn -> Alternative Yeast Nuclear                           
                        vmt -> Vertebrate Mitochondrial                           
                        ymt -> Yeast Mitochondrial                           
                        mmt -> Mold Mitochondrial, Protozoan  
                        ---    Mitochondrial, Coelenterate Mitochondrial,
                        ---    Mycoplasma and Spiroplasma                           
                        imt -> Invertebrate Mitochondrial                           
                        emi -> Echinoderm Mitochondrial and Flatworm
                        ---    Mitochondrial                           
                        enu -> Euplotid Nuclear                           
                        ami -> Ascidian Mitochondrial                           
                        afm -> Alternative Flatworm Mitochondrial                           
                        bma -> Blepharisma Macronuclear                           
                        cmi -> Chlorophycean Mitochondrial                           
                        tmi -> Trematode Mitochondrial                           
                        som -> Scenedesmus obliquus Mitochondrial                           
                        thm -> Thraustochytrium Mitochondrial                           
                      ''')
    opts = parser.parse_args()[0]
    if not opts.outfile or not opts.fastafile:
        exit(parser.print_help())
    return opts


if __name__ == "__main__":
    exit(main())
