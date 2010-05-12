#!/usr/bin/python
#        Author: Francois-Jose Serra
# Creation Date: 2009/11/09 11:18:31

import sys, os, re
import getopt
from commands import getstatusoutput
from Bio import AlignIO, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def run_muscle(prot_records, nuclRecords, inProts, \
               outProts, trim = True, rm_bad = False):
    '''
    simple function that runs Muscle.
    '''
    bads = []
    protFile = open(inProts,'w')
    SeqIO.write(protRecords,protFile,'fasta')
    protFile.close()
    muscleDi = getstatusoutput('muscle -stable -quiet -maxiters 999 '\
                               +'-maxhours 24 -maxtrees 100 -noanchors -in '\
                               + inProts+' -out '+outProts+' -scorefile '+\
                               re.sub('_pep','_score',outProts))[1]
    if 'WARNING' in muscleDi:
        sys.exit("ERROR :P in: "+inProts+'\n')
    if trim:
        keepcols = getstatusoutput('/home/francisco/tools/trimAl/source/trimal -in '\
                                   +outProts+' -colnumbering -out '\
                                   +outProts+'_trim -automated1')[1].split(', ')
        keepcols = map (lambda x: int(re.sub('\n* *','',x)), keepcols)
        bads = list (set (range (0, max (keepcols))) - set(keepcols))
    if rm_bad:
       os.system('/home/francisco/tools/trimAl/source/trimal -in '+outProts+\
                  ' -out '+re.sub('_pep','_net_pep',outProts)+' -resoverlap 0.8 -seqoverlap 90 -cons 100')
       os.system('perl -e "s/  *[0-9][0-9]*  *bp\n/\n/" -p -i ' + re.sub('_pep','_net_pep',outProts))
       muscleDi = getstatusoutput('muscle -stable -quiet -maxiters 999 '\
                                  +'-maxhours 24 -maxtrees 100 -noanchors -in '\
                                  + re.sub('_pep','_net_pep',outProts)+\
                                  ' -out '+outProts+' -scorefile '+\
                                  re.sub('_pep','_score',outProts))[1]
    return bads

def _try_translate(seq, rm_cod = True):
    '''
    translate function that remove gaps and unknown characters, takes
    seq string in argument and returns the same string and its
    translation.
    '''
    seq = Seq(seq)
    errors = []
    prot = ''
    codons = divide(seq, rm_cod=rm_cod)
    to_del = []
    for i in range(0, len(codons)):
        try:
            prot = prot + codons[i].translate()
            if prot[-1]=='*':
                prot = Seq(prot.tostring()[:-1])
                errors.append(codons[i].tostring()+'(STOP)')
                to_del.append(i)
        except:
            to_del.append(i)
            errors.append(codons[i].tostring())
    if len(errors) > 0:
        to_del.reverse()
        map(codons.pop, to_del)
    codons = map(lambda x: x.tostring(),codons)
    return ''.join(codons), prot, codons

def divide(seq, size=3, rm_cod = True):
    '''
    divide a string into a list of string of given size.
    rest will be lose.
    '''
    codons = []
    for i in xrange(0,len(seq), size):
        if not re.match('^[aAtTgGcC]*$',str (seq[i:i+size])) and rm_cod:
            continue
        else :
            codons.append(seq[i:i+size])
    return codons

def alignment2map(nuclRecords, outProts,bads):
    try:
        h = open(outProts,'r')
        alignment = AlignIO.read(h,'fasta')
        h.close()
    except:
        sys.exit(" ERROR:P pep file not created for: "+\
                 re.sub('.*/dir_','',outProts)+'\n')
    protRecords=alignment.get_all_seqs()
    alignment.map2ungap=[]
    for x in protRecords:
        x.map=[]
        try:
            x.codons = divide(
                nuclRecords[\
            map(lambda y: y.id, nuclRecords).index(x.id)\
            ].seq.tostring())
        except ValueError:
            sys.exit('nt and prot not same length for: '+\
                     re.sub('.*/dir_','',outProts)+'\n')
        x.codons.insert(0,'   ')
    for col in range (0, alignment.get_alignment_length()):
        for aa in range(0, len(alignment.get_column(col))):
            if alignment.get_column(col)[aa] == '-' or aa in bads:
                protRecords[aa].map.append(0)
            else :
                try:
                    protRecords[aa].map.append(\
                        max(protRecords[aa].map)+1)
                except ValueError:
                    protRecords[aa].map.append(1)
    for col in range (0, alignment.get_alignment_length()):
        if '-' in alignment.get_column(col) or col in bads:
            for spe in alignment.get_all_seqs():
                spe.codons[spe.map[col]] = spe.codons[\
                    spe.map[col]].lower()
            alignment.map2ungap.append(0)
        else:
            try:
                alignment.map2ungap.append(max(alignment.map2ungap)+1)
            except ValueError:
                alignment.map2ungap.append(1)
    return alignment, protRecords

def _print_map(outMap,alignment,protRecords):
    '''
    for a given alignment object read from file, and the corresponding
    nucleotides sequences object (before alignment), create a map file
    '''
    out = open(outMap,'w')
    out.write('Python maps:\nmap2ungap\t'+\
                 '\t'.join(map(str,alignment.map2ungap[-100:]))+'\n\n')
    for spe in alignment.get_all_seqs():
        out.write('map '+spe.id+'\t')
        out.write('\t'.join(map(str,spe.map))+'\n')
    out.write('Aligned Sequences\' Map with/without Gaps\n\
    (UPPER CASE for columns without gaps)\n\
    (position in alignment: "including gaps"/"original"/"excluding \
gaps")\n\n')
    for record in protRecords:
        out.write(record.name+':'+'\n')
        out.write(' '*9+'1'+(' '*13) + \
                  ''.join\
                  (map(lambda i: str(i+1) + '/' + \
                       re.sub('^0','*',str (record.map[i]))+'/'+\
                       re.sub('^0','*',str (alignment.map2ungap[i]))+\
                       ' '*(20-(2+len (str(record.map[i])) \
                                + len (str(alignment.map2ungap[i]))\
                                + len (str (i+6)))), \
                       range(4,alignment.get_alignment_length(),5)))\
                  +'\n')
        out.write(' '*8+' |'+' '*14+(' |'+' '*18)*\
                  (alignment.get_alignment_length()/5)+'\n')
        out.write(' '*8+' '+'   '.join(map (lambda x: x, record.seq))\
                  +'\n')
        out.write(' '*8+' '.join\
                  (map(lambda x: record.codons[x],record.map))+'\n')
        out.write('\n')
    out.close()

def print_nucleic_alignment(alignment,keepGaps,out,format):
    '''
    print aligned sequences to file, in fasta, or paml format needs
    alignment with map, codons, features, keepGap True/False info,
    fileout info, and format (fasta, paml)  info.
    '''
    for i in range (0,len (alignment.get_all_seqs())):
        if keepGaps:
            alignment[i].seq = \
                             Seq(\
                re.sub(' ','-',''.join\
                       (map(lambda x: alignment[i].codons[x],\
                            alignment[i].map))).upper())
        else:
            alignment[i].seq=Seq(re.sub('[a-z ]','',\
                                        ''.join(alignment[i].codons)))
    if format == 'paml':
        outF = open (out,'w')
        outF.write('  '+str (len (alignment.get_all_seqs()))+' '+\
                   str(alignment.get_alignment_length())+'\n')
        for seq in alignment.get_all_seqs():
            outF.write('>'+seq.description+'\n')
            outF.write(seq.seq.tostring()+'\n')
        outF.close()
    else:
        try:
            h = open(out,'w')
            AlignIO.write([alignment],h,"fasta")
            h.close()
        except:
            sys.exit('ERROR :P problem writing alg: '+out+'\n')
    return alignment.get_alignment_length()

def check_seq(seq, rm_cod = True):
    '''
    check for format sequence
    '''
    seq = seq.upper()
    seq = re.sub('-','',seq)
    prot = ''
    seq, prot, codon_list = _try_translate(seq, rm_cod)
    if len(codon_list[-1])<3: del(codon_list[-1])
    return Seq(''.join(codon_list)), prot, codon_list

def trimalFasta(inFile):
    os.system('/home/francisco/tools/trimAl/source/trimal -in '\
              +inFile+' -fasta | sed "s/ [0-9][0-9]* bp//" > '\
              +inFile+'_ok')
    return inFile+'_ok'

def preserve(name):
    '''
    prevent to overwrite existing file by adding number at the end of new
    file name
    '''
    while os.path.isfile(name):
        if not 'count' in locals(): count = 0
        count += 1
        name = re.sub('_log.*','_log'+str(count),name)
    return name

def remove_extra_species(protRecords,nuclRecords,species):
    new_prots = []
    new_nucls = []
    for rec in protRecords:
        if rec.name in species:
            new_prots.append(rec)
    for rec in nuclRecords:
        if rec.name in species:
            new_nucls.append(rec)
    return new_prots, new_nucls

def fasta2recordLists(seqfile, rm_cod=False):
    '''
    convert fasta file into a list of protein records, and a list
    of nucleotide records.
    rm_cod argument is for removing not ATGC characters.
    '''
    seqFile = open(seqfile,'r')
    protRecords = []
    nuclRecords = []
    iterSeq = SeqIO.parse(seqFile, "fasta")
    for seqRec in iterSeq:
        seqRec.seq, prot, codon_list = \
                    check_seq(seqRec.seq.tostring(), rm_cod)
        for seq in nuclRecords:
            if seqRec.id == seq.id:
                seqRec.id = re.sub(' ', '_',seqRec.description)
                if seqRec.id == seq.id:
                    seqRec.id = seqRec.id + '.'
        protRecords.append(SeqRecord(prot, seqRec.id, seqRec.name,'|'.join(\
            seqRec.description.split('|')[0:4])))
        nuclRecords.append(seqRec)
    seqFile.close()
    return protRecords, nuclRecords


def main(args = ''):
    '''
    '''
    try:
        seqfile, outFile, keepGaps, format, rm_bad , \
                 rm_cod, species, trim  = getOptions(args)
    except TypeError:
        sys.stderr.write('ERROR:P missing infile or bad args\n')
        usage()
    seqfile  = trimalFasta(seqfile)
    outProts = re.sub('_ok', '_pep', seqfile)
    outMap   = re.sub('_ok', '_map', seqfile)
    inProts  = re.sub('_ok', '_tmp', seqfile)
    protRecords, nuclRecords = fasta2recordLists (seqfile, rm_cod)
    if not species == []:
        protRecords, nuclRecords = \
                     remove_extra_species (protRecords,nuclRecords,species)
    try:
        bads = run_Muscle(protRecords, nuclRecords,\
                          inProts, outProts, trim=trim, rm_bad=rm_bad)
    except:
        os.system('mv '+inProts+' '+outProts)
        raise
        sys.exit('ERROR :P running muscle. Moving input file to '\
                 +outProts+'\n')
    alignment, protRecords = alignment2map(nuclRecords, outProts, bads)
    _print_map(outMap,alignment,protRecords)
    length = print_nucleic_alignment(\
        alignment,bool(int(keepGaps)),outFile,format)

def getOptions(args=''):
    format      = "fasta"
    keepGaps    = False
    rm_bad      = False
    rm_cod      = False
    trim        = False
    outFile     = ''
    species     = []
    if args != '':
        sys.argv = ['']+args.split()
    try:                                
        opts,args =getopt.getopt\
                    (sys.argv[1:],"hi:o:krf:XNqs:t",\
                    ["help","inFile=","outFile=","keepsgaps",\
                     "format=","rmbad","Normal","species=","trim"])
    except getopt.GetoptError:
        usage()
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
        elif opt in ("-i","--inFile"):
            inFile = arg
        elif opt in ("-o", "--stdout"):
            outFile = arg
        elif opt in ("-k", "--keepgaps"):
            keepGaps = True
        elif opt in ("-f", "--format"):
            format = arg
        elif opt in ("-X", "--rmbad"):
            rm_bad = True
        elif opt in ("-N", "--Normal"):
            rm_cod = True
        elif opt in ("-s", "--species"):
            species = re.sub('\"','',arg).split('|')
        elif opt in ("-t", "--trim"):
            trim = True
        else:
            sys.stderr.write('... bad args')
            usage()
    try:
        if outFile == '': outFile = inFile+'_ali'
        return inFile, outFile, keepGaps, format, \
               rm_bad, rm_cod, species, trim
    except:
        sys.stderr().write('ERROR:P missing infile or bad arguments')
        usage()

def usage():
    sys.stderr.write('''
muscle_aligner.py v0.2b
bioinfo.cipf.es 09.2009
--------------------------

-h (--help):
     This help
-i (--inFile):
     (absolute) path to fasta input file with sequences to align
-o (--outFile):
      path and name to outfiles.
      eg: "-o lala" will generate:
            - lala_ali.fasta
            - lala_ali.map
            - lala_ali.pep
-k (--keepgaps):
      To keep the gaps generates by muscle in the alignment
-s (--species):
      Specie set you want to include e.g.:
      Homo_sapiens|Pan_troglodytes|Macaca_mulatta|Pongo_pygmaeus|Mus_musculus
-f (--format):
      Default is "fasta" format, but can also be "paml" ("-f paml")
-t (--trim):
      use trimAl to clean alignment with automated1 option
-N (--Normal):
      keep non-ATGC codons.
-X (--rmbad):
      Check if one sequence reduces consequently, the length of the
      ungapped alignment, with trimAl. And remove it (them).


                                              _  _(o)_(o)_  _
                                            ._\`:_ F S M _:' \_,
                                                / (`---'\ `-.
                                             ,-`  _)    (_, 

''')
    sys.exit()

if __name__ == "__main__":
    sys.exit(main())
