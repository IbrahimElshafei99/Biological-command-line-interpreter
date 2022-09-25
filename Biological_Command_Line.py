import sys
import getopt
import Bio
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Seq import UnknownSeq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
##################
####### GC #######
##################
def gc(dna):
  gc_percentage=GC(dna)
  print('GC = {}'.format(gc_percentage))
##########################
####### Transcribe #######
##########################
def Transcribe(DNA):
  dna = Seq(DNA)
  transcribe=dna.transcribe()
  print('{}'.format(transcribe))
###################################
####### Reverse_complementt #######
###################################
def Reverse_comp(DNA):
  dna = Seq(DNA)
  reverse_complementt=dna.reverse_complement()
  print(reverse_complementt)
###########################
####### Calc_nbases #######
###########################
def Calc_nbases(DNA):
  dna=Seq(DNA)
  calc_nbases=dna.count('N')
  print(calc_nbases)
##############################
#######  Filter_nbases #######
##############################
def Filter_nbases(DNA):
  dna=Seq(DNA)
  print(dna.replace('N',''))
'''
def filter_nbases(sequence):
    sequence2=' '
    N_base='nN'
    for base in sequence:
        if base not in N_base:
            sequence2=sequence2+base
    return sequence2
                
print(filter_nbases('ACCNACNN'))
'''
#########################
######## is_valid #######
#########################
def is_valid(sequence,typee):
    dna_sequance='ACGT'
    rna_sequance='ACGU'
    protein_sequance='ABCDEFGHIKLMNPQRSTVWXYZabcdefghiklmnpqrstvwxyz'

    valid=True
    if typee=='DNA' or typee=='dna':
        for base in sequence:
            if base not in dna_sequance: 
                valid=False  
        return valid
        
    elif typee=='RNA' or typee=='rna':
        for base in sequence:
            if base not in rna_sequance:
                valid=False 
        return valid
        
    elif typee=='PROTEIN' or typee=='protein':
        for base in sequence:
            if base not in protein_sequance: 
                valid=False  
        return valid
    print('The type you entered is incorrect')
######################################################
####### pairwise alignment between 2 sequances #######
######################################################
def seq_alignment(seq1,seq2,file_name):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    if len(file_name)==0:
        for alignment in alignments:
            print(pairwise2.format_alignment(*alignment))
    else:
        f = open(file_name, 'a')
        for alignment in alignments:
            f.write(pairwise2.format_alignment(*alignment))
        f.close
#####################################################################
######## pairwise alignment between 2 sequances in text files #######
#####################################################################
def seq_alignment_files(File1Path, File2Path, OutputPath):
    sequance1=''
    sequance2=''
    file1 = open(File1Path, "r")
    for line in file1:
        line=line.rstrip()#this discard the newline at the end (if any)
        if line[0]=='>': #if file is fasta
            sequance1=line[1:]
        else :
            sequance1 = line[0:]
    file1.close
    
    file2 = open(File2Path, "r")
    for line in file2:
        line=line.rstrip()#this discard the newline at the end (if any)
        if line[0]=='>': #if file is fasta
            sequance2 = line[1:]
        else :
            sequance2 = line[0:]
    file2.close

    alignments = pairwise2.align.globalxx(sequance1,sequance2)
    if len(OutputPath)==0:
        for alignment in alignments:
            print(pairwise2.format_alignment(*alignment))
    else:
        f = open(OutputPath, 'a')
        for alignment in alignments:
            f.write(pairwise2.format_alignment(*alignment))
        f.close
################################
####### online_alignment #######
################################
def online_alignment(seq, FileName):
    result_handle = NCBIWWW.qblast('blastn', 'nt', seq)
    blast_record = NCBIXML.read(result_handle)
    threshold = 0.01
    if len(FileName) == 0:
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
              if hsp.expect < threshold:
                print('sequence: ', alignment.title)
                print('length: ', alignment.length)
                print('e value: ', hsp.expect)
                print(hsp.match)
                print(hsp.query)
                print(hsp.sbjct)
    else:
        f=open(FileName, 'a')
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
              if hsp.expect < threshold:
                f.write('sequence: ')
                f.write(alignment.title)
                f.write('\n')
                f.write('length: ')
                f.write(str(alignment.length))
                f.write('\n')
                f.write('e value: ')
                f.write(str(hsp.expect))
                f.write('\n')
                f.write(hsp.match)
                f.write('\n')
                f.write(hsp.query)
                f.write('\n')
                f.write(hsp.sbjct)
        f.close
################################
####### convert_to_fasta #######
################################
def convert(file):
    
   with open(file) as input_handle, open( file+".fasta", "w" ) as output_handle: 
    sequences = SeqIO.parse(input_handle, "genbank") 
    count = SeqIO.write(sequences, output_handle, "fasta") 
    print("Converted %i records" % count)
###########################
####### merge_fasta #######
###########################
def merge ( *filenames, outFile):

    if len(outFile) == 0:
        for names in filenames:
            with open(names) as infile:
                print(infile.read()+"\n")
    else:
        with open(outFile, "w") as outfile:
            for names in filenames:
                with open(names) as infile:
                    outfile.write(infile.read())
                    outfile.write("\n")
def test_getopt():
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], 'o:h:')
    except getopt.GetoptError as err:
        print(err)
        opts=[]

    for ARG in args:
      if ARG == 'gc':
        if len(args) < 2:
          print('Missing parameters')
        else:
          gc(args[1])
          
      elif ARG == 'transcribe': #transcribe
        if len(args) < 2:
          print('Missing parameters')
        else:
          Transcribe(args[1])
          
      elif ARG == 'reverse_complement': #reverse_complement
        if len(args) < 2:
          print('Missing parameters')
        else:
          Reverse_comp(args[1])
          
      elif ARG == 'calc_nbases': #calc_nbases
        if len(args) < 2:
          print('Missing parameters')
        else:
          Calc_nbases(args[1])
          
      elif ARG == 'filter_nbases': #filter_nbases
        if len(args) < 2:
          print('Missing parameters')
        else:
          Filter_nbases(args[1])
          
      elif ARG == 'is_valid': #is_valid
        if len(args) < 3:
          print('Missing parameters')
        else:
          print(is_valid(args[1], args[2]))
        
      elif ARG == 'seq_alignment': #seq_alignment
        if len(args) < 3:
          print('Missing parameters')
        else:
          if len(opts) == 0:
            seq_alignment(args[1], args[2], '')
        
      elif ARG == 'seq_alignment_files': #seq_alignment_files
        if len(args) < 3:
          print('Missing parameters')
        else:
          if len(opts) == 0:
            seq_alignment_files(args[1], args[2], '')
        
      elif ARG == 'online_alignment': #online_alignment
        if len(args) < 2:
          print('Missing parameters')
        else:
          if len(opts) == 0:
            online_alignment(args[1], '')
        
      elif ARG == 'convert_to_fasta': #convert_to_fasta
        if len(args) < 2:
          print('Missing parameters')
        else:
          convert(args[1])
        
      elif ARG == 'merge_fasta': #merge_fasta
        if len(args) < 3:
          print('Missing parameters')
        else:
          if len(opts) == 0:
            merge(*args[1:], outFile='')
        
      else:
        print("Enter the right command name")
      break
      
      
    for opt,ARG in opts:
        if opt in ['-o', '-h']:
            if args[0]=='seq_alignment': #seq_alignment
              if len(ARG)!=0:
                seq_alignment(args[1], args[2], ARG)
              else:
                print('option should be like "output.txt"')
                
            elif args[0]=='seq_alignment_files': #seq_alignment_files
              if len(ARG)!=0:
                seq_alignment_files(args[1], args[2], ARG)
              else:
                print('option should be like "output.txt"')
                
            elif args[0]=='online_alignment': #online_alignment
              if len(ARG)!=0:
                online_alignment(args[1], ARG)
              else:
                print('option should be like "output.txt"')
                
            elif args[0]=='merge_fasta': #merge_fasta
              if len(ARG)!=0:
                merge(*args[1:], outFile=ARG)
              else:
                print('option should be like "output.txt"')
        else:
            assert False, "unhandled option"
    #print('args: {}'.format(args))
    #print('opts: {}'.format(opts))

test_getopt()
