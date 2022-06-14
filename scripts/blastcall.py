import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import os
from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
import scripts.check_file as cf
def call_blast(blhead,oarr,oname,neg_path):
    ontar=[]
    tardb=[]
    my_db = "db/human/hg38_cd34nonoverlap.trans" 
    my_blast_exe = "tools/ncbi-blast-2.13.0+/bin/blastn"
    ofile_path = 'output/blout.out'
    cf.check_file(ofile_path)
    # check if size of file is 0
    cf.check_empty(ofile_path) 
    tfile_path='output/bltest.out'
    out=open(ofile_path,"a")
    for i in range(0,len(blhead)):
        my_query=blhead[i]
        cf.check_file(tfile_path)
        inp=open(tfile_path,"r")
        out=open(ofile_path,"a")
        blast_command = Bio.Blast.Applications.NcbiblastnCommandline(cmd=my_blast_exe,out=tfile_path,task='blastn-short',word_size=7,dust='no',db=my_db,negative_seqidlist=neg_path, outfmt='6 qseqid qstart qend sseqid sstart send pident qcovs evalue',evalue=0.01,perc_identity=85,max_target_seqs=1,num_threads=6)
        stdout, stderr = blast_command(stdin=my_query)
        if os.stat(tfile_path).st_size == 0:
            #print(oname[i])
            out.write(oname[i]+ "\t"+oarr[i] + "\tno target effect\n")
            ontar.append(oarr[i])
            tardb.append("No target effect")
        else:
            #print(oname[i])
            for lines in inp.readlines():
                out.write(lines)
                tardb.append("offtarget found")
        inp.close()
        out.close()  
    os.remove(tfile_path)
    return(tardb)

