import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import os
from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
import scripts.check_file as cf
def call_homo(blhead,oarr,oname):
    ontar=[]
    homdb=[]
    my_db = "db/mouse/M29trans" 
    my_blast_exe = "tools/ncbi-blast-2.13.0+/bin/blastn" 

    ofile_path = 'output/mblout.out' #Blast output file for all ASOs
    cf.check_file(ofile_path)
    cf.check_empty(ofile_path)     # check if size of file is 0
    tfile_path='output/mbltest.out' #Temporary file to store blast output overwritten for each oligo
    for i in range(0,len(blhead)):
        my_query=blhead[i]
        blast_command = Bio.Blast.Applications.NcbiblastnCommandline(cmd=my_blast_exe,out=tfile_path,task='blastn-short',word_size=7,dust='no',db=my_db, outfmt='6 qseqid qstart qend sseqid sstart send pident qcovs evalue',evalue=0.01,perc_identity=85,max_target_seqs=1,num_threads=6)
        cf.check_file(tfile_path)
        inp=open(tfile_path,"r")
        out=open(ofile_path,"a")
        stdout, stderr = blast_command(stdin=my_query)
        if os.stat(tfile_path).st_size == 0:
            out.write(oname[i]+ "\t"+oarr[i] + "\tno homolog\n")
            ontar.append(oarr[i])
            homdb.append("no homolog")
            #print(oname[i],"nohom")
        else:
            for lines in inp.readlines():
                out.write(lines)
                homdb.append("hits found")
                #print(oname[i],"hits")
    
        out.close()
        inp.close()
    os.remove(tfile_path)
    return(homdb)

