
import os
import sys
import pkg_resources
import argparse
import scripts.blastcall as bpc
import scripts.mousehomo as mh
import scripts.hpb as hb
import scripts.cpgis as cpg
import scripts.calcprop as prop
import scripts.id_map as im
import scripts.snp_oligo as so
import scripts.filter as fo
import scripts.check_file as cf
from Bio.Seq import Seq
import pandas as pd
#Split sequence to length 20
def getoligo(seq1):
    oarr=[]
    blhead=[]
    oname=[]
    start=[]
    end=[]
    seqpos=len(seq1)
    for i in range(0,len(seq1)):
        i=int(i)
        j=i+20
        j=int(j)
        sposn=seqpos - i
        eposn= sposn-20
        oarray=seq1[i:j]
        if(len(oarray)==20):
            start.append(sposn)
            end.append(eposn)
            blhead.append(">oligo"+str(i) + '\n' + oarray+'\n')
            oarr.append(oarray)
            oname.append("oligo"+str(i) )
    return(oarr,blhead,oname,start,end)
# Reverse complement to generate antisense oligonucleotide
def findcomp(seq):
    comp = ""
    for base in seq:
        if base == "A" :
            comp += "T"
        elif base == "T" :
            comp += "A"
        elif base == "G" :
            comp += "C"
        elif base == "C" :
            comp += "G"
        else :
            print("check seq")
            comp = None
            break
    rev=comp[::-1]
    return(rev)

def getinp(inp):
    try:
        infile = open(inp, 'r')
        for seq in infile.readlines():
            if not '>' in seq :
                seq=seq.replace(' ','')
                #print(seq)
    except IOError: 
        print("File doesnot exist")
        sys.exit()
    infile.close()
    return(seq)

if __name__ == '__main__':
    dir_name="input" #Directory to store input file  
    ## Input from user
    parser = argparse.ArgumentParser(description='Store input files in folder path input/')
    parser.add_argument('--infile', required=True,help='Fasta file to generate ASOs')
    parser.add_argument('--Exclude', required=False,help='Accession numbers to ignore during off-target search with blast')
    args = parser.parse_args()
    fname=args.infile
    ename=args.Exclude
    inp=os.path.join(dir_name,fname)
    exc_path=os.path.join(dir_name,ename) ## Excluding Accession nos of the transcript to which ASOs are designed
    ## To get input sequence from file to design oligo and split it to length of 20 bp
    seq=getinp(inp) 
    ## Reverse complementing sequence 
    revseq=findcomp(seq)
    ## Design Oligos
    oarr,blhead,oname,start,end= getoligo(revseq) # blhead- Oligos in fasta format to use in blast
    Es=[535837,536243,537013,541346] # Exon start position
    Ed=[535923,536340,537147,542087] # Exon end position
    ## Calculate gcpercentage & cpg island of Asos
    islarr,gcp = cpg.getcpg(oarr)
    ## Calculate Molecular weight and melting temperature
    mw,tm=prop.getprop(oarr)
    ## Predict folding
    stru,mfe=hb.getstr(oarr)
    ## Off-target identification
    ontar=bpc.call_blast(blhead,oarr,oname,exc_path)
    ## Mouse homology
    hom=mh.call_homo(blhead,oarr,oname)
    ## ID mapping for SNP analysis
    narr,oarr,gs,ge,s2,e2=im.imap(oarr,Es,Ed,seq)
    ## Population frequency
    
    oname,csnp,pop=so.getsnp(oname,narr,oarr,gs,ge,s2,e2)
    
    df = pd.DataFrame()
    df['Name']=oname
    df['ASO']=oarr
    df['start']=start
    df['end']=end
    df['cpg island']=islarr
    df['GC percentage']=gcp #G/C content: a 45 to 65% or greater G/C content will facilitate stronger ASO binding
    df['Molecular weight']=mw
    df['Melting Temperature']=tm
    df['Base associablity'] = stru #small hairpin loop of four unpaired bases can initiate antisense activity and  a segment  of four consecutive bases is entirely single-stranded. ye Ding 2002.
    df['BA score']=mfe
    df['Target effect'] = ontar
    df['Mouse homology']=hom
    df['snp']=csnp
    df['population frequency']=pop
    print(df)
    ##Filtering antisense oligonucleotides
    fo.filteroligo(df)
    ## Dimer formation
    
    ## Dataframe to store oligo properties 
    dfile_path='output/oligo_analysis_master.out'
    cf.check_file(dfile_path)
    df.to_csv(dfile_path, sep='\t')
    df_filtered=fo.filteroligo(df)
    farr=list(df_filtered['ASO'])
    #Check for dimer formation of ASOs
    sarr2,stru2,mfe2=hb.getstr2(farr)
    dimer_path='output/dimer_filtered.out'
    cf.check_file(dimer_path)
    out1=open(dimer_path,"w")
    for d in range(len(sarr2)):
        print(str(sarr2[d])+"\t"+str(stru2[d])+"\t"+str(mfe2[d])+"\n")
        out1.write(str(sarr2[d])+"\t"+str(stru2[d])+"\t"+str(mfe2[d])+"\n")
    out1.close()
    ## Filtered dataframe
    dfiltered_path='output/oligo_analysis_filtered.out'
    cf.check_file(dfiltered_path)
    df.to_csv(dfiltered_path, sep='\t')



    


  

    

    
