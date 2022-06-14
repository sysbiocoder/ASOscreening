# ENST00000606065_CDC34-004_19_535837-542087
import re
import scripts.check_file as cfs
def getsnp(oname,narr,oarr,gs,ge,s2,e2):
    csnp=[]
    pop=[]
    for i in range(len(gs)):
        inp_path="db/dbsnp/CDC34.vcf"
        cfs.check_files(inp_path)
        inp2=open(inp_path,"r")
        #print("\n",oname[o], oarr[o])
        flag=0
        m=gs[i]
        n=ge[i]
        snppos,ref,alt,resgp=[],[],[],[]
        cstart,cend,gstart,gend=[],[],[],[]
        Es=[535837,536243,537013,541346] # Exon start position
        Ed=[535923,536340,537147,542087] #Exon end position
        cs= [1,88,186,723]# Exon start position
        cd= [87,185,722,1043] # Exon end position
        for s in inp2.readlines():
            if not('##' in s):
                col1= s.split("\t")[0] #chr
                col2 = s.split("\t")[1] #pos
                col3= s.split("\t")[2] #snpid
                col4= s.split("\t")[3] #ref
                col5= s.split("\t")[4] #alt
                #print(m,n,col2)
                if(int(col2) in range(int(m),int(n))): #Check if snp lies within ASO
                    if re.search(r'(CAF=\d*\.\d*\,\d*\.\d*);', s):
                        res=re.search(r'(CAF=\d*\.\d*\,\d*\.\d*);', s)
                        #print(m,n,col2)
                        flag=flag+1
                        snppos.append(col2)
                        ref.append(col4)
                        alt.append(col5)
                        resgp.append(res.group(1))
                        cstart.append(s2[i])
                        cend.append(e2[i])
                        gstart.append(gs[i])
                        gend.append(ge[i])
                    else:
                        #print(m,n,'CAF not found',col1,col2,col3,col4,col5,';')
                        flag=flag+0
                else:
                    #print('pos not found in vcf',col1,col2,col3,col4,col5,m,n,';')
                    flag=flag+0
        if flag >= 1:
            #print("\n",oname[i], oarr[i],gstart,gend,cstart,cend, resgp,snppos, ref,alt)
            csnpstr=[]
            k=1
            for l in range(len(ref)):
                for r in range(len(cs)):
                    if(cstart[l] in range(cs[r],cd[r])):
                        k=r
                        break

                #if(cstart[l] in range(1,87) ):
                #    k=0
                #if(cstart[l] in range(88,185) ):
                #    k=1
                #if(cstart[l] in range(186,722) ):
                #    k=2
                #if(cstart[l] in range(723,1043) ):
                #    k=3
                if ((int(gstart[l]) >= int(Es[k])) & (int(gend[l]) <= int(Ed[k])) & (int(snppos[l]) <= int(Ed[k]))):
                    csnppos=cend[l]-(gend[l]-int(snppos[l]))
                    csnpstr.append("c."+ref[l]+str(csnppos)+alt[l])
                elif ((int(gstart[l]) <= int(Es[k])) & (int(gend[l]) <= int(Ed[k]) ) & (int(snppos[l]) >= int(Es[k]))) :
                    csnppos=int(cend[l])-(gend[l]-int(snppos[l]))
                    csnpstr.append("c."+ref[l]+str(csnppos)+alt[l])
                            #print(str(oname[i])+"\t"+ str(cend[l])+"\t"+ str(cstart[l])+"\t"+ str(gstart[l]) +"\t"+ str(gend[l]) +"\t"+str(snppos[l]))
                            #64 83 535904 536242  536154 535837-535923 536243
                elif ((int(snppos[l]) <= int(Es[k])) & (int(gstart[l]) <= int(Es[k])) & int(gend[l])< int(Ed[k])):
                    csnppos=int(cstart[l])+(int(snppos[l])-int(gstart[l]))
                    csnpstr.append("c."+ref[l]+str(csnppos)+alt[l])
                else:
                    print("no snp")
            pop.append(resgp)
            csnp.append(csnpstr)
            #print(csnp)
            flag=0
        else:
            csnp.append('NA')
            pop.append('NA')
            #print(oname,csnp,pop)
            #print(csnp)
    ofile_path="output/snpout.txt"
    out=open(ofile_path,"w")
    cfs.check_file(ofile_path)
    for t in range(len(csnp)) : 
        #print(len(csnp))      
        #print(oname[t],csnp[t],pop[t])
        out.write(str(oname[t])+"\t"+str(oarr[t])+"\t"+str(csnp[t])+"\t"+str(pop[t])+"\n")
    out.close()
    #out.write(narr,oarr,gsarray,gearray,s2array,e2array)
    return(oname,csnp,pop)
         
