import oligo as og
import array
import scripts.check_file as cf
def imap(oarr,Es,Ed,aseq): 
    sarr,alhead,aname,astart,aend= og.getoligo(aseq)
    arr=sarr[::-1]
    gln=len(Es)-1 # Exonic genomic positions array length
    oln=len(aseq) # oligos  array length
    gend = Ed[gln]
    gstart=gend-20 
    ge = Ed[gln] + 1
    
    e2array=[]
    s2array=[]
    gsarray=[]
    gearray=[]
    
    for i in range(len(arr)):
        e2=oln-i
        e2array.append(e2)
        s2=e2-20
        s2array.append(e2-20+1)
        ## Testing oligo length w.r.t genomic positions
        oglen=gend-gstart
        #print(gstart,Es[gln],gln)
        if gstart > Es[gln]:
            ge=ge-1
            gearray.append(ge)
            gs=ge-20 +1
            gsarray.append(gs)
            #print(arr[i],oarr[i],gs[i],ge[i],s2[i],e2[i])
            gstart=gs
        else:
            oglen2 = ge - gstart
            if gln > 0:
                gen=Ed[gln-1]
                if oglen2 > 0:
                    gs= gen -20 +oglen2 +1 
                    ge = Es[gln] + oglen2-1
                    gearray.append(ge)
                    gsarray.append(gs)
                    #print(arr[i],oarr[i],gs,ge,s2,e2)
                    oglen2=oglen2-1
                else:
                    gln=gln-1
                    ge=Ed[gln]
                    gs=ge-20 + 1
                    gearray.append(ge)
                    gsarray.append(gs)
                    gstart=gs
        #print(arr[i],oarr[i],gsarray[i],gearray[i],s2array[i],e2array[i])
    ofile_path="output/Idmapped.txt"
    out=open(ofile_path,"w")
    cf.check_file(ofile_path)

    for t in range(len(gsarray)):
        out.write(str(aname[t]) + "\t"+str(arr[t]) + "\t"+ str(oarr[t]) + "\t"+ str(gsarray[t]) +"\t"+str(gearray[t])+"\t"+str(s2array[t])+"\t"+str(e2array[t]) +"\n")
    out.close()
 
    return(arr,oarr,gsarray,gearray,s2array,e2array)



        
