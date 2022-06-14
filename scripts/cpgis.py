def getcpg(oarr):
    gcp=[]
    cpg=[]
    isl=[]
    for i in range(0,len(oarr)):
        c=oarr[i].count("C")
        t=oarr[i].count("T")
        g=oarr[i].count("G")
        a=oarr[i].count("A")
        cg=c*g
        ocpg=oarr[i].count("CG") #Observed CPGs
        if cg >1:
            ecpg= c*g/20 # Expected CPGs
            cpg=(ocpg/ecpg)*100
        else:
            cpg=0
        gcp.append((c+g)/(a+t+g+c)*100)
        isl.append(cpg)
       
    return(isl,gcp)