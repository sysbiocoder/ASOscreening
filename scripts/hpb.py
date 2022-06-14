import RNA
def getstr2(sarr2):
    settings = RNA.md()
    stru1=[]
    scr1=[]
    stru2=[[]]
    scr2=[[]]
    for i in range(0,len(sarr2)-1):
        for j in range(i+1,len(sarr2) ):
            seqstr=sarr2[i] + '&' + sarr2[j]
            print(seqstr)
            fc_obj2 = RNA.cofold(seqstr)
            (structure2,mfe2) = fc_obj2
            if mfe2 < 0:
                stru1.append(structure2)
                scr1.append(mfe2)
        stru2.append(stru1)
        scr2.append(scr1)
    return(sarr2,stru2,scr2)
def getstr(sarr):
    settings = RNA.md()
    stru=[]
    scr=[]
    for i in range(0,len(sarr)):
        fc_obj = RNA.fold_compound(sarr[i])
        structure, mfe = fc_obj.mfe()
        #print(structure,mfe)
        stru.append(structure)
        scr.append(mfe)
    return(stru,scr)