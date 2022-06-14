def getprop(oarr):
    #Calculate molecular weight and melting temperature
    #Molecular Weight = (An x 329.21) + (Un x 306.17) + (Cn x 305.18) + (Gn x 345.21) + 159.0

    mw=[]
    tm=[]
    for i in range(len(oarr)):
        An=0
        Tn=0
        Cn=0
        Gn=0
        for base in oarr[i]:
            base=str(base)
            if base == "A" :
                An += 1
            elif base == "T" :
                Tn += 1
            elif base == "G" :
                Gn += 1
            elif base == "C" :
                Cn += 1
            else :
                print("Cannot calculate check oligos bases")
                break
            
        wt=((An * 329.21) + (Tn * 306.17) + (Cn * 305.18) + (Gn * 345.21) + 159.0)
        Tm= 64.9 +41*(Gn+Cn-16.4)/(An+Tn+Gn+Cn)
        mw.append(wt)
        tm.append(Tm)
    return(mw,tm)