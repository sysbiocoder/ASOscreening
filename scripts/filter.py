
import pandas as pd
def filteroligo(df):
    df_fil1 = df[df['GC percentage'].isin(range(45,65))] 
    #print(df_fil1)
    df_fil2=df_fil1[df_fil1['Target effect']== 'No target effect'] 
    print(df_fil2)
    df_fil3=df_fil2[df_fil2['Mouse homology']== 'hits found']
    print(df_fil3)
    import re
    special = r"(((("  
    special2 = r"."  
    #Atleast 4 consecutive single strand bases
    df_fil4=df_fil3[df_fil3['Base associablity'].str.count(re.escape(special)) >= 1 ]
    print(df_fil4.shape)
    # hairpin loops
    df_fil5=df_fil4[df_fil4['Base associablity'].str.count(re.escape(special2))>=4]
    print(df_fil5)
    df_fil6=df_fil5[df_fil5['BA score']<0 ] #negative indicates stability
    print(df_fil6) 
    return(df_fil6)
'''
if __name__ == '__main__':
    df=pd.read_csv('../output/oligo_analysis_master.out', sep='\t')
    fil=filteroligo(df)

'''
    