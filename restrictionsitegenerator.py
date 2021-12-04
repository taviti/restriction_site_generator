from Bio.Seq import Seq
import numpy as np
import pandas as pd
from pandas import ExcelWriter





dic_rd  = {'AclI': 'AACGTT', 'ApaI': 'GGGCCC', 'HindIII': 'AAGCTT', 'SspI': 'AATATT', 'MluCI': 'AATT', 'PciI': 'ACATGT', 'AgeI': 'ACCGGT', 'MluI': 'ACGCGT', 'HpyCH4IV': 'ACGT', 'HpyCH4III': 'ACNGT', 'SpeI': 'ACTAGT', 'BglII': 'AGATCT', 'AfeI': 'AGCGCT', 'AluI': 'AGCT', 'StuI': 'AGGCCT', 'ScaI': 'AGTACT', 'ClaI BspDI': 'ATCGAT', 'NsiI': 'ATGCAT', 'AseI': 'ATTAAT', 'SwaI': 'ATTTAAAT', 'MfeI': 'CAATTG', 'Nb.BssSI': 'CACGAG', 'PmlI': 'CACGTG', 'PvuII': 'CAGCTG', 'NdeI': 'CATATG', 'NlaIII': 'CATG', 'NcoI': 'CCATGG', 'SmaI': 'CCCGGG', 'SacII': 'CCGCGG', 'MspI HpaII': 'CCGG', 'NciI': 'CCSGG', 'AvrII': 'CCTAGG', 'Nb.BbvCI': 'CCTCAGC', 'SbfI': 'CCTGCAGG', 'PvuI': 'CGATCG', 'BstUI': 'CGCG', 'EagI': 'CGGCCG', 'BsiWI': 'CGTACG', 'BsmBI-v2': 'CGTCTC', 'BfaI': 'CTAG', 'XhoI PaeR7I': 'CTCGAG', 'PstI': 'CTGCAG', 'AflII': 'CTTAAG', 'Nb.BsmI': 'GAATGC', 'EcoRI': 'GAATTC', 'AatII': 'GACGTC', 'Eco53kI': 'GAGCTC', 'EcoRV': 'GATATC', 'MboI Sau3AI DpnII': 'GATC', 'Nb.BsrDI': 'GCAATG', 'Nb.BtsI': 'GCAGTG', 'SphI': 'GCATGC', 'SrfI': 'GCCCGGGC', 
'NaeI': 'GCCGGC', 'AsiSI': 'GCGATCGC', 'HinP1I': 'GCGC', 'BssHII': 'GCGCGC', 'NotI': 'GCGGCCGC', 'NheI': 'GCTAGC', 'BamHI': 'GGATCC', 'HaeIII': 'GGCC', 'FseI': 'GGCCGGCC', 'SfoI': 'GGCGCC', 'AscI': 'GGCGCGCC', 'PspOMI': 'GGGCCC', 'KpnI': 'GGTACC', 'CviQI': 'GTAC', 'BstZ17I': 'GTATAC', 'SalI': 'GTCGAC', 'ApaLI': 'GTGCAC', 'HpaI': 'GTTAAC', 'PmeI': 'GTTTAAAC', 'SnaBI': 'TACGTA', 'BspHI': 'TCATGA', 'BspEI': 'TCCGGA', 'TaqI-v2': 'TCGA', 'NruI': 'TCGCGA', 'XbaI': 'TCTAGA', 'BclI': 'TGATCA', 'HpyCH4V': 'TGCA', 'FspI': 'TGCGCA', 'MscI': 'TGGCCA', 'BsrGI': 'TGTACA', 'MseI': 'TTAA', 'PacI': 'TTAATTAA', 'PsiI-v2': 'TTATAA', 'BstBI': 'TTCGAA', 'DraI': 'TTTAAA'}

# to translate given seq 



def translate1():

    
    print('****************************')
    userinpt = input('Enter input seq = ').upper()
    prompt1 = userinpt[1:(len(userinpt))]
    prompt2 = userinpt[2:(len(userinpt))]

    aa1 = Seq(userinpt).translate()
    aa2 = Seq(prompt1).translate()
    aa3 = Seq(prompt2).translate()

    translate1.a1=aa1
    translate1.a2=aa2
    translate1.a3=aa3

    print('****************************')
    print('ORF1:',aa1,'\n','ORF2:',aa2,'\n','ORF3:',aa3 )

    return userinpt, aa1,aa2,aa3

    

def chooseorf(aa1, aa2, aa3):

   
    while True:
        try:
            print('****************************')
            userinput2 = int(input('Please choose the peptide (Enter 1, or 2 or 3) = ')) 
            if userinput2 == 1:
                aa = aa1
                break          
            elif userinput2 == 2:
                aa = aa2
                break
            elif userinput2 == 3:
                aa = aa3
                break
            else: 
                print("Entered orf number is not Valid")
                break
        except:
            continue

  

    chooseorf.az = aa 

   
    print('choosen aa = ', aa)
    return aa   #userinpt

# to create Alanine/glycine mutation in translated seq
def muatate_aa(aa):
    print('****************************')
    for num, char in enumerate(aa):
        print(num, char)
    
    print('****************************')

    aanum = int(input("Enter the aminoacid number where you want to mutate or/and create Restriciotn site (Enter number only from +2 to -2 of sequence) = "))
    print('****************************')
    print("1: Alanine \n 2: Glycine,\n 3: None")
    print('****************************')
    mutateTo = int(input("Choose the aminoacid to replace (Enter number) =  "))

    
    if mutateTo == 1:
        aaseq = list(str(aa))
    
        mutation = ['A']
        aaseq[aanum] = mutation[0]  # to assign the 'A' mutation in place of input number
    
        #print(aanum) # index no of AA being mutated

        print ('Mutated aa Seq = ', aaseq)
        muatate_aa.mutatedseq = aaseq
        return aaseq, aanum
    
    elif mutateTo == 2:
        aaseq = list(str(aa))
        mutation = ['G']
        aaseq[aanum] = mutation[0]  # to assign the 'A' mutation in place of input number
    

        print(aanum) # index no of AA being mutated

        #print ('Mutated aa Seq = ', aaseq)
        muatate_aa.mutatedseq = aaseq
        return aaseq, aanum

    elif  mutateTo == 3:

        aaseq = list(str(aa))
        
        muatate_aa.mutatedseq = aaseq
        return aaseq, aanum


def list_chop(maaseq,aanum):

  
    fiveaalist = maaseq[aanum-2:aanum+3]
    
    #print('chopped list= ', fiveaalist)
    
    return fiveaalist


# to extract ntd seq for each list
def extract(lst):

    mydic = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
    
    extract = []
       
    for i in lst:
        for key, value in mydic.items():
            if value == i:
                extract.append(key)

   
    return extract # extracted ntd seq from given aa seq


# to get combinations of codons for 5 aa peptide
def arraycombination(five_aalist,aanum, userinptseq):

  

    a = extract(five_aalist[0]) # this will extract list of codons for a indexed aminoacid, for e.g 'S'
    b = extract(five_aalist[1])
    c = extract(five_aalist[2])
    d = extract(five_aalist[3])
    e = extract(five_aalist[4])


    newarray = np.array(np.meshgrid(a,b,c,d,e)).T.reshape(-1,5) # generates a array for all possible combinations of codon from each extracted list
   
  
    arraytolist = newarray.tolist()
   
    
      
    print(len(arraytolist)) # prints total no. of combinations 

  
    listofstrings = [] # this will create list with each seq as string

    for sublist in arraytolist:
           listofstrings.append(''.join(sublist))

    

    sequences = []
    simi_scores = []
    frequency_scores =[]
    sites = []

    aan = (aanum*3)


    def fullstring(strg):
        if chooseorf.az == translate1.a1:
            full_strg = userinptseq[:aan-6]+ strg +userinptseq[aan+9:]
            return full_strg

        elif chooseorf.az== translate1.a2:
            full_strg = userinptseq[:aan-5]+ strg +userinptseq[aan+10:]
            return full_strg

        elif chooseorf.az== translate1.a3:
            full_strg = userinptseq[:aan-4]+ strg +userinptseq[aan+11:]
            return full_strg



    for strg in listofstrings:
        
        if find_site(strg):
            #print("Sequence = ", strg)
            full_strg = fullstring(strg)
            sequences.append(full_strg)           
            simi_scores.append(seqsimilarity(userinptseq, full_strg)[0])
            frequency_scores.append(seqsimilarity(userinptseq, full_strg)[1])
            sites.append(find_site(strg))       
   
    return strg, sequences, simi_scores,  frequency_scores, sites

    
def seqsimilarity(seq1, seq2):
     

    def mismatch(s1,s2):
   
        counter = 0

        for p in range(len(s1)):
            if s1[p] != s2[p]:
                counter += 1

        return counter


    mismatch_frequency = mismatch(seq1, seq2)


    similarityscore = (mismatch_frequency / len(seq1))*100

    

    #print("Seq_Similarity||", similarityscore, "%")

   

    return similarityscore, mismatch_frequency



# dictionary for restriction enzymes and sites
def find_site(strng):
     
     # re dictionary imported as dic_rd           
    re_sites_list = []
    enzy = {}
    
        
    for key, value in dic_rd.items():
        if value in str(strng):
            #print("||Enzyme and Cleavage site: ", key,",", value)
            re_sites_list.append(value)
            enzy[key] = value
            
    
    #print(re_sites_list) 
        
    # to count no. of restriction sites for each enzyme       
    for site in re_sites_list:
        if site in strng:
            #print("No. of Restriction sites for above each Enzyme = ", strng.count(site))
            pass
                              
        return  enzy 


def towriteexcel(column1, column2,column3, column4):

    dat = { 'Degenerate Sequences': column1, 'Distance %': column2, 'Number of bases changed': column3, 'Restriction Enzymes and sites': column4  }


    df = pd.DataFrame(data=dat)
    pd.set_option('colheader_justify', 'center')

    dfsort = df.sort_values(by = ['Distance %'], ascending = True)

    dfreset = dfsort.reset_index(drop=True)

    
    with ExcelWriter('results.xlsx') as writer:
        dfreset.to_excel(writer)



    return dfreset 



def main():

   userinpt1_seq, aa_seq1, aa_seq2, aa_seq3 = translate1()

   choosen_aa_seq = chooseorf( aa_seq1, aa_seq2, aa_seq3)

   mutated_aa_seq, userinpt3_aanum = muatate_aa(choosen_aa_seq)

   fiveaalist = list_chop(mutated_aa_seq, userinpt3_aanum )

   strng, column1, column2, column3, column4 = arraycombination(fiveaalist,userinpt3_aanum, userinpt1_seq)

   find_site(strng)

   towriteexcel(column1, column2, column3, column4)

   print("Sucsessful, please find the excel file 'results.xlsx'in same foleder")
   print('****************************')




main()



    

    






