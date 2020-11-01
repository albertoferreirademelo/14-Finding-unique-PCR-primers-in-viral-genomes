#Python 3.3.5 (32 bits)
#Libraries installed: - Biopython (biopython.org) and numpy

from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from datetime import datetime

def Open_file(name_of_file, type_of_file):
    """
    The Open_file function will open a file with the genetic code using biopython and it will return the sequence
    INPUT: file containing a genetic cod (eg. "example.fasta" and what kind of file it is (eg. "fasta")
    OUTPUT: The sequence of the genetic code from the file
    """     
    for seq_record in SeqIO.parse(name_of_file, type_of_file):
        #ID = (seq_record.id)        
        Sequence = (seq_record.seq)                
        Length = (len(seq_record))
        #print (Length)
                            
    return (Sequence, Length)

def Cut_sequence(sequence, size, forward):

    """
    The Cut_sequence function will open a sequence and it will do primers of size "size", then it will return all possible
    primers of that size
    INPUT: A genetic code (eg. "ATTACTACGCTGCA") and the size of the parts you want to cut (eg. 3)
    OUTPUT: All possible primers from this genetic code, in the example above if I would cut in the size 3 the output would look like
    this: [ATT, TTA, TAC, ..., TGC, GCA]
    """
    dictionary_all_primers = {}    
    all_primers = []
    location = 1
    test = []
    test.append(sequence)        
    for i in sequence:
        temp = sequence[location-1:location-1+size]
        all_primers.append(temp)
        dictionary_all_primers[temp] = location
        '''
        if forward == 1:
            dictionary_all_primers[temp] = location
        elif forward == 0:
            if i == "GCACCTATCGTGCCCTTATG":
                print ("LOCATION!!! !!!! !!! : ",location)
            dictionary_all_primers[temp] = length_sequence - location
        '''        
        location += 1        
    for i in range(size-1):
        del all_primers[-1]    
    return all_primers, dictionary_all_primers        

def Filter_unique(primers_list):
    """
    The Filter_unique function will compare all the primers from the primers list and return a list with primers 
    INPUT: All possible primers (can be used after the function Cut_sequence)
    OUTPUT: All primers that have no exact duplicates
    """
    unique = []    
    primers_dict = dict(Counter(primers_list))
    for k, v in primers_dict.items():
        if v == 1:
            unique.append(k)
    return (unique)
    
def Filter_Temp_and_GC(listan, temperature, GC_max, GC_min):
    """
    The Filter_temp_and_GC function will filter out all primers that does not have the specific temperature, this
    function can also filter out the GC quantity in % and it will filter out all primers that does not contain 
    a GC quantity between 40 and 60%
    INPUT: A list with the possible primers and a desired melting temperature, where A and T add 2 Celsius
    to the melting temperature and G and C add 4 Celsius.
    OUTPUT: A list of all primers that exactly at the melting temperature desired
    """        
    filtered_codons = []    
    for i in listan:
        temp_list = []        
        for j in i:
            if j == 'A':
                temp_list.append(2)                
            elif j == 'T':
                temp_list.append(2)
            elif j == 'G':
                temp_list.append(4)
            elif j == 'C':
                temp_list.append(4)         
        sum_temp = sum(temp_list)        
        if sum_temp == temperature:           
            GC_Percent = Filter_GC_percent(temp_list)           
            if (GC_Percent >= GC_min) and (GC_Percent <= GC_max):
                filtered_codons.append(i)               
    return (filtered_codons) 

def Filter_GC_percent(value):
    """
    This function will check GC content
    INPUT: list with temperatures (where A and T value is 2 and G and C value is 4)
    OUTPUT: GC content in percent
    """                
    return (Counter(value)[2]/float(len(value)))

def Filter_endC(listan):
    """
    This function will see if at the 3' end of the primer will be GC or CG because
    that will give a better PCR product
    INPUT: list with primers
    OUTPUT: list with primers that have CG or GC at 3' end
    """
    temp = []
    for i in listan:
        if ((i[-2::] == 'GC') or (i[-2::] == 'CG')):
            temp.append(i)
    return (temp)
            
def Filter_check_repetition(listan, max_allowed_aa_repetition = 4):
    """
    This function will filter out all primers that have aa repetitions (ex. aaaa or tttt)
    INPUT: list with primers, number of max allowed aa repetitions
    OUTPUT: list of primers without repetitions
    """
    max_allowed_aa_repetition = max_allowed_aa_repetition+1
    filtered_list = []
    for i in listan:
        if (("A"*max_allowed_aa_repetition in i) == False) and (("T"*max_allowed_aa_repetition in i) == False) and (("C"*max_allowed_aa_repetition in i) == False) and (("G"*max_allowed_aa_repetition in i) == False):
            filtered_list.append(i)  
    return (filtered_list)

def Filter_CG_5base_end(listan):
    """
    This function will check if it is not to many C or G at the last 5 nucleotides of the primer 
    since to many of C or G will give a to strong bond between the primer and the DNA strand
    INPUT: list with primers
    OUTPUT: list of primers without to many C or G at 3' end  
    """
    #nr = 0
    filtered_CG_5base = []
    for cod in listan:
        nr = 0 
        for i in cod[-5::]:
            if i == "C" or i == "G":
                nr += 1
        if nr <= 4:
            filtered_CG_5base.append(cod)
    return filtered_CG_5base

def Filter_hamming_codTemp(all_sequence, list_of_primers, forward, annealing_temp = 15):   
    """
    This function will compare the primers with the DNA with the help of the hamming distance
    algorithm. If the difference of a nucleotide is A or T the temperature of annealing will 
    rise 2 degree Celsius but if the nucleotide is C or G the temperature will rise 4 degre instead.
    If the total sum of the difference is less than the annealing temperature, then the primer will
    be removed from the list since the annealing can bind to another place.
    INPUT: List with all sequences of the DNA, list of primers, if the sequence is forward or reverse
    and the min allowed difference in the annealing temperature
    OUTPUT: Filtered list with primers that are unique in the annealing temperature
    """ 
    filtered_hamming = list(list_of_primers)        
    for cod in all_sequence:
        for i in list_of_primers:            
            total_sum = 0
            compare_cod = zip(cod, i)
            for ch1, ch2 in compare_cod:
                if total_sum > annealing_temp:
                    break                
                if ch1 != ch2:
                    if ch2 == "A" or ch2 == "T":                        
                        total_sum = total_sum+2
                    else:                        
                        total_sum = total_sum+4
            if forward == 1:
                if total_sum < annealing_temp and total_sum != 0:
                    if i in filtered_hamming: filtered_hamming.remove(i)
            else:
                if total_sum < annealing_temp:
                    if i in filtered_hamming: filtered_hamming.remove(i)                
    return (filtered_hamming)

def save_to_file(listan, name_of_file):
    """
    This function will save a list to a file
    INPUT: A list and the name of the file
    OUTPUT: A file with the list on it
    """
    text_file = open(name_of_file, "w")
    for i in listan:
        output = str(i+"\n")
        text_file.write(output)
    text_file.close()
    
def save_dic_to_file(dictionary, name_of_file):
    """
    This function will save a dictionary to a file
    INPUT: A dictionary and the name of the file
    OUTPUT: A file with the dictionary on it
    """
    text_file = open(name_of_file, "w")
    for i in dictionary:
        output = str("Sequence: "+str(i)+" location: "+str(dictionary[i])+" Complementary: "+Do_complement(i)+"\n")
        text_file.write(output)
    text_file.close()
    
def hamming(cod, complement):
    """
    This function will check the sum of the difference of two sequences
    (also known as the hamming distance). The sequences need to have the
    same size
    INPUT: Two sequences of the same size
    OUTPUT: Sum of difference of the two sequences
    """
    total_sum = 0
    comparissom = zip(cod,complement)
    for ch1,ch2 in comparissom:
        if ch1 != ch2:
            total_sum += 1
    return (total_sum)    
    
def Filter_hairpin(listan):
    """
    This function will check if the primer can bind to itself
    INPUT: List with primers
    OUTPUT: Primers that will not bind to itself
    """
    filtered_hairpin = []
    for cod in listan:
        complement = Do_complement(cod)[::-1]
        test = hamming(cod, complement)
        if test > len(cod)-3: 
            filtered_hairpin.append(cod)
    return filtered_hairpin        
    
def Do_complement(cod):
    """
    This function will do a complement of a sequence
    INPUT: A DNA sequence
    OUTPUT: The complement of the DNA sequence
    """
    complement = []
    for i in cod:
        if i == "A":
            complement.append("T")
        if i == "T":
            complement.append("A")
        if i == "C":
            complement.append("G")
        if i == "G":
            complement.append("C")
    complement = "".join(complement)
    return complement

def Do_reverse(listan):
    """
    This function will do the reverse of a list of sequences
    INPUT: List with sequences
    OUTPUT: List with the reverse of the sequences
    """
    reversed_list = []
    for i in listan:
        reversed_list.append(i[::-1])
    return (reversed_list)

def Do_all_filters(Sequence,specie,Forward_or_Reverse, lenght, temp, length_sequence, GC_max = 60, GC_min = 40, Pair_max = 1500, Pair_min = 300,max_allowed_aa_repetition = 4, annealing_temp = 15):
    """
    This function will run all filters
    """
    print (Forward_or_Reverse)
    Possible_primers_forward, dictionary_all_primers = Cut_sequence(str(Sequence), lenght, 1)    
    #save_to_file(Possible_primers_forward, specie+"_AllPrimers_"+Forward_or_Reverse+".txt")    
    Possible_primers_reverse = Cut_sequence(Do_complement(str(Sequence[::-1])), lenght, 0)
    print ("List of all possible primers:",len(Possible_primers_forward))    
    All_uniques = Filter_unique(Possible_primers_forward)        
    Filter_temp_GC = Filter_Temp_and_GC(All_uniques,temp, GC_max, GC_min)    
    Filter_CG_end = Filter_endC(Filter_temp_GC)
    print ("After Temp and GC filter:",len(Filter_CG_end))
    Filtered_with_hairpin = Filter_hairpin(Filter_CG_end)
    print ("After hairpin filter:", len(Filtered_with_hairpin))
    Filter_rep = Filter_check_repetition(Filtered_with_hairpin, max_allowed_aa_repetition)
    print ("After check of repetitions:",len(Filter_rep))    
    Filter_5GC_end = Filter_CG_5base_end(Filter_rep)    
    print ("After max 5 GC at end:",len(Filter_5GC_end))
    #Filter_5GC_end = Filter_rep 
    hamming_filter_forward = Filter_hamming_codTemp(Possible_primers_forward, Filter_5GC_end, 1, annealing_temp)
    print ("Nr of Forward primers after forward filter:",len(hamming_filter_forward))
    hamming_filter_reverse = Filter_hamming_codTemp(Possible_primers_reverse, hamming_filter_forward, 0, annealing_temp)
    print ("Nr of Forward primers after reverse filter:",len(hamming_filter_reverse))
    #save_to_file(hamming_filter_reverse,specie+"_"+Forward_or_Reverse+"_primers.txt")
    
    return hamming_filter_reverse, dictionary_all_primers
    '''
    if Forward_or_Reverse == "Forward":
        return hamming_filter_reverse, dictionary_all_primers
    else:
        return (Do_reverse(hamming_filter_reverse))
    '''

def Index_primers(dict_all_primers, list_of_primers, forward, length_sequence, length_primer):
    """
    This function will index each primer with the location
    """
    result_dict = {}
    if forward == 1:
        for i in list_of_primers:                    
            result_dict[(i)] = dict_all_primers[i]
    elif forward == 0:
        for i in list_of_primers:
            result_dict[(i)] = length_sequence-dict_all_primers[i]-length_primer+2
    return result_dict
    
    #result_dict = {}
    #for i in hamming_filter_reverse:
        #result_dict[Do_complement(i)] = dictionary_all_primers[i]
    #save_dic_to_file(result_dict, "test_"+Forward_or_Reverse+".txt")
    #print (result_dict)
    #return(result_dict)
    
def Pair_dic_primers(name_of_file,Forward, Reverse, lenght, Pair_min = 300, Pair_max = 1500):
    """
    This function will pair the forward and the reverse primer and save into a file as a result
    """
    output = open(name_of_file+'_result.txt', 'a')
    for i in Forward:        
        for j in Reverse:
            if (Reverse[j]-(Forward[i]+lenght) > Pair_min) and (Reverse[j]-(Forward[i]+lenght) < Pair_max):
                result = str(str(i)+" "+str(Forward[i])+"--->"+str(j)+" "+str(Reverse[j])+"\n")
                output.write(result)

def Run_the_primer_program():
    pass

if __name__ == "__main__":
    """
    Here is where the program will run
    """
    time_start = datetime.now()    
    
    #Needed variables
    #name_of_file = "NC_000866-Enterobact-phage-T4-complete-genome.fasta"
    name_of_file = "NC_001895-Enterobact-phage-P2-complete-genome.fasta"
    #name_of_file = "ecoli.fasta"    
    temp = 60
    lenght = 20
    GC_max = 60
    GC_max = GC_max*0.01
    GC_min = 40
    GC_min = GC_min*0.01
    Pair_max = 1500
    Pair_min = 300
    max_allowed_aa_repetition = 4
    annealing_temp = 15
    
    #Run the program
    Forward_Sequence, length_sequence = Open_file(name_of_file, "fasta")    
    Result_Forward, dictionary_all_primers_forward = Do_all_filters(Forward_Sequence,name_of_file, "Forward", lenght, temp, length_sequence, GC_max, GC_min, Pair_max, Pair_min, max_allowed_aa_repetition, annealing_temp)
    Result_index_Forward = Index_primers(dictionary_all_primers_forward, Result_Forward, 1, length_sequence, lenght)
    
    Reverse_Sequence = Do_complement(Forward_Sequence[::-1])
    Result_Reverse, dictionary_all_primers_reverse = Do_all_filters(Reverse_Sequence,name_of_file, "Reverse", lenght, temp, length_sequence, GC_max, GC_min, Pair_max, Pair_min, max_allowed_aa_repetition, annealing_temp)
    Result_index_Reverse = Index_primers(dictionary_all_primers_reverse, Result_Reverse, 0, length_sequence, lenght)
    
    Pair_dic_primers(name_of_file,Result_index_Forward, Result_index_Reverse, lenght)
    #print (Result_index_Reverse)

    
       
    time_end = datetime.now()
    time_to_run = time_end - time_start
    print ("Time to run the program:",time_to_run)