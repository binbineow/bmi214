#!/usr/bin/env python
#BMI214 Project 1
#Binbin Chen
#bchen45

#the sum of all number in the numberList
#Index specificaiton
#the first sequence will be represented by the columns in the match matrix (i)
#the second sequence will be represented by the columns in the match matrix (j)
#three difference betweeen local and global alignments
#local never records a score below 0, so turn all negative scores to 0
#local looks for the best s ore anywhere in the matrix
#local requries negative scores in the matrix

##Ends-free alignment
#we exclude all and trailing insertions or deletions, basically once i 
#or j hits 0, stop the tracing (output exclude trailing or leading) 
#and not penalize the score

#import necessary packages
import numpy as np
import sys
from collections import defaultdict

#check if two numbers are the same
def is_same(val1,val2):
    delta0 = 1e-7
    if abs(val1-val2) < delta0:
        return True
    else:
        return False

#function: check values of a list, if any of them is negative, replace it with zero
#input: a list
#output: a list with no negative scores
def give_zero(list0):
    list_out = []
    for x in list0:
        if x<0:
            list_out.append(0)
        else:
            list_out.append(x)
    return list_out


#define input specification class
class input0():
    '''class to store input specficiation
        seqA: the first sequence to be aligned (A)
        seqB: the second sequence to be aligned (B)
        is_local: whether the user resires a local (True) or global alignment (False)
        gap_open_A: gap opening penalty for A
        gap_extention_A: gap extension penalty for A
        gap_open_B: gap opening penalty for B
        gap_extention_B: gap extension penalty for B
        symbol_A: a string of all symbols used in A in order
        symbol_B: a string of all symbols used in B in order
        symbol_length: how many different symbols are in the A or B (should be the same for A&B)
        match_dict: a dictioanry containing match matrix information 
                    format: a length two string symbol from A + symbol from B -> match score
    '''
    #intitaize the input paramteres 
    def __init__(self):
        self.seqA = ''
        self.seqB = ''
        self.is_local = True
        self.gap_open_A = 0
        self.gap_extention_A = 0
        self.gap_open_B = 0
        self.gap_extention_B = 0
        self.symbol_A = ''
        self.symbol_B = ''
        self.symbol_length = 0
        self.match_dict = dict()
    
    #function get sequence a and b
    #input: opened input file
    #output seqA and seqB    
    def get_seq(self,input_file0):
        self.seqA = input_file0.readline().rstrip()
        self.seqB = input_file0.readline().rstrip()
    
    #function get sequence a and b similar
    #input: opened input file
    #output seqA and seqB 
    def get_local_or_global(self,input_file0):
        self.is_local = input_file0.readline().rstrip() == '1'
        
    #function get sequence a and b similar
    #input: opened input file
    #output seqA and seqB 
    def get_gap_penalty(self,input_file0):
        penalty_values = input_file0.readline().rstrip().split(' ')
        self.gap_open_A = float(penalty_values[0])
        self.gap_extention_A = float(penalty_values[1])
        self.gap_open_B = float(penalty_values[2])
        self.gap_extention_B = float(penalty_values[3])
        
    #function get sequence a and b similar
    #input: opened input file
    #output seqA and seqB 
    def get_symbols(self,input_file0):
        #length of the symbol length is not necessary 
        #read in A
        _ = input_file0.readline() 
        self.symbol_A = input_file0.readline().rstrip()
        #length of the symbol length is not necessary 
        #read in B
        _ = input_file0.readline()
        self.symbol_B = input_file0.readline().rstrip()
        if len(self.symbol_A) == len(self.symbol_B):
            self.symbol_length = len(self.symbol_A)
        else:
            #in case things went unexpectedly
            print('I <3 Emily.')
                  
    #function get sequence a and b similar
    #input: opened input file
    #output seqA and seqB
    #this function can be only called after get_symbols is called
    #assuming no information is after match matrix
    def create_match_dict(self,input_file0):
        for line0 in input_file0:
            if len(line0) > 1:
                line0 = line0.rstrip().split(' ')
                self.match_dict[line0[2]+line0[3]] = float(line0[4])
    
    #this method returns the match score given i and j
    def get_s(self,i,j):
        #print(i,j)
        symA = self.seqA[i-1]
        symB = self.seqB[j-1]
        #print(str(symA)+str(symB)+str(self.match_dict[symA+symB]))
        return self.match_dict[symA+symB]

                


#function: the main function of this project similar
#input: read in the input file with specification and sequence information
#output: perform Smith-Waterman alignment (global or local), output best matching scores
#        and all possible best alignments
class map0():
    '''class to store alignment maps including both dynamic programming and paths'''
    def __init__(self,seqA,seqB):
        #i, length of A
        len_seqA = len(seqA)
        #j, length of B
        len_seqB = len(seqB)
        #set the first number i representing seqA
        self.matrix_s = np.zeros([len_seqA+1,len_seqB+1])
        self.matrix_M = np.zeros([len_seqA+1,len_seqB+1])
        #dynamic program values when we add a gap to B
        self.matrix_Ix = np.zeros([len_seqA+1,len_seqB+1])
        #dynamic program values when we add a gap to A
        self.matrix_Iy = np.zeros([len_seqA+1,len_seqB+1])
        #a default dictionary to store the path 
        #three parts joint by _
        #first  M, Ix, Iy + 
        #second i
        #third j
        #point to tuble (-1,-1) M, (i-1,j) Ix, (i,j-1) Iy
        self.path_dict = defaultdict(list)
    
    #based on if it's a local or global alignment
    #get the appropriate best score and starting M
    #because we ignore trailing sequence, so the sequence will never end at Ix and Iy
    def output_best_socre(self,is_local,lenA,lenB):
        if not is_local:
            #to dealing with trailing situation, we should also include the possibility 
            #of adding gap at the end the shorter sequence
            best_score = 0
            #assuming the sequence ends at the sequence B
            for i in range(1,lenA+1):
                if self.matrix_M[i,lenB] > best_score:
                    best_score = self.matrix_M[i,lenB]
                    best_i = i
                    best_j = lenB
            #assuming the sequence ends at the sequence A
            for j in range(1,lenB+1):
                if self.matrix_M[lenA,j] > best_score:
                    best_score = self.matrix_M[lenA,j]
                    best_i = lenA
                    best_j = j                
        else:
            #at the local setting, looking for the best score
            best_score = 0
            for i in range(1,lenA+1):
                for j in range(1,lenB+1):
                    if self.matrix_M[i,j]> best_score:
                        best_score = self.matrix_M[i,j]
                        best_i = i
                        best_j = j
                    
        return best_score, best_i, best_j
    
    def output_path(self,project_input,type0,i,j,seqA0,seqB0,list0):  
        #type0 indicates the current matching stratergy M, Ix, Iy
        #boundary condition
        if is_same(i,1) or is_same(j,1) or (is_same(getattr(self,'matrix_'+type0)[i,j],0) and project_input.is_local):
            if not type0 == 'M':
                list0.append([seqA0,seqB0])
            else:
                adding_A = project_input.seqA[i-1]
                adding_B = project_input.seqB[j-1]
                list0.append([adding_A+seqA0,adding_B+seqB0])
        #trace back based on path_dict
        else:
            #get new characters for seqA and seqB based on type0
            if type0 == 'M':
                adding_A = project_input.seqA[i-1]
                adding_B = project_input.seqB[j-1]
            elif type0 == 'Ix':
                adding_A = project_input.seqA[i-1]
                adding_B = '_' 
            elif type0 == 'Iy':
                adding_A = '_'
                adding_B = project_input.seqB[j-1]            
            index0 = '_'.join([type0,str(i),str(j)])
            #checking if the path info is in the path dictionry, mainly for debugging purpose
            if not index0 in self.path_dict:
                assert(False)
            for path0 in self.path_dict[index0]:
                #based on what path0 is, trace back recursively
                if path0 == 'M':
                    self.output_path(project_input,'M',i-1,j-1,adding_A+seqA0,adding_B+seqB0,list0)
                if path0 == 'Ix':
                    self.output_path(project_input,'Ix',i-1,j,adding_A+seqA0,adding_B+seqB0,list0)
                if path0 == 'Iy':
                    self.output_path(project_input,'Iy',i,j-1,adding_A+seqA0,adding_B+seqB0,list0)


def get_best_map(sw_map,project_input,lenA,lenB):
    #the M, Ix, Iy have been initiated 
    #compute the matrix
    for i in range(1,lenA+1):
        for j in range(1,lenB+1):
            #convert i and j into string for storage purpose
            str_i = str(i)
            str_j = str(j)
            #compute M matrix based on Smith_watermen algorithm
            val_M_1 = sw_map.matrix_M[i-1,j-1]+project_input.get_s(i,j)
            val_M_2 = sw_map.matrix_Ix[i-1,j-1]+project_input.get_s(i,j)
            val_M_3 = sw_map.matrix_Iy[i-1,j-1]+project_input.get_s(i,j)
            #pick the best
            max0 = max(val_M_1,val_M_2,val_M_3)
            sw_map.matrix_M[i,j] = max0
            if i == lenA and j == lenB:
                print ('here')
                print max0
                print(sw_map.matrix_M[i,j])
            #store the possible path for each M matrix element
            if is_same(max0,val_M_1):
                sw_map.path_dict['_'.join(['M',str_i,str_j])].append('M')
            if is_same(max0,val_M_2):
                sw_map.path_dict['_'.join(['M',str_i,str_j])].append('Ix') 
            if is_same(max0,val_M_3):
                sw_map.path_dict['_'.join(['M',str_i,str_j])].append('Iy')   
            #compute Ix matrix based on Smith_watermen algorithm
            #Ix: adding _ into B, so use penalty scores from B
            #print(i)
            #print(j)
            #print(sw_map.matrix_M[i-1,j])
            #print(project_input.gap_open_B)
            val_M_1 = sw_map.matrix_M[i-1,j] - project_input.gap_open_B
            val_M_2 = sw_map.matrix_Ix[i-1,j] - project_input.gap_extention_B
            #pick the best
            max0 = max(val_M_1,val_M_2)
            sw_map.matrix_Ix[i,j] = max0
            #store the possible path for each Ix matrix element
            if is_same(max0,val_M_1):
                sw_map.path_dict['_'.join(['Ix',str_i,str_j])].append('M')
            if is_same(max0,val_M_2):
                sw_map.path_dict['_'.join(['Ix',str_i,str_j])].append('Ix') 
            #compute Iy matrix based on Smith_watermen algorithm
            #Ix: adding _ into A, so use penalty scores from A
            val_M_1 = sw_map.matrix_M[i,j-1] - project_input.gap_open_A
            val_M_2 = sw_map.matrix_Iy[i,j-1] - project_input.gap_extention_A
            #pick the best
            max0 = max(val_M_1,val_M_2)
            sw_map.matrix_Iy[i,j] = max0
            #store the possible path for each Ix matrix element
            if is_same(max0,val_M_1):
                sw_map.path_dict['_'.join(['Iy',str_i,str_j])].append('M')
            if is_same(max0,val_M_2):
                sw_map.path_dict['_'.join(['Iy',str_i,str_j])].append('Iy') 
    


#function: the main function of this project similar
#input: read in the input file with specification and sequence information
#output: perform Smith-Waterman alignment (global or local), output best matching scores
#        and all possible best alignments
def get_best(sw_map,project_input,type0,i,j):

    #if the recurrsion reaches the boundary condition
    if ('_'.join([type0,str(i),str(j)]) in sw_map.path_dict):
        max0 = getattr(sw_map,'matrix_'+type0)[i,j]
    elif is_same(i,1) or is_same(j,1):        
        #the boundary condition, ignore the leading match 
        if not type0 == 'M':
            max0 = 0
        else:
            max0 = project_input.get_s(i,j)
    else:
        if type0 == 'M':
            #without adding gap 
            val1 = get_best(sw_map,project_input,'M',i-1,j-1) + project_input.get_s(i,j)
            val2 = get_best(sw_map,project_input,'Ix',i-1,j-1) + project_input.get_s(i,j)
            val3 = get_best(sw_map,project_input,'Iy',i-1,j-1) + project_input.get_s(i,j)
            #if this is a local alignment 
            if project_input.is_local:
                [val1,val2,val3] = give_zero([val1,val2,val3])
            max0 = max(val1,val2,val3)
            sw_map.matrix_M[i,j] = max0
            #update path
            #depending on the best score
            #multiple path possible by using append
            if is_same(val1,max0):
                sw_map.path_dict['_'.join(['M',str(i),str(j)])].append('M')
            if is_same(val2,max0):
                sw_map.path_dict['_'.join(['M',str(i),str(j)])].append('Ix')    
            if is_same(val3,max0):
                sw_map.path_dict['_'.join(['M',str(i),str(j)])].append('Iy')
        elif type0 == 'Ix':
            #adding a gap into B
            val1 = get_best(sw_map,project_input,'M',i-1,j) - project_input.gap_open_B
            val2 = get_best(sw_map,project_input,'Ix',i-1,j) - project_input.gap_extention_B
            max0 = max(val1,val2)
            #if this is a local alignment 
            if project_input.is_local:
                [val1,val2] = give_zero([val1,val2])
            #update path
            sw_map.matrix_Ix[i,j] = max0
            if is_same(val1,max0):
                sw_map.path_dict['_'.join(['Ix',str(i),str(j)])].append('M')
                #print '_'.join(['Ix',str(i),str(j)])
                #print sw_map.path_dict['_'.join(['Ix',str(i),str(j)])]
            if is_same(val2,max0):
                sw_map.path_dict['_'.join(['Ix',str(i),str(j)])].append('Ix')
        else:
            #adding a gap into A
            val1 = get_best(sw_map,project_input,'M',i,j-1) - project_input.gap_open_A
            val2 = get_best(sw_map,project_input,'Iy',i,j-1) - project_input.gap_extention_A
            max0 = max(val1,val2)
            #if this is a local alignment 
            if project_input.is_local:
                [val1,val2] = give_zero([val1,val2])
            #update path
            sw_map.matrix_Iy[i,j] = max0
            if is_same(val1,max0):
                sw_map.path_dict['_'.join(['Iy',str(i),str(j)])].append('M')
            if is_same(val2,max0):
                sw_map.path_dict['_'.join(['Iy',str(i),str(j)])].append('Iy')      
    return max0       
        


#function: the main function of this project similar
#input: read in the input file with specification and sequence information
#output: perform Smith-Waterman alignment (global or local), output best matching scores
#        and all possible best alignments
def sw_alignment(project_input):
    #the Smith-Waterman map
    sw_map = map0(project_input.seqA,project_input.seqB)
    #recursively get the best alignment back proprogating from M(N,M)
    lenA = len(project_input.seqA)
    lenB = len(project_input.seqB)
    get_best_map(sw_map,project_input,lenA,lenB)
    return sw_map
    
    

#function: the main function of this project similar
#input: read in the input file with specification and sequence information
#output: perform Smith-Waterman alignment (global or local), output best matching scores
#        and all possible best alignments
#score for the best alignment (rounded to the first decimal place i.e. 3.1415->3.1).
def write_alignment_output(best_score,list_alignment,output_file_name):
    file_out = open(output_file_name,'w+')
    #rounding best score to a decimal 
    best_score = round(best_score,1)
    file_out.write(str(best_score)+'\n\n')
    for two_lines in list_alignment:
        for seq0 in two_lines:
            file_out.write(seq0+'\n')
        file_out.write('\n')
    file_out.close()

def print_path_dict(dict0,lenA,lenB):
    print('M map')
    for i0 in range(1,lenA+1):
        out0 = ''
        i = str(i0)
        for j0 in range(1,lenB+1):
            j = str(j0)
            if '_'.join(['M',i,j]) in dict0:
                out0 = out0+','+'|'.join(dict0['_'.join(['M',i,j])])
            else:
                out0 = out0+',0'
        print out0
    print('Ix map')
    for i0 in range(1,lenA+1):
        out0 = ''
        i = str(i0)
        for j0 in range(1,lenB+1):
            j = str(j0)
            #print('_'.join(['Ix',i,j]))
            if '_'.join(['Ix',i,j]) in dict0:
                out0 = out0+','+'|'.join(dict0['_'.join(['Ix',i,j])])
            else:
                out0 = out0+',0'
        print out0
    
    print('Iy map')
    for i0 in range(1,lenA+1):
        out0 = ''
        i = str(i0)
        for j0 in range(1,lenB+1):
            j = str(j0)
            if '_'.join(['Iy',i,j]) in dict0:
                out0 = out0+','+'|'.join(dict0['_'.join(['Iy',i,j])])
            else:
                out0 = out0+',0'
        print out0

#function: the main function of this project
#input: read in the input file with specification and sequence information
#output: perform Smith-Waterman alignment (global or local), output best matching scores
#        and all possible best alignments
def main():
    #initinize files and get I/O file names
    project_input = input0()
    #input_file_name = sys.argv[1]
    #output_file_name = sys.argv[2]
    #these two inputs are for local machine debugging only
    input_file_name=input_name
    output_file_name=output_name
    #open the input file
    input_file0 = open(input_file_name,'r')
    #get specific input parameters using input0 methods
    project_input.get_seq(input_file0)
    project_input.get_local_or_global(input_file0)
    project_input.get_gap_penalty(input_file0)
    project_input.get_symbols(input_file0)
    project_input.create_match_dict(input_file0)
    #the main alignment algorithm
    sw_map = sw_alignment(project_input)
    #get the size of the input sequence
    lenA = len(project_input.seqA)
    print(lenA)
    lenB = len(project_input.seqB)
    print(lenB)
    #get results including both the best score and all possible good alignments
    [best_score, best_i, best_j] = sw_map.output_best_socre(project_input.is_local,lenA,lenB)
    print('best')
    print(best_score,best_i,best_j)
    print_path_dict(sw_map.path_dict,lenA,lenB)
    print(sw_map.matrix_M)
    list_alignment = []
    sw_map.output_path(project_input,'M',best_i,best_j,'','',list_alignment)
    #write output into the file
    write_alignment_output(best_score,list_alignment,output_file_name)

#the following line is for debugging in my local machine only
path0 = '/Users/binbineow2/Documents/Medical_school/PhD/BMI214/p1-2/examples/'
input_name = path0+'alignment_example3.input'
output_name = path0+'alignment_example3.output_user'

if __name__ == '__main__':
    main()        
        
        
        