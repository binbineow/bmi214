#!/usr/bin/env python
#BMI214 Project 3
#Binbin Chen
#bchen45

#Function

#For our purposes, you should only find the non-bonded pairs once - at the
# beginning of the simulation using the
# input files - and use these as a constant set throughout the simulation.

#For each pair of interacting atoms, 
#you only need to compute the bond or non-bond energy ONCE per iteration.

# Summary of Steps for Coding the Simulation
# 
# Parse an input file containing atom, coordinate, initial velocity, connectivity, and parameter information.
# Iterate through 1000 time steps, and for each time step do the following:
# Update the velocities (for t+dt/2) on each atom (the forces should be initialized to zero
#                                                   when you do this in the first time step, 
#                                                   so your velocities won't change yet)


#import necessary packages
import numpy
from scipy.spatial.distance import squareform, pdist
import argparse
from collections import defaultdict

#input: opened input file
#output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery 
#round a float values ensuring it has two decimial places
def magic_round(x,n):
    return format(x,'.'+str(n)+'f')

#function take a list and covnert into a line
#input: opened input file
#output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery 
def give_line(list0,sep0,n=-1):
    if len(list0) >0:
        if n > -1:
            element0 = magic_round(list0[0], n)
        else:
            element0 = str(list0[0])
        line0 = str(element0)
        for element0 in list0[1:]:
            if n > -1:
                element0 = magic_round(element0, n)
            line0 = line0 + sep0 + str(element0)
        return line0    
    else:
        return ''
    
#function: parse the inputline specified and set the paramters into the default values
#input: opened input file
#output eight values specified below  
# Update the positions of each atom
# If you are keeping global variables for forces and energy then make sure that they are reset to zero.
# For each interaction pair (bonded and non-bonded):
# Calculate the potential energy and update the total potential energy 
def parse_input0():
    program_name = 'MyMD'
    parser0 = argparse.ArgumentParser(prog=program_name)
    #I did not write a loop here because it does not save me time
    #set parser
    parser0.add_argument('--iF',default='')
    parser0.add_argument('--kB',default=40000.0,type=float)
    parser0.add_argument('--kN',default=400.0,type=float)
    parser0.add_argument('--nbCutoff',default=0.50,type=float)
    parser0.add_argument('--m',default=12.0,type=float)
    parser0.add_argument('--dt',default=0.001,type=float)
    parser0.add_argument('--n',default=1000,type=int)
    parser0.add_argument('--out',default='')
    #parse the values from the input line
    input_values = parser0.parse_args()
    #set output file name correctly
    if input_values.out == '':
        input_values.out = input_values.iF.split('.')[0]
    #return
    return input_values


#function: the main class of this program
#input: read in the input file with values
#output: the kinetic energy in the second column, the bond potential energy in the third column, 
#     the nonbond potential energy in the fourth column, and the total energy in the last column.
#      All values should be written out to exactly 1 decimal place.  
class Atoms_info:
    def __init__(self,Input_values):
        self.pos_matrix = []
        self.vel_matrix = []
        self.dist_matrix = []
        self.dict_bond = defaultdict(list)
        self.dict_nonbond = defaultdict(list)
        self.input_values = Input_values
        ###get input position and distance
        self.get_pos0()
        self.dist0_matrix = self.get_dist(self.pos_matrix)
        self.dist_matrix = self.get_dist(self.pos_matrix)
        self.bond_pair = self.get_pair(self.dict_bond)
        self.get_nonbond()
        self.nonbond_pair = self.get_pair(self.dict_nonbond)
        #set accerlations
        self.a = numpy.zeros((len(self.pos_matrix),3))
        #prepare output
        
        
        
    #function   energy should be only calculated when key0 < x
                #getting the list of pairs
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery  
    def prepare_output(self):
        self.file_energy = open(self.input_values.out + '_out.erg','w+')
        self.file_energy.write(give_line(['#','step','E_k','E_b','E_nB','E_tot'],'\t')+'\n')
        self.file_rvc =  open(self.input_values.out + '_out.rvc','w+')
    
    #function   energy should be only calculat
                #getting the list of pairs
    #input: a list of energy
    #output:
    #he RVC output file will look very similar to the input RVC file, but is tab-delimited. 
    #During simulation, you should output all the atoms' IDs (1-indexed), positions, velocities and
    #connectivities every 10 time steps, including at t=0, concatenated together into one RVC output y 
    def write_energy(self,list0):
        self.file_energy.write(str(list0[0])+'\t'+give_line(list0[1:],'\t',n=1)+'\n')
    
    #function   energy should be only calculated when key0 < x
                #getting the list of pairs
    #input: the first line of rvc output and lines
    #output:
    #10 steps, after 20 steps, etc.)  You should include the original input RVC file and you should use 
    #1000 as the default number of total steps (feel free to perform longer simulations if you have the time and 

    def write_rvc(self,line0):
        self.file_rvc.write(line0)
        for i in range(0,len(self.pos_matrix)):
            line0_out = str(i+1)+'\t' + give_line(numpy.concatenate([self.pos_matrix[i],self.vel_matrix[i]]),'\t',n=4) 
            line0_out = line0_out + '\t'+ give_line([int(x)+1 for x in self.dict_bond[i]],'\t',n=0)
            self.file_rvc.write(line0_out+'\n')
    
    #function   print the warning when the energy overflows
    #input: iteration number when the program stops
    #output:     (inclination), so you should have about 100 frames in your output file. Positions and velocities should be
    #     written out to exactly 4 decimal places.
    def write_overflow(self,iteration0):
        self.file_energy.write('#The system becomes unstable at the iteration '+str(iteration0)+'. The program terminates')
        self.file_rvc.write('#The system becomes unstable at the iteration '+str(iteration0)+'. The program terminates')
        
    #function  close all functions
    #input: N/A
    #output: N/A  
    def close_output(self):
        self.file_energy.close()
        self.file_rvc.close()
    
    #function   energy should be only calculated when key0 < x
                #getting the list of pairs
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery    
    def get_pos0(self):
        n_atom = 0
        for line0 in open(self.input_values.iF,'rU'):
            if not '#' in line0:
                line0 = line0.rstrip().split()
                pos_info = [float(x) for x in line0[1:4]]
                vel_info = [float (x) for x in line0[4:7]]
                bond_info = [int(x)-1 for x in line0[7:]]
                #apppend
                self.pos_matrix.append(pos_info)
                self.vel_matrix.append(vel_info)
                self.dict_bond[n_atom] = bond_info
                #update atom number
                n_atom += 1
            else:
                self.first_line = line0
        #convert the matrix to numpy array
        self.pos_matrix = numpy.array(self.pos_matrix)
        self.vel_matrix = numpy.array(self.vel_matrix)
        
    #function   energy should be only calculated when key0 < x
                #getting the list of pairs
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery  
    def get_dist(self,pos_matrix):
        dist_matrix = squareform(pdist(pos_matrix,'euclidean'))
        #pickle.dump(dist_matrix,open('matrix.pic','w+'))
        return numpy.array(dist_matrix)
 
    #function   energy should be only calculated when key0 < x
                #getting the list of pairs
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery
    def get_pair(self,dict0):
        list_x = []
        list_y = []
        for key0, val0 in dict0.iteritems():
            for x in val0:
                if key0 < x:
                    list_x.append(key0)
                    list_y.append(x)
        return numpy.array([list_x,list_y])
    
    #functiong get nonbond pairs
    #input: self
    #output:  
    # The ERG file should contain your energy values at every 10th time step.  
    # The lines should be tab-delimited, with the time step number in the first column, 
    def get_nonbond(self):
        for i in range(0,len(self.pos_matrix)):
            list_nonbond = []
            for j in range(0,len(self.pos_matrix)):
                if self.dist0_matrix[i,j] < self.input_values.nbCutoff and i < j and not j in self.dict_bond[i]:
                    list_nonbond.append(j)
            self.dict_nonbond[i] = list_nonbond
     
    #function get sequence a and b
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery  
    def cal_kinetic(self):
        kinetic_energy = 0.5*self.input_values.m*numpy.square(self.vel_matrix)
        kinetic_energy = numpy.sum(kinetic_energy)
        return kinetic_energy

    #function convert a list of atom pairs into nergey 
    #input: atom pairs
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery  
    def cal_potential(self,list_to_cal,constant0):
        #print(self.dist_matrix.shape)
        delta_matrix = self.dist_matrix[list_to_cal[0],list_to_cal[1]]-self.dist0_matrix[list_to_cal[0],list_to_cal[1]]
        potential_energy = 0.5*constant0*numpy.square(delta_matrix)
        #print(potential_energy.shape)
        #we are double counting the energy, so /2
        return numpy.sum(potential_energy)

    #function calculate the energy 
    #input: input file value specified in self
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery  
    def cal_energy(self):
        kinetic_energy = self.cal_kinetic()
        bond_pot_energy = self.cal_potential(self.bond_pair,self.input_values.kB)
        nonbond_pot_energy = self.cal_potential(self.nonbond_pair,self.input_values.kN)
        total_energery = sum([kinetic_energy,bond_pot_energy,nonbond_pot_energy])
        return [kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery]
    
 
    #function calculate the force 
    #input: input file value specified in self
    #output: matrixes of force 
    def update_force(self):
        #bonded force
        for atom_i,atom_j in zip(self.bond_pair[0],self.bond_pair[1]):
            b = self.dist_matrix[atom_i,atom_j]
            force0 = self.input_values.kB * (b - self.dist0_matrix[atom_i,atom_j])
            force_vec = force0/b*(self.pos_matrix[atom_j]-self.pos_matrix[atom_i])
            #print(len(force_vec))
            self.force[atom_i] += force_vec
            self.force[atom_j] -= force_vec
        #unbonded force
        for atom_i,atom_j in zip(self.nonbond_pair[0],self.nonbond_pair[1]):
            b = self.dist_matrix[atom_i,atom_j]
            force0 = self.input_values.kN * (b - self.dist0_matrix[atom_i,atom_j])
            force_vec = force0/b*(self.pos_matrix[atom_j]-self.pos_matrix[atom_i])
            self.force[atom_i] += force_vec
            self.force[atom_j] -= force_vec

    
    #function main update function
    #input: force and position
    #output: new a, force and postion
    # connectivities every 10 time steps, including at t=0, concatenated together into one RVC output file,
    # like the example. (In other words, you should write to the output before any steps have been taken, after 
    def update_velocity_verlet(self):
        #clear up force
        self.force = numpy.zeros((len(self.pos_matrix),3))
        self.vel_matrix += 0.5*self.a*self.input_values.dt
        self.pos_matrix = self.pos_matrix + self.vel_matrix*self.input_values.dt
        self.dist_matrix = self.get_dist(self.pos_matrix)
        self.update_force()
        self.a = 1.0/self.input_values.m*self.force
        self.vel_matrix +=  0.5*self.a*self.input_values.dt

    


#function: the main function of this project
#input: read in the input file with k value
#output: generate energy and frame files
def main():
    #initinize files and get I/O file names
    # input values
    #n is the number of iteration
    input_values = parse_input0()
    #set initial conditions
    protein = Atoms_info(input_values)
    protein.prepare_output()
    protein.write_rvc(protein.first_line)
    #main loop
    #in this example, atom index starts at 0 except during the ouput 
    #remember to do i+1 during output
    total_energy0 = protein.cal_energy()[-1]
    for iteration0 in range(1,input_values.n+1):
        #list1 = protein.cal_energy()
        protein.update_velocity_verlet()
        #calculate energy
        list1 = protein.cal_energy()
        if list1[-1] > 10 * total_energy0 or list1[-1] < total_energy0/10.0:
                protein.write_overflow(iteration0)
                break
        if numpy.mod(iteration0,10) == 0:
            protein.write_rvc('#At time step '+str(iteration0)+',energy = '+str(list1[-1])+'kJ\n')
            protein.write_energy([iteration0]+list1)
    #close the output file
    protein.close_output()
        


#run the main function
if __name__ == '__main__':
    main()        
        
        
        