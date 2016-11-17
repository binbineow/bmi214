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
# Update the positions of each atom
# If you are keeping global variables for forces and energy then make sure that they are reset to zero.
# For each interaction pair (bonded and non-bonded):
# Calculate the potential energy and update the total potential energy
# Calculate the force and update the total forces in each dimension for each atom
# Update the velocities (for t+dt) on each atom and calculate the kinetic energy
# Write the appropriate output every 10 time steps. Do this for all file types: 
# (rvc, crd, files needed to generate euc)


#import necessary packages
import numpy
from scipy.spatial.distance import squareform, pdist
import argparse
from collections import defaultdict
import cPickle as pickle






#function: ouptut performance to knn.out
#input: assignment of each points to clusters
#output: write values into the kmeans.out
def write_output(file_name0,list0):
    fileout = open(file_name0,'w+')
    index0 = 1
    for cluster0 in list0:
        fileout.write(str(index0)+'\t'+str(cluster0+1)+'\n')
        index0 +=1

#function: parse the inputline specified and set the paramters into the default values
#input: opened input file
#output seqA and seqB    
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


#function: the main function of this project
#input: read in the input file with k value
#output: clusters of points based on centroid provided or generated       
class Atoms_info:
    def __init__(self,Input_values):
        self.pos_matrix = []
        self.vel_matrix = []
        self.dist_matrix = []
        self.dict_bond = defaultdict(list)
        self.dict_nonbond = defaultdict(list)
        self.input_values = Input_values
        ###get input
        self.get_pos0()
        self.dist0_matrix = self.get_dist(self.pos_matrix)
        self.dist_matrix = self.get_dist(self.pos_matrix)
        ##calculate distance matrix at step0
        #self.dist_matrix_xyz0 = []
        #self.dist_matrix_xyz0.append(self.get_dist(self.pos_matrix[:,0:1]))
        #print(self.get_dist(self.pos_matrix[:,0:1]).shape)
        #self.dist_matrix_xyz0.append(self.get_dist(self.pos_matrix[:,1:2]))
        #self.dist_matrix_xyz0.append(self.get_dist(self.pos_matrix[:,2:3]))
        # bond_pair is a list of tuble list convenient for calculating 
        self.bond_pair = self.get_pair(self.dict_bond)
        self.get_nonbond()
        self.nonbond_pair = self.get_pair(self.dict_nonbond)
        self.a = numpy.zeros((len(self.pos_matrix),3))

     
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
    
    #function   energy should be only calculated when key0 < x
                #getting the list of pairs
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery
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

    #function get sequence a and b
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery  
    def cal_potential(self,list_to_cal,constant0):
        #print(self.dist_matrix.shape)
        delta_matrix = self.dist_matrix[list_to_cal[0],list_to_cal[1]]-self.dist0_matrix[list_to_cal[0],list_to_cal[1]]
        potential_energy = 0.5*constant0*numpy.square(delta_matrix)
        #print(potential_energy.shape)
        #we are double counting the energy, so /2
        return numpy.sum(potential_energy)

    #function get sequence a and b
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery  
    def cal_energy(self):
        kinetic_energy = self.cal_kinetic()
        bond_pot_energy = self.cal_potential(self.bond_pair,self.input_values.kB)
        nonbond_pot_energy = self.cal_potential(self.nonbond_pair,self.input_values.kN)
        total_energery = sum([kinetic_energy,bond_pot_energy,nonbond_pot_energy])
        return kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery
    
    #function get sequence a and b
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery 
    def update_force_old(self):
        #calculate distance based on x,y,z axises 
        self.dist_matrix_xyz = []
        self.dist_matrix_xyz.append(self.get_dist(self.pos_matrix[:,0:1]))
        self.dist_matrix_xyz.append(self.get_dist(self.pos_matrix[:,1:2]))
        self.dist_matrix_xyz.append(self.get_dist(self.pos_matrix[:,2:3]))
        #bonded force
        for atom_i,atom_j in zip(self.bond_pair[0],self.bond_pair[1]):
            for dimention0 in range(0,3):
                delta_r = self.dist_matrix_xyz[dimention0][atom_i,atom_j] - self.dist_matrix_xyz0[dimention0][atom_i,atom_j]
                vec_j2i = self.pos_matrix[atom_j,dimention0] - self.pos_matrix[atom_i,dimention0]
                self.force[atom_i,dimention0] += self.input_values.kB*delta_r*vec_j2i/abs(vec_j2i)
        #unbonded force
        for atom_i,atom_j in zip(self.nonbond_pair[0],self.nonbond_pair[1]):
            for dimention0 in range(0,3):
                delta_r = self.dist_matrix_xyz[dimention0][atom_i,atom_j] - self.dist_matrix_xyz0[dimention0][atom_i,atom_j]
                vec_j2i = self.pos_matrix[atom_j,dimention0] - self.pos_matrix[atom_i,dimention0]
                #print(delta_r,vec_j2i)
                self.force[atom_i,dimention0] += self.input_values.kN*delta_r*vec_j2i/abs(vec_j2i)
        #function get sequence a and b
    
    #function get sequence a and b    
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery 
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

    
    #function get sequence a and b    
    #input: opened input file
    #output: a list of kinetic_energy,bond_pot_energy,nonbond_pot_energy,total_energery 
    def upate_velocity_verlet(self):
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
#output: clusters of points based on centroid provided or generated
def main():
    #initinize files and get I/O file names
    # input values
    #n is the number of iteration
    input_values = parse_input0()
    #set initial conditions
    protein = Atoms_info(input_values)
    #main loop
    #in this example, atom index starts at 0 except during the ouput 
    #remember to do i+1 during output
    list1 = protein.cal_energy()
    print(list1)
    for iteration0 in range(1,input_values.n+1):
        #list1 = protein.cal_energy()
        protein.upate_velocity_verlet()
        if numpy.mod(iteration0,10) == 0:
            list1 = protein.cal_energy()
            print(iteration0,list1)
        
    #output
#     The ERG file should contain your energy values at every 10th time step.  
#     The lines should be tab-delimited, with the time step number in the first column, 
#     the kinetic energy in the second column, the bond potential energy in the third column, 
#     the nonbond potential energy in the fourth column, and the total energy in the last column.
#      All values should be written out to exactly 1 decimal place.

#     The RVC output file will look very similar to the input RVC file, but is tab-delimited. 
#     During simulation, you should output all the atoms' IDs (1-indexed), positions, velocities and
#      connectivities every 10 time steps, including at t=0, concatenated together into one RVC output file,
#       like the example. (In other words, you should write to the output before any steps have been taken, after 
#     10 steps, after 20 steps, etc.)  You should include the original input RVC file and you should use 
#     1000 as the default number of total steps (feel free to perform longer simulations if you have the time and 
#     inclination), so you should have about 100 frames in your output file. Positions and velocities should be
#      written out to exactly 4 decimal places.
    # #At time step 610,energy = 273.911kJ is an example of in-between frame
    




if __name__ == '__main__':
    main()        
        
        
        