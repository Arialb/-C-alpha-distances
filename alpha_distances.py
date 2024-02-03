import sys
import numpy as np
from scipy.stats import ttest_ind, ttest_rel

# Function parsing the PDB file
def parse_pdb(filename, chain):
  pdb_coords={}
  f=open(filename,'r')
  for line in f:
      line=line.rstrip()
      # Check the field ATOM at the beginning of the line
      # Check the correct chain in column 22
      if line[:4]!='ATOM' or line[21]!=chain: continue
      resn=int(line[22:26].strip())
      atom=line[12:16].strip()
      x=float(line[30:38].strip())
      y=float(line[38:46].strip())
      z=float(line[46:54].strip())
      coord=[x,y,z]
      # Initialize the dictionary with atoms coordinates of residue resn
      pdb_coords[resn]=pdb_coords.get(resn,{})
      # Add the atom's coordinates
      pdb_coords[resn][atom]=coord
  return pdb_coords


# Function calculating the distance between to points
def get_distance(coord1,coord2):
  return np.sqrt((coord1[0]-coord2[0])**2+\
                 (coord1[1]-coord2[1])**2+\
                 (coord1[2]-coord2[2])**2)


# Function returning the diatnces between all consecutives CA
def get_ca_dist(pdb_coords):
  ca_dists=[]
  keys=list(pdb_coords.keys())
  keys.sort()
  n=len(keys)
  for i in range(n-1):
    if int(keys[i+1])-int(keys[i])==1:
      if "CA" in pdb_coords[keys[i]] and "CA" in pdb_coords[keys[i+1]]:
        dist=get_distance(pdb_coords[keys[i]]['CA'],pdb_coords[keys[i+1]]['CA'])
        ca_dists.append(dist)
  return ca_dists

def count_CA(filename):
  f=open(filename,'r')
  atoms=0
  CA_atoms=0
  for line in f:
    if line[:4]=='ATOM':
      atoms=atoms+1
      if "CA" in line:
        CA_atoms=CA_atoms+1
  return atoms, CA_atoms


import os
# get a list of all file in a path
apath = "./structures/b-pos/"
alpha_files = os.listdir( apath )
bpath = "./structures/c-pos/"
beta_files = os.listdir( bpath )


ca_a_dists=[]
ca_b_dists=[]
chaina=["P", 'E', 'E','E','P','E'] #chains for b-state
chainb = ['E','P','R','P','P','P'] #chains for c-state
for a in range(len(alpha_files)):
  pdb_coords=parse_pdb(apath+alpha_files[a],chaina[a]) #obtain a dictionary pdb_coords[resn][atom]=coord
  ca_a_dists=ca_a_dists+get_ca_dist(pdb_coords) #obtain a list of distances added to the previous ones
for b in range(len(beta_files)):
  pdb_coords=parse_pdb(bpath+beta_files[b],chainb[b]) #obtain a dictionary pdb_coords[resn][atom]=coord
  ca_b_dists=ca_b_dists+get_ca_dist(pdb_coords) #obtain a list of distances added to the previous ones

  
print (len(ca_a_dists))
print ('Dist CA mean in b-pos:',np.mean(ca_a_dists),"std dev:",np.std(ca_a_dists))
print (len(ca_b_dists))
print ('Dist CA mean in c-pos:',np.mean(ca_b_dists),"std dev:",np.std(ca_b_dists))
print(ttest_ind(ca_a_dists, ca_b_dists))