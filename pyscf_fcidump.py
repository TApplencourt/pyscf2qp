#!/usr/bin/env python

# Import
import numpy
from pyscf.pbc import gto, df, scf, tools

#  __                                                     
# /__ |  _  |_   _. |   ._   _. ._ _. ._ _   _ _|_  _  ._ 
# \_| | (_) |_) (_| |   |_) (_| | (_| | | | (/_ |_ (/_ |  
#                       |                                 
#
use_density_fiting = True #If unset, the RHF will have an hard time to converge.
int_threshold = 1E-8 # The integral will be not printed in they are bellow that

# ___                                          
#  |  ._  o _|_ o  _. | o _   _. _|_ o  _  ._  
# _|_ | | |  |_ | (_| | | /_ (_|  |_ | (_) | | 
#                                              

# Cell creation
# TODO: What is the difference between gto.M and gto.Cell ?!
L=1
cell = gto.M(
    atom = '''He 0. 0. 0.
	      He 0.45 0. 0.
	      He 0. 0.5 0.
	      He 0.5 0.5 0.''',
    basis = 'sto-3g',
    a =  numpy.diag([L,L,L]), # Cell dimension
    gs = [10]*3,
    verbose= 2)


natom = len(cell.atom_coords())
print 'n_atom',   natom
print 'num_elec', cell.nelectron
print 'nucl_num', len(cell.atom_coords())

# Only Restricted HF is avalaible in pysf
assert (cell.spin == 0) # cell.spin == 2*S aka  multiplicity-1.

# We handle only gamma point for now.
# WARNING: But pyscf have a strange behaviour in ERI for  gama point (The dimmesion is reduced by (1) )
# 	   So we shift the gamma point.

kpts_abs = [ [1E-5, 0., 0.] ]
kpts = numpy.array([cell.get_abs_kpts(k) for k in kpts_abs])
# More the general but with the bug
#nk = [1,1,1] # Number of Kpoint in each direction
#kpts = cell.make_kpts(nk) #kpts in absolute value (unit 1/Bohr).

print kpts

nkpt = len(kpts)
#
# /__  _  ._   _  ._ _. _|_ o  _  ._    |\/| / \
# \_| (/_ | | (/_ | (_|  |_ | (_) | |   |  | \_/
#
# SCF on all the kpts. 
# Note: This will not compute the Energy. It will just return you an engine
kmf = scf.KRHF(cell, kpts, exxdiv=None) #Note that on Restictred HF are avalaible.

# 2-electron integrals may not be accurate enough for
# all-electron calculation.  It's recommended to use MDF (mixed density fitting)
# technique to improve the accuracy.

use_density_fiting = True
if use_density_fiting:
	mydf = df.MDF(cell) # MDF seem beter than DF. No idea why. 
	mydf.kpts = kpts
	mydf.auxbasis = 'weigend' # Note: No idea what is this. But the use is all the time
	kmf.with_df = mydf

# Do thr SCF and the MO. 
e_scf = kmf.kernel()
print 'e_scf',e_scf

print ''
mo_coeff = kmf.mo_coeff # List of mo_coeff for each k-point
nmo = mo_coeff.shape[1]

print 'mo_tot_num',nmo

# Wrote all the parameter need to creat a dummy EZFIO folder who will containt the integral after.
# More an implentation detail than a real thing
with open('param','w') as f:
        f.write(' '.join(map(str,(cell.nelectron, nmo, natom))))
#                             _                             
# |\ |      _ |  _   _. ._   |_)  _  ._      |  _ o  _  ._  
# | \| |_| (_ | (/_ (_| |    | \ (/_ |_) |_| | _> | (_) | | 
#                                    |                      
e_nuc = kmf.energy_nuc() # We need it, because of the evalf sumation not implemtented in QP. 
print 'nucl_repul', e_nuc
with open('e_nuc','w') as f:
        f.write(str(e_nuc))

# ___                                              
#  |  ._ _|_  _   _  ._ _. |  _   |\/|  _  ._   _  
# _|_ | | |_ (/_ (_| | (_| | _>   |  | (_) | | (_) 
#                 _|                              
from itertools import product

def gen_mono_MO(mo_coeff,l_int,idx_kpt=0):
	# 2Id transfortion Transformation. For now we handle only one K point.
        l_int_mo = reduce(numpy.dot, (mo_coeff[idx_kpt].T, l_int[idx_kpt], mo_coeff[idx_kpt])) #This formula is only right for one kpt.
	for i,j in product(range(nmo), repeat=2):
                int_ = l_int_mo[i,j].real # Only to handle the gamma-point bug.
		if abs(int_) > int_threshold:
			yield (i+1,j+1, int_)

# A tuple with the name, and the value of the integrals for each k-pts
t_ao = ('kinetic', cell.pbc_intor('cint1e_kin_sph',kpts=kpts)) #  Note .intor function computes integrals without PBC
v_ao = ('nuclear', [kmf.with_df.get_nuc(kpt) for kpt in kpts]) # Ref:  `get_hcore` in pbc/scf file

#Print
for idx_kpt in range(nkpt): #Prepare the jon when we will have many k point
	for name, ao in (t_ao,v_ao):
		with open('%s_mo' % name,'w') as f:
			for mono in gen_mono_MO(mo_coeff,ao,idx_kpt):
				f.write('%s %s %s\n'% mono)

# ___                              _    
#  |  ._ _|_  _   _  ._ _. |  _   |_) o 
# _|_ | | |_ (/_ (_| | (_| | _>   |_) | 
#                 _|                    
#

def gen_eri_MO(kernel,kpts, qkpt_idx):
	     
	     # Compute the ERI
	     eri_kpt_flatish = kernel.with_df.ao2mo([kernel.mo_coeff[i] for i in qkpt_idx],
                                            	    [kpts[i] for i in qkpt_idx])

	     nmo = kernel.mo_coeff.shape[1]
	     # Putin in a good shape
             eri_kpt = eri_kpt_flatish.reshape((nmo,)*4)
	
	     # Warning: We need to handle the 8-fold symetrie. If not this will bug in the QP size.
	     # If you don't remenber it: http://vergil.chemistry.gatech.edu/notes/permsymm/permsymm.pdf	
	     for i in range(nmo):
		for j in range(0,i+1):
			for k in range(0,i+1):
				for l in range(0,k+1):
			 		v =  eri_kpt[i,j,k,l].real #TODO THIS IS ONLY A HACK FOR THE GAMA POINT BUG (aka we need to a perturbation)
								   # We need to do the cafaral perturbation one day
                         		if (i*j >= k*l) and abs(v) > int_threshold:
						yield (i+1,k+1,j+1,l+1,v)


def get_kpt_idx(kpt,kpts):
	return next(i for i,j in enumerate(kpts) if all( k==l for k,l in zip(j,kpt)))

class Conveniant_notation(object):

	def __init__(self, q_kpts,q_kpts_inv):
		self.q_kpts = q_kpts
		self.q_kpts_inv = q_kpts_inv
	
	def get_kpt(self,i):
		array = self.q_kpts if i>0 else self.q_kpts_inv
		return array[abs(i)-1]	

	def get(self,t1,t2):
		k1, k3 = t1
		k2, k4 = t2

		q_kpts =  map(self.get_kpt, [k1,k2,k3,k4])
		qkpt_idx = [get_kpt_idx(kpt, kpts) for kpt in q_kpts] 		


		l_mo_coeff = [kernel.mo_coeff[i] for i in qkpt_idx]
	
                # Compute the ERI
                eri_kpt_flatish = kernel.with_df.ao2mo(l_mo_coeff,q_kpts)

      	        nmo = kernel.mo_coeff.shape[1]
                # Putin in a good shape
                return eri_kpt_flatish.reshape((nmo,)*4)
	
		
def netherland_obliteration(q_kpts,q_kpts_inv,q_sigma):
	#Chimist notation
	#s = +1 or -1
	print q_kpts
	print q_kpts_inv

	s1,s2,s3,s4 = q_sigma
	p = Conveniant_notation(q_kpts,q_kpts_inv)

	zero_sigma  =     p.get((1,2),(3,4))
	print zero_sigma
	print '--------'

	one_sigma  =  s1*p.get((-1,2),(3,4)) + s2*p.get((1,-2),(3,4)) + \
		      s3*p.get((1,2),(-3,4)) + s4*p.get((1,2),(3,-4))

	two_sigma = s1*s2*p.get((-1,-2),(3,4)) + s1*s3*p.get((-1,2),(-3,4)) + \
		    s1*s4*p.get((1,-2),(3,-4)) + s2*s3*p.get((1,-2),(-3,4)) + \
		    s2*s4*p.get((1,-2),(3,-4)) + s3*s4*p.get((1,2),(-3,-4))

        three_sigma  =  s1*s2*s3*p.get((-1,-2),(-3,4)) + s1*s2*s4*p.get((-1,-2),(3,-4)) + \
                        s1*s3*s4*p.get((-1,2),(-3,-4)) + s2*s3*s4*p.get((1,-2),(-3,-4))
			
	four_sigma = s1*s2*s3*s4*p.get((-1,-2),(-3,-4))
 

	l_nu = [1 if i == 1 else 1j for i in  [s1,s2,s3,s4]]
	return (zero_sigma*two_sigma*three_sigma*four_sigma)/(16*numpy.prod(l_nu)) 


# Print all the eri
# Note: Never tester with more than one K point.  
kconserv = tools.get_kconserv(cell, kpts) # Get the momentum conservation array for a set of k-points.
with open('bielec_mo','w') as f:
	for kp,kq,kr in product(range(nkpt),repeat=3): # Get all the kp,kq,kr possible
            ks = kconserv[kp,kq,kr] # Get the ks
	    for eri in gen_eri_MO(kmf,kpts,(kp,kq,kr,ks)):
                  f.write('%s %s %s %s %s\n' % eri)
	    

# TODO: CLEANING AND PRINTING
#kernel = kmf
#print netherland_obliteration([kpts[0]]*4,[kpts[1]]*4,[1,1,1,1])

