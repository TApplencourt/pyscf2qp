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
L=3.5
cell = gto.M(
    atom = '''He 0. 0. 0.''',
    basis = '6-31g',
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
# 	   So please shift the gamma point by epsilon.

kpts_rel     = [ [0.5, 0.5,  0.] ]
kpts_rel_inv = [ [-1*i for i in k] for k in kpts_rel ]

kpts_abs    = [cell.get_abs_kpts(k) for k in kpts_rel]
kpts_abs_inv= [cell.get_abs_kpts(k) for k in kpts_rel_inv]

kpts = numpy.array(kpts_abs+kpts_abs_inv)

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
	mydf.auxbasis = 'weigend' # Note: No idea what 'weigend' is this. But the use it all the time in the example.
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


# ~!
# Now we will compute the needed integral
# In Quantum Package we can only use REAL integral, so we will do some rotation
# Note. We will rotate the integral and note the orbital. This is really inefficient!
# ~!
def get_kpt_idx(kpt,kpts):
        return next(i for i,j in enumerate(kpts) if all( k==l for k,l in zip(j,kpt)))

# ___                                              
#  |  ._ _|_  _   _  ._ _. |  _   |\/|  _  ._   _  
# _|_ | | |_ (/_ (_| | (_| | _>   |  | (_) | | (_) 
#                 _|                              

def mono_netherland_obliteration(kpt,kpt_inv,sigma,l_int):
	#By symetrie only k1 | k2 is null if k1 != k2
	#This is why we use kpt and do d_kpt 

        def ao2mo(kpt):
                idx_kpt = get_kpt_idx(kpt,kpts)
                return reduce(numpy.dot, (mo_coeff[idx_kpt].T, l_int[idx_kpt], mo_coeff[idx_kpt]))

        int_ = ao2mo(kpt)
        int_inv = ao2mo(kpt_inv)

	print int_
	print int_inv
	print ''
        nu = 1 if sigma == 1 else 1j
        l_int = (int_+ sigma*int_inv) / ( 2*nu)
	# Need to check if this is real
	
	return numpy.real(l_int)

from itertools import product

def print_MO(l_int_mo,padding=0):
     '''Print erything do not care about he symetrie'''
     for i,j in product(range(nmo), repeat=2):
            int_ = l_int_mo[i,j]
            if abs(int_) > int_threshold:
                    yield (padding+i+1,padding+j+1, int_)

# A tuple with the name, and the value of the integrals for each k-pts
t_ao = ('kinetic', cell.pbc_intor('cint1e_kin_sph',kpts=kpts)) #  Note .intor function computes integrals without PBC
v_ao = ('nuclear', [kmf.with_df.get_nuc(kpt) for kpt in kpts]) # Ref:  `get_hcore` in pbc/scf file


for kpt, kpt_inv in zip(kpts_abs,kpts_abs_inv):

	for name, l_int_ao in (t_ao,v_ao):	
		l_mo = mono_netherland_obliteration(kpt,kpt_inv, 1,l_int_ao)
	
		with open('%s_mo' % name,'w') as f:
			l_mo = mono_netherland_obliteration(kpt,kpt_inv, 1,l_int_ao)
                        for mono in print_MO(l_mo, padding=nmo*get_kpt_idx(kpt,kpts)):
                                f.write('%s %s %s\n'% mono)

	                l_mo = mono_netherland_obliteration(kpt,kpt_inv, -1,l_int_ao)
                        for mono in print_MO(l_mo, padding=nmo*nkpt*get_kpt_idx(kpt,kpts)):
                                f.write('%s %s %s\n'% mono)

# ___                              _    
#  |  ._ _|_  _   _  ._ _. |  _   |_) o 
# _|_ | | |_ (/_ (_| | (_| | _>   |_) | 
#                 _|                    
#

class Conveniant_notation(object):

	def __init__(self, q_kpts,q_kpts_inv):
		self.q_kpts = q_kpts
		self.q_kpts_inv = q_kpts_inv

		self.cache = dict()	
	def get_kpt(self,i):
		array = self.q_kpts if i>0 else self.q_kpts_inv
		return array[abs(i)-1]	



	def get(self,t1,t2):
		
		if (t1,t2) in self.cache:
			return self.cache[ (t1,t2) ]

		k1, k3 = t1
		k2, k4 = t2

		q_kpts =  map(self.get_kpt, [k1,k2,k3,k4])
		qkpt_idx = [get_kpt_idx(kpt, kpts) for kpt in q_kpts] 		


		l_mo_coeff = [kmf.mo_coeff[i] for i in qkpt_idx]
	
                # Compute the ERI
                eri_kpt_flatish = kmf.with_df.ao2mo(l_mo_coeff,q_kpts)

      	        nmo = kmf.mo_coeff.shape[1]
                # Putin in a good shape
		eri = eri_kpt_flatish.reshape((nmo,)*4)

		self.cache[ (t1,t2) ] = eri
                return eri

def bi_netherland_obliteration(p,q_sigma):
	#Chimist notation
	#s = +1 or -1

	s1,s2,s3,s4 = q_sigma

	zero_sigma  =     p.get((1,2),(3,4))

	one_sigma  =  s1*p.get((-1,2),(3,4)) + s2*p.get((1,-2),(3,4)) + \
		      s3*p.get((1,2),(-3,4)) + s4*p.get((1,2),(3,-4))

	two_sigma = s1*s2*p.get((-1,-2),(3,4)) + s1*s3*p.get((-1,2),(-3,4)) + \
		    s1*s4*p.get((1,-2),(3,-4)) + s2*s3*p.get((1,-2),(-3,4)) + \
		    s2*s4*p.get((1,-2),(3,-4)) + s3*s4*p.get((1,2),(-3,-4))

        three_sigma  =  s1*s2*s3*p.get((-1,-2),(-3,4)) + s1*s2*s4*p.get((-1,-2),(3,-4)) + \
                        s1*s3*s4*p.get((-1,2),(-3,-4)) + s2*s3*s4*p.get((1,-2),(-3,-4))
			
	four_sigma = s1*s2*s3*s4*p.get((-1,-2),(-3,-4))
 

	l_nu = [1 if i == 1 else 1j for i in  q_sigma]
	return (zero_sigma+two_sigma+three_sigma+four_sigma)/(16*numpy.prod(l_nu)) 


def gen_eri_MO(eri_kpt, q_padding):

             # Warning: We need to handle the 8-fold symetrie. If not this will bug in the QP size.
             # If you don't remenber it: http://vergil.chemistry.gatech.edu/notes/permsymm/permsymm.pdf
             for i in range(nmo):
                for j in range(0,i+1):
                        for k in range(0,i+1):
                                for l in range(0,k+1):
                                        v =  eri_kpt[i,j,k,l] 
                                        if (i*j >= k*l) and abs(v) > int_threshold:
                                                yield (q_padding[0]+i+1,
						       q_padding[2]+k+1,
						       q_padding[1]+j+1,
						       q_padding[3]+l+1,
						       v)

def get_padding(kpt,i):
	return nmo*kpt if i == 1 else nmo*(nkpt/2+kpt)

kconserv = tools.get_kconserv(cell, numpy.array(kpts_abs))
with open('bielec_mo','w') as f:
	for kp,kq,kr in product(range(len(kpts_abs)),repeat=3): # Get all the kp,kq,kr possible
            ks = kconserv[kp,kq,kr] # Get the ks


	    q_kpts = [kpts[kp],
		      kpts[kq],
		      kpts[kr],
		      kpts[ks] ]
	    q_kpts_inv = [ [-1*i for i in k] for k in q_kpts]


	    print q_kpts
	    print q_kpts_inv
	    print ''
	    engine = Conveniant_notation(q_kpts,q_kpts_inv)
	    #Now do the permutation
	    for p in [1, -1]:
		p_pad = get_padding(kp,p)
		for q in [1, -1]:
			q_pad = get_padding(kq,q)
			for r in [1, -1]:
				r_pad = get_padding(kr,r)
				for s in [1, -1]:
					s_pad = get_padding(ks,s)

					q_sigma = [p,q,r,s]
					print q_sigma
					print (p_pad,q_pad,r_pad,s_pad)
					l_int = bi_netherland_obliteration(engine,q_sigma)
					
					for eri in gen_eri_MO(l_int,(p_pad,q_pad,r_pad,s_pad)):
						 f.write('%s %s %s %s %s\n' % eri) 
						 print eri
				
					print ''
