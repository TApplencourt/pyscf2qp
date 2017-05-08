# pyscf2qp
Compute the Molecular integrals with PBC by PySCF and save it for Quantum Package.

Remarks:
- No Pseudo for now
- Only gamme point calculation supported
- Only clossed shell are supported.

# Requiremenet

- ![pyscf](https://github.com/sunqm/pyscf)
- QP_Plugin `read_integral` installed  (`qp_module.py install read_integral`) 

# How to

 - Modify `pyscf_fcidump.py` with the correct cell
 - Run `gen.py $ezfio` 
 - Done ! You can do your CI calculuation as usuall (`qp_run fci_zmq $ezfio`  for example)


## Hacking.
   1. `gen.py` will run `pyscf_fcidump.py`.
   2.  `pyscf_fcidump.py` will generate integrals and other information (MO_tot_num, number of electron...)
   3.  `create_ezfio.py` is called to generated a dummy EZFIO folder
   4.  `read_integrals_mo` is called to population the EZFIO folder with the integral
