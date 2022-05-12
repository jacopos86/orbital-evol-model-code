import numpy as np
from pyspinorbitevol.crystal_field_class import CrystalFieldHamilt
from pyspinorbitevol.unit_cell_class import UnitCell
from pyspinorbitevol.KPoints_class import KPoints
from pyspinorbitevol.atomic_sites_class import AtomicSiteList
from pyspinorbitevol.phys_constants import hbar
from pyspinorbitevol.utility_functions import MatrixEntry
from pyspinorbitevol.orbital_operators_class import OrbitalMomentumOperators
from pyspinorbitevol.spin_operators_class import SpinMomentumOperators
from pyspinorbitevol.onsite_energies import set_onsite_energies
from pyspinorbitevol.crystal_field_gradient_class import CrystalFieldHamiltGradient
from pyspinorbitevol.ground_state_calc import GroundState
from pyspinorbitevol.expect_values import spin_expect_val, orbital_expect_val
import matplotlib.pyplot as plt
#
#uc = UnitCell([[1.,1.,-1.],[-1.,1.,1.],[1.,-1.,1.]])
uc = UnitCell([[2.,0.,0.],[0.,2.,0.],[0.,0.,2.]])
uc.set_volume()
uc.set_rec_vectors()
uc.set_rec_versors()
#
kg = KPoints(3, [3,3,3])
kg.set_kgrid(uc)
kg.set_kpts_weights()
#
sites_list = AtomicSiteList()
sites_list.add_site_to_list('Fe', 1.0, [0., 0., 0.], [0., 0., 0.], [0, 1, 2], 8)
#sites_list.add_site_to_list('Gd', 1.0, [1., 0., 0.], [0., 0., 0.], [0, 1, 2, 3])
sites_list.set_number_of_sites()
sites_list.set_number_of_states()
sites_list.set_total_mass()
sites_list.set_center_of_mass_position()
sites_list.set_center_of_mass_velocity()
sites_list.set_relative_positions()
sites_list.set_relative_velocities()
sites_list.set_lattice_orbital_momentum()
print(sites_list.latt_orbital_mom)
#
uc.set_structure(sites_list.Atomslist, kg)
uc.set_nn_atoms(sites_list.Atomslist, kg)
t = {('Fe','Fe') : {'ss' : 0.1, 'sp' : 0.2, 'ps' : 0.2, 'pp' : [0.2,0.3], 'sd' : 0.05, 'ds' : 0.05, 'pd' : [0.1, 0.2], 'dp' : [0.1, 0.2], 'dd' : [0.05, 0.06, 0.06]}}
#
#t = {('Fe','Gd') : {'ss' : 0.1, 'sp' : 0.2, 'ps' : 0.2, 'pp' : [0.2,0.3], 'sd' : 0.05, 'ds' : 0.05, 'pd' : [0.1, 0.2], 'dp' : [0.1, 0.2], 'dd' : [0.05, 0.06, 0.06], 'sf' : 0.0, 'pf' : [0.05, 0.05], 'df' : [0.01, 0.02, 0.005]}}
#atomic_energies = {'Fe' : {0 : -1.0, 1 : 12.0, 2 : -2.0}, 'Gd' : {0 : -1.0, 1 : 12.0, 2 :-2.0, 3 : -5.0}}
atomic_energies = {'Fe' : {0 : -1.0, 1 : 6.0, 2 : -2.0}}
CrystalFieldH = CrystalFieldHamilt(t, sites_list, kg, uc, MatrixEntry)
#
orbital_operators = OrbitalMomentumOperators(sites_list, kg, MatrixEntry)
spin_operators = SpinMomentumOperators(sites_list, kg, MatrixEntry)
print(spin_operators.S.shape)
set_onsite_energies(CrystalFieldH, sites_list, kg, MatrixEntry, atomic_energies)
#
#gt = {('Fe','Gd') : {'ss' : [0.1, 0.01]}}
#CrystalFieldGrad = CrystalFieldHamiltGradient(gt, sites_list, kg, uc, MatrixEntry)
#
GS = GroundState(CrystalFieldH, sites_list, kg, uc, MatrixEntry)
GS.plot_band_structure(kg, "./bands-plot.dat")
GS.set_wfc_occupations(sites_list, kg, 0.03)
GS.set_density_matrix(sites_list, kg)
#R = np.matmul(CrystalFieldH.H0, orbital_operators.L0[:,:,2]) - np.matmul(orbital_operators.L0[:,:,2], CrystalFieldH.H0)
#for i in range(sites_list.Nst):
#	for j in range(sites_list.Nst):
#		if abs(R[i,j]) > 1.E-7:
#			print(i,j,R[i,j])
dos = GS.compute_elec_DOS(kg, sites_list, "./DOS", [0., 10.], 0.05, 0.01)
Sexp = spin_expect_val(spin_operators, GS.rho_e, sites_list, kg, MatrixEntry)
orbital_expect_val(orbital_operators, GS.rho_e, sites_list, kg, MatrixEntry)
