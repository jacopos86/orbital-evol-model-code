#
#  This module defines the system's
#  unit cell class
#  with number of periodic dimensions : D
#  lattice vectors
#  reciprocal vectors
#
import numpy as np
from pymatgen.core.lattice import Lattice
from pymatgen.core.structure import Structure
from pymatgen.analysis.local_env import CrystalNN
from pyspinorbitevol.utility_functions import norm_realv
from pyspinorbitevol.phys_constants import eps
import warnings
#
class UnitCell:
	def __init__(self, primitive_vectors):
		# primitive vectors must be in angstrom
		self.lattice = Lattice(primitive_vectors)
		self.a1 = np.array(primitive_vectors[0])
		self.a2 = np.array(primitive_vectors[1])
		self.a3 = np.array(primitive_vectors[2])
		self.prim_vecs = np.array([self.a1, self.a2, self.a3])
	def set_volume(self):
		a23 = np.cross(self.a2, self.a3)
		self.volume = np.dot(self.a1, a23)
		# Ang^3
	def set_rec_vectors(self):
		self.rec_lattice = self.lattice.reciprocal_lattice
		# b1 = 2pi a2 x a3 / V (ang^-1)
		self.b1 = 2.*np.pi*np.cross(self.a2, self.a3) / self.volume
		# b2 = 2pi a3 x a1 / V
		self.b2 = 2.*np.pi*np.cross(self.a3, self.a1) / self.volume
		# b3 = 2pi a1 x a2 / V
		self.b3 = 2.*np.pi*np.cross(self.a1, self.a2) / self.volume
		self.rec_vecs = np.array([self.b1, self.b2, self.b3])
	def set_rec_versors(self):
		self.rcv = [None]*3
		self.rcv[0] = self.b1 / norm_realv(self.b1)
		self.rcv[1] = self.b2 / norm_realv(self.b2)
		self.rcv[2] = self.b3 / norm_realv(self.b3)
	def set_structure(self, Atomslist, kgrid):
		# set up the full atomic structure
		species = []
		for Site in Atomslist:
			elem = Site.element
			species.append(elem)
		coords = []
		for Site in Atomslist:
			R = Site.R0
			coords.append(R)
		# set up the structure
		self.struct = Structure(lattice=self.lattice, species=species, coords=coords,
			charge=0, validate_proximity=True, coords_are_cartesian=True)
	def set_nn_atoms(self, Atomslist, kgrid):
		self.NNlist = []
		# periodic dimension 0
		if kgrid.D == 0:
			# suppress warnings
			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				NNstruct = CrystalNN()
				for i in range(len(Atomslist)):
					nndata = NNstruct.get_nn_info(self.struct, i)
					# store ONLY first nn in list (weight = 1)
					# check if the atom is in the structure
					nndata2 = []
					for data in nndata:
						R0 = Atomslist[data['site'].index].R0
						d = data['site'].coords - R0
						if norm_realv(d) < eps:
							#print(i, data['site'].coords, data['site'].index)
							nndata2.append(data)
					# store ONLY first nn in list (weight = 1)
					self.NNlist.append(nndata2)
		# periodic dimension 1
		elif kgrid.D == 1:
			i0 = np.where(np.array(kgrid.nkpts) > 0)[0][0]
			# lattice vector
			L = self.prim_vecs[i0]
			# suppress crystalNN warnings
			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				NNstruct = CrystalNN()
				for i in range(len(Atomslist)):
					nndata = NNstruct.get_nn_info(self.struct, i)
					# store ONLY first nn in list (weight = 1)
					# check if the atom is in the structure
					# accounting for D < 3
					nndata2 = []
					for data in nndata:
						R0 = Atomslist[data['site'].index].R0
						for k in [-2, -1, 0, 1, 2]:
							d = data['site'].coords - R0 + k * L
							if norm_realv(d) < eps:
								#print(i, data['site'].coords, data['site'].index, d)
								nndata2.append(data)
					self.NNlist.append(nndata2)
		# periodic dimension 2
		elif kgrid.D == 2:
			[i0, i1] = np.where(np.array(kgrid.nkpts) > 0)[0]
			# lattice vectors
			L0 = self.prim_vecs[i0]
			L1 = self.prim_vecs[i1]
			# crystalNN calculation
			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				NNstruct = CrystalNN()
				for i in range(len(Atomslist)):
					nndata = NNstruct.get_nn_info(self.struct, i)
					# store ONLY first nn in list (weight = 1)
					# check if the atom is in the structure
					# accounting for D < 3
					nndata2 = []
					for data in nndata:
						R0 = Atomslist[data['site'].index].R0
						for k0 in [-2, -1, 0, 1, 2]:
							for k1 in [-2, -1, 0, 1, 2]:
								d = data['site'].coords - R0 + k0 * L0 + k1 * L1
								if norm_realv(d) < eps:
									nndata2.append(data)
					self.NNlist.append(nndata2)
		# periodic dimension 3
		elif kgrid.D == 3:
			# crystalNN calculation
			with warnings.catch_warnings():
				warnings.simplefilter("ignore")
				NNstruct = CrystalNN()
				for i in range(len(Atomslist)):
					nndata = NNstruct.get_nn_info(self.struct, i)
					# store ONLY first nn in list (weight = 1)
					# check if the atom is in the structure
					# accounting for D < 3
					nndata2 = []
					for data in nndata:
						nndata2.append(data)
					self.NNlist.append(nndata2)
