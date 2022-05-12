#
#   This module defines the orbital momentum
#   class operators -
#   it implements
#   1) the atom localized orbital momentum operator (L0)
#   2) the conserved orbital momentum operator (COM)
#   3) the conserved COM with phonons (COMph)
#
import numpy as np
import sys
from pyspinorbitevol.orbital_mtxel_functions import L0x_mtxel, L0y_mtxel, L0z_mtxel
#
class OrbitalMomentumOperators:
	def __init__(self, siteslist, kg, MatrixEntry, phonons_calc=False):
		# set orbital operators
		self.set_orbital_operators(siteslist, kg)
		# set L0
		self.set_L0(siteslist, kg, MatrixEntry)
		# set conserved COM
		# if phonons -> compute COM_ph
	def set_orbital_operators(self, siteslist, kg):
		Nst = siteslist.Nst
		# initialize L0 -> atom localized OM
		# initialize COM -> COM
		# initialize COM_ph -> COM + phonon orbital momentum
		if kg.D == 0:
			self.L0 = np.zeros((Nst, Nst, 3), dtype=np.complex128)
			self.COM = np.zeros((Nst, Nst, 3), dtype=np.complex128)
			self.COM_ph = np.zeros((Nst, Nst, 3), dtype=np.complex128)
		elif kg.D == 1:
			[nk1, nk2, nk3] = kg.nkpts
			if nk1 != 0:
				self.L0 = np.zeros((Nst, Nst, nk1, 3), dtype=np.complex128)
				self.COM= np.zeros((Nst, Nst, nk1, 3), dtype=np.complex128)
				self.COM_ph= np.zeros((Nst, Nst, nk1, 3), dtype=np.complex128)
			elif nk2 != 0:
				self.L0 = np.zeros((Nst, Nst, nk2, 3), dtype=np.complex128)
				self.COM= np.zeros((Nst, Nst, nk2, 3), dtype=np.complex128)
				self.COM_ph= np.zeros((Nst, Nst, nk2, 3), dtype=np.complex128)
			elif nk3 != 0:
				self.L0 = np.zeros((Nst, Nst, nk3, 3), dtype=np.complex128)
				self.COM= np.zeros((Nst, Nst, nk3, 3), dtype=np.complex128)
				self.COM_ph= np.zeros((Nst, Nst, nk3, 3), dtype=np.complex128)
			else:
				print("wrong n. k-pts for D=1")
				sys.exit(1)
		elif kg.D == 2:
			[nk1, nk2, nk3] = kg.nkpts
			if nk1 == 0:
				self.L0 = np.zeros((Nst, Nst, nk2, nk3, 3), dtype=np.complex128)
				self.COM= np.zeros((Nst, Nst, nk2, nk3, 3), dtype=np.complex128)
				self.COM_ph= np.zeros((Nst, Nst, nk2, nk3, 3), dtype=np.complex128)
			elif nk2 == 0:
				self.L0 = np.zeros((Nst, Nst, nk1, nk3, 3), dtype=np.complex128)
				self.COM= np.zeros((Nst, Nst, nk1, nk3, 3), dtype=np.complex128)
				self.COM_ph= np.zeros((Nst, Nst, nk1, nk3, 3), dtype=np.complex128)
			elif nk3 == 0:
				self.L0 = np.zeros((Nst, Nst, nk1, nk2, 3), dtype=np.complex128)
				self.COM= np.zeros((Nst, Nst, nk1, nk2, 3), dtype=np.complex128)
				self.COM_ph= np.zeros((Nst, Nst, nk1, nk2, 3), dtype=np.complex128)
			else:
				print("wrong n. k-pts for D=2")
				sys.exit(1)
		elif kg.D == 3:
			[nk1, nk2, nk3] = kg.nkpts
			self.L0 = np.zeros((Nst, Nst, nk1, nk2, nk3, 3), dtype=np.complex128)
			self.COM= np.zeros((Nst, Nst, nk1, nk2, nk3, 3), dtype=np.complex128)
			self.COM_ph= np.zeros((Nst, Nst, nk1, nk2, nk3, 3), dtype=np.complex128)
	# compute L0 here
	def set_L0(self, siteslist, kg, MatrixEntry):
		for i in range(siteslist.Nsites):
			site = i+1
			for l1 in siteslist.Atomslist[i].OrbitalList:
				for ml1 in range(-l1, l1+1):
					for l2 in siteslist.Atomslist[i].OrbitalList:
						for ml2 in range(-l2, l2+1):
							for ms in [-0.5, 0.5]:
								row = MatrixEntry(siteslist.Atomslist, site, l1, ml1, ms)
								col = MatrixEntry(siteslist.Atomslist, site, l2, ml2, ms)
								if kg.D == 0:
									self.L0[row,col,0] = L0x_mtxel(l1, ml1, l2, ml2)
									self.L0[row,col,1] = L0y_mtxel(l1, ml1, l2, ml2)
									self.L0[row,col,2] = L0z_mtxel(l1, ml1, l2, ml2)
								elif kg.D == 1:
									self.L0[row,col,:,0] = L0x_mtxel(l1, ml1, l2, ml2)
									self.L0[row,col,:,1] = L0y_mtxel(l1, ml1, l2, ml2)
									self.L0[row,col,:,2] = L0z_mtxel(l1, ml1, l2, ml2)
								elif kg.D == 2:
									self.L0[row,col,:,:,0] = L0x_mtxel(l1, ml1, l2, ml2)
									self.L0[row,col,:,:,1] = L0y_mtxel(l1, ml1, l2, ml2)
									self.L0[row,col,:,:,2] = L0z_mtxel(l1, ml1, l2, ml2)
								elif kg.D == 3:
									self.L0[row,col,:,:,:,0] = L0x_mtxel(l1, ml1, l2, ml2)
									self.L0[row,col,:,:,:,1] = L0y_mtxel(l1, ml1, l2, ml2)
									self.L0[row,col,:,:,:,2] = L0z_mtxel(l1, ml1, l2, ml2)
								else:
									print("Error: wrong k grid size")
									sys.exit(1)
