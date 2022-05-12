import numpy as np
#
#   This module implements
#   the atomic site class
#
class AtomicSite:
	def __init__(self, Element, Mass, R0, V0, OrbitalList, nel):
		self.element = Element
		# mass in eV fs^2/Ang^2
		self.mass = Mass
		# R0 in Ang
		self.R0 = np.array(R0)
		# V0 in Ang/fs
		self.V0 = np.array(V0)
		# delta R
		self.dR0 = np.array([0., 0., 0.])
		# delta V
		self.dV0 = np.array([0., 0., 0.])
		# orbital list
		self.OrbitalList = OrbitalList
		# n. electrons on site
		self.nel = nel
	def update_position(self, R0):
		self.R0 = np.array(R0)
	def update_velocity(self, V0):
		self.V0 = np.array(V0)
	def set_orbital_momentum(self):
		self.L0 = self.mass * np.cross(self.R0, self.V0)
		# eV fs units
	def set_relative_position(self, Rcm):
		self.dR0 = self.R0 - Rcm
	def set_relative_velocity(self, Vcm):
		self.dV0 = self.V0 - Vcm
#
#   AtomicSiteList class
#
class AtomicSiteList:
	def __init__(self):
		self.Atomslist = []
		self.Nsites = 0
		self.Nst = 0
		self.Rcm = np.zeros(3)
		self.Vcm = np.zeros(3)
		self.M = 0.
		self.latt_orbital_mom = np.zeros(3) 
	def add_site_to_list(self, Element, Mass, R0, V0, OrbitalList, nel):
		# set atomic site
		site = AtomicSite(Element, Mass, R0, V0, OrbitalList, nel)
		site.set_orbital_momentum()
		# append site to list
		self.Atomslist.append(site)
	def set_number_of_sites(self):
		self.Nsites = len(self.Atomslist)
	def set_number_of_states(self):
		# compute total number of states in list
		Nst = 0
		for i in range(self.Nsites):
			llist = self.Atomslist[i].OrbitalList
			for j in range(len(llist)):
				l = llist[j]
				Nst = Nst + 2*(2*l+1)
		self.Nst = Nst
	def set_total_mass(self):
		for site in self.Atomslist:
			self.M = self.M + site.mass
	def set_center_of_mass_position(self):
		rcm = np.zeros(3)
		for site in self.Atomslist:
			rcm[:] = rcm[:] + site.mass * site.R0[:]
		rcm[:] = rcm[:] / self.M
		self.Rcm = rcm
	def set_center_of_mass_velocity(self):
		vcm = np.zeros(3)
		for site in self.Atomslist:
			vcm[:] = vcm[:] + site.mass * site.V0[:]
		vcm[:] = vcm[:] / self.M
		self.Vcm = vcm
	def set_relative_positions(self):
		for site in self.Atomslist:
			site.set_relative_position(self.Rcm)
	def set_relative_velocities(self):
		for site in self.Atomslist:
			site.set_relative_velocity(self.Vcm)
	def set_lattice_orbital_momentum(self):
		self.latt_orbital_mom[:] = 0.
		for site in self.Atomslist:
			self.latt_orbital_mom = self.latt_orbital_mom + site.L0
