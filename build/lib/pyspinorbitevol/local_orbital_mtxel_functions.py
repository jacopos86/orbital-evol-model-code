#
#   set of isolated atom orbital
#   momentum matrix elements:
#   1)  <l1, ml1|L+|l2, ml2>
#   2)  <l1, ml1|L-|l2, ml2>
#   3)  <l1, ml1|Lx|l2, ml2>
#   4)  <l1, ml1|Ly|l2, ml2>
#   5)  <l1, ml1|Lz|l2, ml2>
#
import numpy as np
from pyspinorbitevol.utility_functions import delta
#
#   function (1)
#
def L0plus_mtxel(l1, ml1, l2, ml2):
	r = np.sqrt((l2 - ml2) * (l2 + ml2 + 1)) * delta(l1, l2) * delta(ml1, ml2+1)
	return r
#
#   function (2)
#
def L0minus_mtxel(l1, ml1, l2, ml2):
	r = np.sqrt((l2 + ml2) * (l2 - ml2 + 1)) * delta(l1, l2) * delta(ml1, ml2-1)
	return r
#
#   function (3)
#
def L0x_mtxel(l1, ml1, l2, ml2):
	rp = L0plus_mtxel(l1, ml1, l2, ml2)
	rm = L0minus_mtxel(l1, ml1, l2, ml2)
	r = (rp + rm) / 2.
	return r
#
#   function (4)
#
def L0y_mtxel(l1, ml1, l2, ml2):
	rp = L0plus_mtxel(l1, ml1, l2, ml2)
	rm = L0minus_mtxel(l1, ml1, l2, ml2)
	r = (rp - rm) / (2.*1j)
	return r
#
#   function (5)
#
def L0z_mtxel(l1, ml1, l2, ml2):
	r = ml1 * delta(l1, l2) * delta(ml1, ml2)
	return r
