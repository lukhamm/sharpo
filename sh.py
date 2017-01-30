# ===========================================================
#
#   Module for sharpo
#   sh.py :: returns grid values for the spherical harmonic 
#            functions Ylm
#
#   License
#   :: Copyright (c) 2016 by Lukas Hammerschmidt
#      l.hammerschmidt@auckland.ac.nz
#
#   :: This module is part of sharpo. Sharpo is free software: you
#      can redistribute it and/or modify it under the terms of the
#      GNU Lesser General Public License as published by the Free
#      Software Foundation, either version 3 of the License, or any
#      later version.
#
#   :: This program is distributed in the hope that it will be useful,
#      but WITHOUT ANY WARRANTY; without even the implied warranty of
#      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#      GNU Lesser General Public License for more details. You should
#      have received a copy of the GNU Lesser General Public License
#      along with sharpo.  If not, see <http://www.gnu.org/licenses/>.

#   Parameters
#   [int]         :: l     :: angular momentum quantum number
#   [int]         :: m     :: magnetic quantum number
#   [numpy array] :: theta :: contains the theta angle grid
#   [numpy array] :: phi   :: contains the phi angle grid
#
#  Comments
#  :: Couldn't be bothered to figure out how to implement the
#     recursion formulars for Plm in python. So all spherical
#     harmonics are explicitely typed in. Change at some point
#     maybe.
#  :: lmax == 6 at the moment!!!
#
#  Last edited
#  13.09.2016 :: Descriptions and finalization
#
# ===========================================================

import numpy as np

def sh(l,m,theta,phi):
  #-----------------------------------
  # l = 0 (s-type) spherical harmonics
  #-----------------------------------
  if l == 0 and m == 0:
    return (1./2.)*np.sqrt(1./np.pi)
  #-----------------------------------
  # l = 1 (p-type) spherical harmonics
  #-----------------------------------
  if l == 1 and m == -1:
    return (1./2.) * np.sqrt(3./(2.*np.pi)) * np.exp(-1.j*phi) * np.sin(theta)
  if l == 1 and m == 0:
    return (1./2.)*np.sqrt(3./np.pi) * np.cos(theta)
  if l == 1 and m == 1:
    return -(1./2.) * np.sqrt(3./(2.*np.pi)) * np.exp(1.j*phi) * np.sin(theta)
  #-----------------------------------
  # l = 2 (d-type) spherical harmonics
  #-----------------------------------
  if l == 2 and m == -2:
    return (1./4.) * np.sqrt(15./(2.*np.pi)) * np.exp(-2.j*phi) * np.sin(theta)**2.
  if l == 2 and m == -1:
    return (1./2.) * np.sqrt(15./(2.*np.pi)) * np.exp(-1.j*phi) * np.sin(theta)*np.cos(theta)
  if l == 2 and m == 0:
    return (1./4.) * np.sqrt(5./np.pi) * (3.*np.cos(theta)**2. - 1.)
  if l == 2 and m == 1:
    return -(1./2.) * np.sqrt(15./(2.*np.pi)) * np.exp(1.j*phi) * np.sin(theta)*np.cos(theta)
  if l == 2 and m == 2:
    return (1./4.) * np.sqrt(15./(2.*np.pi)) * np.exp(-2.j*phi) * np.sin(theta)**2.
  #-----------------------------------
  # l = 3 (f-type) spherical harmonics
  #-----------------------------------
  if l == 3 and m == -3:
    return (1./8.) * np.sqrt(35./np.pi) * np.exp(-3.j*phi) * np.sin(theta)**3.
  if l == 3 and m == -2:
    return (1./4.) * np.sqrt(105./(2.*np.pi)) * np.exp(-2.j*phi) * np.sin(theta)**2. * np.cos(theta)
  if l == 3 and m == -1:
    return (1./8.) * np.sqrt(21./np.pi) * np.exp(-1.j*phi) * np.sin(theta) * (5.*np.cos(theta)**2. - 1.)
  if l == 3 and m ==  0:
    return (1./4.) * np.sqrt(7./np.pi) * (5.*np.cos(theta)**3. - 3.*np.cos(theta))
  if l == 3 and m ==  1:
    return -(1./8.) * np.sqrt(21./np.pi) * np.exp(1.j*phi) * np.sin(theta) * (5.*np.cos(theta)**2. - 1.)
  if l == 3 and m ==  2:
    return (1./4.) * np.sqrt(105./(2.*np.pi)) * np.exp(2.j*phi) * np.sin(theta)**2. * np.cos(theta)
  if l == 3 and m ==  3:
    return -(1./8.) * np.sqrt(35./np.pi) * np.exp(3.j*phi) * np.sin(theta)**3.
  #-----------------------------------
  # l = 4 (g-type) spherical harmonics
  #-----------------------------------
  if l == 4 and m == -4:
    return (3./16.) * np.sqrt(35./(2.*np.pi)) * np.exp(-4.j*phi) * np.sin(theta)**4.
  if l == 4 and m == -3:
    return (3./8.) * np.sqrt(35./np.pi) * np.exp(-3.j*phi) * np.sin(theta)**3. * np.cos(theta)
  if l == 4 and m == -2:
    return (3./8.) * np.sqrt(5./(2.*np.pi)) * np.exp(-2.j*phi) * np.sin(theta)**2. * (7.*np.cos(theta)**2. - 1.)
  if l == 4 and m == -1:
    return (3./8.) * np.sqrt(5./np.pi) * np.exp(-1.j*phi) * np.sin(theta) * (7.*np.cos(theta)**3. - 3.*np.cos(theta))
  if l == 4 and m == 0:
    return (3./16.) * np.sqrt(1./np.pi) * (35. * np.cos(theta)**4. - 30.*np.cos(theta)**2. + 3.)
  if l == 4 and m == 1:
    return -(3./8.) * np.sqrt(5./np.pi) * np.exp(1.j*phi) * np.sin(theta) * (7.*np.cos(theta)**3. - 3.*np.cos(theta))
  if l == 4 and m == 2:
    return (3./8.) * np.sqrt(5./(2.*np.pi)) * np.exp(2.j*phi) * np.sin(theta)**2. * (7.*np.cos(theta)**2. - 1.)
  if l == 4 and m == 3:
    return -(3./8.) * np.sqrt(35./np.pi) * np.exp(3.j*phi) * np.sin(theta)**3. * np.cos(theta)
  if l == 4 and m == 4:
    return (3./16.) * np.sqrt(35./(2.*np.pi)) * np.exp(4.j*phi) * np.sin(theta)**4.
  #-----------------------------------
  # l = 5 (h-type) spherical harmonics
  #-----------------------------------
  if l == 5 and m == -5:
    return (3./32.) * np.sqrt(77./np.pi) * np.exp(-5.j*phi) * np.sin(theta)**5.
  if l == 5 and m == -4:
    return (3./16.) * np.sqrt(385./(2.*np.pi)) * np.exp(-4.j*phi) * np.sin(theta)**4. * np.cos(theta)
  if l == 5 and m == -3:
    return (1./32.) * np.sqrt(385./np.pi) * np.exp(-3.j*phi) * np.sin(theta)**3. * (9.*np.cos(theta)**2. - 1.)
  if l == 5 and m == -2:
    return (1./8.) * np.sqrt(1155./(2.*np.pi)) * np.exp(-2.j*phi) * np.sin(theta)**2. * (3.*np.cos(theta)**3. - np.cos(theta))
  if l == 5 and m == -1:
    return (1./16.) * np. sqrt(165./(2.*np.pi)) * np.exp(-1.j*phi) * np.sin(theta) * (21.*np.cos(theta)**4. - 14.*np.cos(theta)**2. + 1.)
  if l == 5 and m == 0:
    return (1./16.) * np.sqrt(11./np.pi) * (63.*np.cos(theta)**5. - 70.*np.cos(theta)**3. + 15.*np.cos(theta))
  if l == 5 and m == 1:
    return -(1./16.) * np. sqrt(165./(2.*np.pi)) * np.exp(1.j*phi) * np.sin(theta) * (21.*np.cos(theta)**4. - 14.*np.cos(theta)**2. + 1.)
  if l == 5 and m == 2:
    return (1./8.) * np.sqrt(1155./(2.*np.pi)) * np.exp(2.j*phi) * np.sin(theta)**2. * (3.*np.cos(theta)**3. - np.cos(theta))
  if l == 5 and m == 3:
    return -(1./32.) * np.sqrt(385./np.pi) * np.exp(3.j*phi) * np.sin(theta)**3. * (9.*np.cos(theta)**2. - 1.)
  if l == 5 and m == 4:
    return (3./16.) * np.sqrt(385./(2.*np.pi)) * np.exp(4.j*phi) * np.sin(theta)**4. * np.cos(theta)
  if l == 5 and m == 5:
    return -(3./32.) * np.sqrt(77./np.pi) * np.exp(5.j*phi) * np.sin(theta)**5.
  #-----------------------------------
  # l = 6 (i-type) spherical harmonics
  #-----------------------------------
  if l == 6 and m == -6:
    return (1./64.) * np.sqrt(3003./np.pi) * np.exp(-6.j*phi) * np.sin(theta)**6.
  if l == 6 and m == -5:
    return (3./32.) * np.sqrt(1001./np.pi) * np.exp(-5.j*phi) * np.sin(theta)**5. * np.cos(theta)
  if l == 6 and m == -4:
    return (3./32.) * np.sqrt(91./(2.*np.pi)) * np.exp(-4.j*phi) * np.sin(theta)**4. * (11.*np.cos(theta)**2. - 1.)
  if l == 6 and m == -3:
    return (1./32.) * np.sqrt(1365./np.pi) * np.exp(-3.j*phi) * np.sin(theta)**3. * (11.*np.cos(theta)**3. - 3.*np.cos(theta))
  if l == 6 and m == -2:
    return (1./64.) * np.sqrt(1365./np.pi) * np.exp(-2.j*phi) * np.sin(theta)**2. * (33.*np.cos(theta)**4. - 18.*np.cos(theta)**2. + 1.)
  if l == 6 and m == -1:
    return (1./16.) * np.sqrt(273./(2.*np.pi)) * np.exp(-1.j*phi) * np.sin(theta) * (33.*np.cos(theta)**5. - 30.*np.cos(theta)**3. + 5.*np.cos(theta))
  if l == 6 and m == 0:
    return (1./32.) * np.sqrt(13./np.pi) * (231. * np.cos(theta)**6. - 315. * np.cos(theta)**4. + 105.*np.cos(theta)**2. - 5.)
  if l == 6 and m == 1:
    return -(1./16.) * np.sqrt(273./(2.*np.pi)) * np.exp(1.j*phi) * np.sin(theta) * (33.*np.cos(theta)**5. - 30.*np.cos(theta)**3. + 5.*np.cos(theta))
  if l == 6 and m == 2:
    return (1./64.) * np.sqrt(1365./np.pi) * np.exp(2.j*phi) * np.sin(theta)**2. * (33.*np.cos(theta)**4. - 18.*np.cos(theta)**2. + 1.)
  if l == 6 and m == 3:
    return -(1./32.) * np.sqrt(1365./np.pi) * np.exp(3.j*phi) * np.sin(theta)**3. * (11.*np.cos(theta)**3. - 3.*np.cos(theta))
  if l == 6 and m == 4:
    return (3./32.) * np.sqrt(91./(2.*np.pi)) * np.exp(4.j*phi) * np.sin(theta)**4. * (11.*np.cos(theta)**2. - 1.)
  if l == 6 and m == 5:
    return -(3./32.) * np.sqrt(1001./np.pi) * np.exp(5.j*phi) * np.sin(theta)**5. * np.cos(theta)
  if l == 6 and m == 6:
    return (1./64.) * np.sqrt(3003./np.pi) * np.exp(6.j*phi) * np.sin(theta)**6.
  #-----------------------------------
  # l = 7 (j-type) spherical harmonics
  #-----------------------------------
  # ...



#-----------------------------------
# Real linear combinations: (for test purposes only!)
#-----------------------------------
# possible extention
# ...

