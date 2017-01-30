# ===========================================================
#
#   Module for sharpo
#   mo_coefficients.py :: intended to write out the absolute 
#                         values of the MO-coefficients
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

#   Contains Functions:
#   :: make_sure_path_exists(path)
#      Parameters
#      [string]   :: path  :: folder to write into

#   :: Get_MO_Coefficients(qc, print_coeff_only, spin_polarized, pointsPerEv, sigma)
#      Parameters
#      [orbkit dictionary] :: qc               :: the wave function is saved in here
#      [bool]              :: print_coeff_only :: true -> program stops after printing coeff
#      [bool]              :: spin_polarized   :: true -> alpha and beta are evaluated seperately
#      [int]               :: pointsPerEv      :: #of grid points per eV for smearing
#      [float]             :: sigma            :: smearing factor

#  Comments
#  :: no comments
#
#  Last edited
#  26.01.2017 :: Optimizations / adaption of emax and emin
#  20.01.2017 :: Descriptions and finalization
#
# ===========================================================
from __future__ import print_function

import numpy
import errno
import time
import os
import sys

#--------------------------------------------
# Make sure a path exists
#--------------------------------------------
def make_sure_path_exists(path):
  try:
    os.makedirs(path)
  except OSError as exception:
    if exception.errno != errno.EEXIST:
      raise


#--------------------------------------------
# Write out MO coefficients
#--------------------------------------------
def Get_MO_Coefficients(qc, print_coeff_only, spin_polarized, 
                        pointsPerEv, sigma, E_fermi, emin, emax, selected_MO):
  

  lmax = 5                                                      #: s,p,d,f,g (can be extended in future)
  no_of_atoms = qc.geo_spec.size/3
  coeff_array = numpy.zeros((len(selected_MO),no_of_atoms,lmax)) #: coeff_array[mo][atom][s,p,d,f,g:|coefficient|^2]
  xmin = emin                                                   #: as set in arguments
  xmax = emax                                                   #: as set in arguments
  ipolpt = int((xmax - xmin)*pointsPerEv)+1                     #: nr of grid points
  xstep = (xmax - xmin)/float(ipolpt)                           #: step from one to the other grid point
  x = numpy.zeros(ipolpt+1)                                     #: x[x-value]



  
  #---------------------------------------------------------
  #: orbkit saves the AOs and MOs in the dictionary
  #  qc.ao_spec and qc.mo_spec. Such that:
  #  qc.mo_spec(i,'energy') gives the energy
  #  of the i-th MO. For more infos please refer to the
  #  orbkit documentation.
  #: we use those data set to extract in the MO coefficients
  #  for the individual atomic orbitals in the following
  #---------------------------------------------------------
  start_time = time.time()
  
  print ('')
  print ('Reading and writing out smeared MO-coefficients between %.2f and %.2f eV' % (xmin, xmax))
  
  for mo in range(0,len(selected_MO)):                           #: scans through all MOs
    counter=0
    for i in range(0,len(qc.ao_spec)):                          #: scans through all atoms
      
      if (qc.ao_spec[i]['type'] == 's'):
        coeff_array[mo, qc.ao_spec[i]['atom'], 0] += numpy.sqrt(numpy.power(qc.mo_spec[selected_MO[mo]-1]['coeffs'][counter],2.0))
        counter += 1
        
      elif (qc.ao_spec[i]['type'] == 'p'):
        for l in xrange(0,3):
          coeff_array[mo, qc.ao_spec[i]['atom'], 1] += numpy.sqrt(numpy.power(qc.mo_spec[selected_MO[mo]-1]['coeffs'][counter],2.0))
          counter += 1
          
      elif (qc.ao_spec[i]['type'] == 'd'):
        for l in xrange(0,5):
          coeff_array[mo, qc.ao_spec[i]['atom'], 2] += numpy.sqrt(numpy.power(qc.mo_spec[selected_MO[mo]-1]['coeffs'][counter],2.0))
          counter += 1
        if (qc.ao_spherical == None):
          coeff_array[mo, qc.ao_spec[i]['atom'], 2] += numpy.sqrt(numpy.power(qc.mo_spec[selected_MO[mo]-1]['coeffs'][counter],2.0))
          counter += 1
          
      elif (qc.ao_spec[i]['type'] == 'f'):
        for l in xrange(0,7):
          coeff_array[mo, qc.ao_spec[i]['atom'], 3] += numpy.sqrt(numpy.power(qc.mo_spec[selected_MO[mo]-1]['coeffs'][counter],2.0))
          counter += 1
        if (qc.ao_spherical == None):
          for l in xrange(0,3):
            coeff_array[mo, qc.ao_spec[i]['atom'], 3] += numpy.sqrt(numpy.power(qc.mo_spec[selected_MO[mo]-1]['coeffs'][counter],2.0))
            counter += 1
      
      elif (qc.ao_spec[i]['type'] == 'g'):
        for l in range(0,9):
          coeff_array[mo, qc.ao_spec[i]['atom'], 4] += numpy.sqrt(numpy.power(qc.mo_spec[selected_MO[mo]-1]['coeffs'][counter],2.0))
          counter += 1
        if (qc.ao_spherical == None):
          for l in range(0,5):
            coeff_array[mo, qc.ao_spec[i]['atom'], 4] += numpy.sqrt(numpy.power(qc.mo_spec[selected_MO[mo]-1]['coeffs'][counter],2.0))
            counter += 1
            

                                                              #: initialises arrays for Gaussian smearing
  if not spin_polarized:
    gx     = numpy.zeros((ipolpt+1, no_of_atoms, lmax))       #: gx[x-value][no_of_atoms][s,p,d,f,g:|coefficient|^2]
    gx_all = numpy.zeros((ipolpt+1, lmax))                    #: gx[x-value][s,p,d,f,g:|coefficient|^2]
  else:
    gx_alpha = numpy.zeros((ipolpt+1, no_of_atoms, lmax))
    gx_beta  = numpy.zeros((ipolpt+1, no_of_atoms, lmax))
    gx_all_a = numpy.zeros((ipolpt+1, lmax))                  #: gx[x-value][s,p,d,f,g:|coefficient|^2]
    gx_all_b = numpy.zeros((ipolpt+1, lmax))
  
  #-----------------------------------------------------
  # Calculates Gaussian smearing (like DOS in sol.state)
  #---
  for i in range(0,ipolpt+1):
    x[i] = xmin + xstep*float(i)
  
  for l in xrange(0,lmax):
    for atn in xrange(0, no_of_atoms):
      for i_mos in xrange(0,len(selected_MO)):
        if not spin_polarized:
          gx[:,atn,l] += coeff_array[i_mos,atn,l] * (1./(sigma*numpy.sqrt(2.*numpy.pi))) * numpy.exp(-(1./(2.*sigma**2.))*(x[:] - 27.21138602*qc.mo_spec[selected_MO[i_mos]-1]['energy'])**2.)
          gx_all[:,l] += coeff_array[i_mos,atn,l] * (1./(sigma*numpy.sqrt(2.*numpy.pi))) * numpy.exp(-(1./(2.*sigma**2.))*(x[:] - 27.21138602*qc.mo_spec[selected_MO[i_mos]-1]['energy'])**2.)
        else:
          if "_a" in qc.mo_spec[selected_MO[i_mos]-1]["sym"]:
            gx_alpha[:,atn,l] += coeff_array[i_mos,atn,l] * (1./(sigma*numpy.sqrt(2.*numpy.pi))) * numpy.exp(-(1./(2.*sigma**2.))*(x[:] - 27.21138602*qc.mo_spec[selected_MO[i_mos]-1]['energy'])**2.)
            gx_all_a[:,l]     += coeff_array[i_mos,atn,l] * (1./(sigma*numpy.sqrt(2.*numpy.pi))) * numpy.exp(-(1./(2.*sigma**2.))*(x[:] - 27.21138602*qc.mo_spec[selected_MO[i_mos]-1]['energy'])**2.)
          if "_b" in qc.mo_spec[selected_MO[i_mos]-1]["sym"]:
            gx_beta[:,atn,l] += coeff_array[i_mos,atn,l] * (1./(sigma*numpy.sqrt(2.*numpy.pi))) * numpy.exp(-(1./(2.*sigma**2.))*(x[:] - 27.21138602*qc.mo_spec[selected_MO[i_mos]-1]['energy'])**2.)
            gx_all_b[:,l]    += coeff_array[i_mos,atn,l] * (1./(sigma*numpy.sqrt(2.*numpy.pi))) * numpy.exp(-(1./(2.*sigma**2.))*(x[:] - 27.21138602*qc.mo_spec[selected_MO[i_mos]-1]['energy'])**2.)
   
  #--------------------------------------------
  # Print MO-coefficients
  #    for each atom a different file is created... this may be changed in the future
  #    one file for all atoms is written out
  #---
  
  #--------------------------------------------
  # Discrete values are written to files here
  #---
  bla = time.time()
  make_sure_path_exists('./pldos_data')
  for na in range(0,no_of_atoms+1):
    if (na < no_of_atoms):
      filestream = open('./pldos_data/'+str(na)+'.dat', "w+")
      print ("#     E / a.u.         E / eV               c-s             c-p             c-d            c-f           c-g", file=filestream)
      for mo in range(0,len(selected_MO)):
        print ("%16.8f" % (qc.mo_spec[selected_MO[mo]-1]['energy']), end="", file=filestream)
        print ("  %16.8f" % (qc.mo_spec[selected_MO[mo]-1]['energy']*27.21138602), end="", file=filestream)
        print ("  %13.8f" % coeff_array[mo,na,0], end="", file=filestream)
        print ("  %13.8f" % coeff_array[mo,na,1], end="", file=filestream)
        print ("  %13.8f" % coeff_array[mo,na,2], end="", file=filestream)
        print ("  %13.8f" % coeff_array[mo,na,3], end="", file=filestream)
        print ("  %13.8f" % coeff_array[mo,na,4], file=filestream)
    elif (na == no_of_atoms):
      filestream = open('./pldos_data/all_atoms.dat', "w+")
      print ("#     E / a.u.         E / eV               c-s             c-p             c-d            c-f           c-g", file=filestream)
      for mo in range(0,len(selected_MO)):
        print ("%16.8f" % (qc.mo_spec[selected_MO[mo]-1]['energy']), end="", file=filestream)
        print ("  %16.8f" % (qc.mo_spec[selected_MO[mo]-1]['energy']*27.21138602), end="", file=filestream)
        print ("  %13.8f" % numpy.sum(coeff_array[mo,:,0]), end="", file=filestream)
        print ("  %13.8f" % numpy.sum(coeff_array[mo,:,1]), end="", file=filestream)
        print ("  %13.8f" % numpy.sum(coeff_array[mo,:,2]), end="", file=filestream)
        print ("  %13.8f" % numpy.sum(coeff_array[mo,:,3]), end="", file=filestream)
        print ("  %13.8f" % numpy.sum(coeff_array[mo,:,4]), file=filestream)


  #--------------------------------------------
  # Smeared values are written to files here
  #---
  
  def writeOutSmearedCoef(output_filename, np_out):
    headerstring = ' Energy [eV]      C-s         C-p         C-d         C-f         C-g'
    numpy.savetxt(output_filename, np_out, fmt='%16.8f  %.8f  %.8f  %.8f  %.8f  %.8f', header=headerstring)
  
  if not spin_polarized:
    for na in range(0,no_of_atoms):
      writeOutSmearedCoef('./pldos_data/'+str(na)+'.smeared.dat', numpy.column_stack((x-E_fermi,gx[:,na,:])))
    writeOutSmearedCoef('./pldos_data/all_atoms.smeared.dat', numpy.column_stack((x-E_fermi,gx_all[:,:])))
      
  elif spin_polarized:
    for na in range(0,no_of_atoms):
      writeOutSmearedCoef('./pldos_data/'+str(na)+'.smeared.a.dat', numpy.column_stack((x-E_fermi,gx_alpha[:,na,:])))
      writeOutSmearedCoef('./pldos_data/'+str(na)+'.smeared.b.dat', numpy.column_stack((x-E_fermi,gx_beta[:,na,:])))
    writeOutSmearedCoef('./pldos_data/all_atoms.smeared.a.dat', numpy.column_stack((x-E_fermi,gx_all_a[:,:])))
    writeOutSmearedCoef('./pldos_data/all_atoms.smeared.b.dat', numpy.column_stack((x-E_fermi,gx_all_b[:,:])))
         
  #-------------------------------------------
  # Sharpo stops if it was only to do the MO-c
  #-------------------------------------------
  if (print_coeff_only == True):
    print ("Elapsed time for coefficients: %.2f s" % (time.time()-start_time))
    print ()
    sys.exit()
  else:
    print ("Sweet. Done with coefficients.")