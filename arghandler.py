# ===========================================================
#
#   Module for sharpo
#   arghandler.py :: handels the command-line arguments
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

#  Last edited
#  20.01.2017 :: added E-Fermi shift option
#  15.09.2016 :: added #ofPoints argument
#  13.09.2016 :: completed
#
# ===========================================================


from __future__ import print_function
import sys
import os.path
import argparse

prog_description = '''
This is sharpo (Spherical HARmonics - Projected Orbitals). 
Copyright (C) 2016 Lukas Hammerschmidt. This program comes with ABSOLUTELY
NO WARRANTY. This is free software, and you are welcome to redistribute it
under certain conditions. Type "-l" for more details.
'''

prog_license = '''
sharpo
Lukas Hammerschmidt

MacDiarmid Institute for Advanced Materials and Nanotechnology,
Department of Physics, University of Auckland, Auckland, New Zealand

sharpo is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as 
published by the Free Software Foundation, either version 3 of 
the License, or any later version.

sharpo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public 
License along with sharpo.  If not, see <http://www.gnu.org/licenses/>.
'''

input_types = [ "molden", "gaussian.log", "aomix" ]
input_units = [ "A", "au" ]
global args

#------------------------------------------------
# Here command-line arguments are defined and set
# argparse is used for that
#------------------------------------------------
def get_parser(i_version, prog_info):
  
  parser = argparse.ArgumentParser(prog="sharpo",description=prog_description)

  parser.add_argument("-l", "--license", action="store_true", help="information about the license")
  parser.add_argument('-v', '--version', action='version', version='%(prog)s '+i_version)
  parser.add_argument("INPUT", type=str, nargs='?', help="Provide an input name.")
  parser.add_argument("-r", "--radius", type=float, default="5.", help="Integration cut-off radius in [units (default=[A])]. Default=5.")
  parser.add_argument('-c', '--center', nargs=3, type=float, default=[0., 0., 0.],
                      metavar=("x","y","z"), help='Sets the center for the SH expansion')
  parser.add_argument('-u', '--units', type=str, choices=input_units, default='A', help="Set length units to Ang or a.u. Default=[A]")
  parser.add_argument("-t", "--type", type=str, nargs='?', choices=input_types,
                    help="Provide the input type: " + ", ".join(input_types), 
                    metavar='TYPE', default="molden")
  parser.add_argument("-o", "--output", type=str, default='porb', help='Define output file name')
  parser.add_argument("-s", "--sigma", type=float, default="0.1", help="Sigma. Smearing factor for Gaussian smearing of projected states")
  parser.add_argument("-np", "--npoints", type=float, default="50.", help="#ofPoints/eV for the smeared output. (Default = 50p/eV)")
  parser.add_argument('-nef', "--no-efermi", dest='efermi_shift', action='store_false', help='Prevents Sharpo from shifting the orbital energies by E-E_HOMO (default=True).')
  parser.add_argument("--e-range", type=float, nargs = 2, default=[0., 0.], dest='e_range', metavar=("Emin", "Emax"), help='in [eV] - Lets you determine an energy range between which MOs are considered in the expansion. If 0. 0. is chosen, the whole energy range is scanned (default). The energy scale is set such that E_HOMO = 0 (E_F = 0) except if "--no-efermi" is set. Then the absolute energy values are taken.')
  parser.add_argument("--mo-energies-only", action = 'store_true', dest='moe_only', help='Prints a list of the MO energies and stops.')
  group = parser.add_mutually_exclusive_group()
  group.add_argument('--coeff', dest='print_coeff', action='store_true', help='Print MO-coefficients')
  #parser.add_argument('--coeff', dest='print_coeff', action='store_true', help='Print MO-coefficients')
  group.add_argument('--coeff-only', dest='print_coeff_only', action='store_true', help='Print MO-coefficients and stop program')
  parser.add_argument("-p", "--proc", type=int, help="Number of subprocess within orbkit",
                      metavar="NUM", default = 1)
  parser.add_argument("--info", action="store_true", 
                      help="Prints important information about the compatibility of sharpo")
  parser.set_defaults(print_coeff=False)

  args = parser.parse_args()

  #----------------------
  # license print
  if args.license:
    print(prog_license)
    sys.exit(1)
  
  #----------------------
  # Program info
  if args.info:
    print(prog_info)
    sys.exit(1)

  #----------------------
  # No input > print help
  if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
  
  #----------------------
  # if file is not a file
  if not os.path.isfile(args.INPUT):
    print ("INPUT provided is not a file or doesn't exist!")
    print
    parser.print_help()
    sys.exit(1)

  #-------------------------------
  # check emin, emax inputs
  if (args.e_range[0] > args.e_range[1]):
    tmpE = args.e_range[0]
    args.e_range[0] = args.e_range[1]
    args.e_range[1] = tmpE
  if (abs(args.e_range[0] - args.e_range[1]) < 0.1 and abs(args.e_range[0]) > 0.001 and abs(args.e_range[1]) > 0.001):
    raise argparse.ArgumentTypeError("Emin and Emax are too close to each other. Please choose larger energy range!")
    sys.exit(1)
  if abs(args.e_range[0]) < 0.001 and abs(args.e_range[1]) < 0.001:
    energy_range_is_set = False
  else:
    energy_range_is_set = True
 
  
  return_list = [args.INPUT, str(args.radius), args.type, args.proc,
                 args.center[0], args.center[1], args.center[2], args.output, 
                 args.units, args.sigma, args.npoints, args.print_coeff, 
                 args.print_coeff_only, args.efermi_shift, args.e_range[0], 
                 args.e_range[1], energy_range_is_set, args.moe_only]
  return (return_list)
  
