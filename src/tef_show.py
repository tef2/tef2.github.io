#!/usr/bin/env python


"""
 This file is part of TEF.

    TEF is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TEF is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TEF.  If not, see <http://www.gnu.org/licenses/>.


(C) 2018 Stratmann D, Pathmanathan JS, Postic G, Rey J, Chomilier J
Contact: dirk.stratmann@sorbonne-universite.fr
"""

import string
import colorsys

def color_obj(nobj, rainbow=0, jmol=0):

	"""

AUTHOR 

        Gareth Stockwell

USAGE

        color_obj(rainbow=0)

        This function colours each object currently in the PyMOL heirarchy
        with a different colour.  Colours used are either the 22 named
        colours used by PyMOL (in which case the 23rd object, if it exists,
        gets the same colour as the first), or are the colours of the rainbow

        """

	# Process arguments
	rainbow = int(rainbow)
	jmol = int(jmol)
	
	colors = []
	if rainbow:

		#print "\nColouring objects as rainbow\n"

		# Create colours starting at blue(240) to red(0), using intervals
		# of 240/(nobj-1)
		for j in range(nobj):
			hsv = ((240-float(j)*240/(nobj-1))/float(360), 1, 1)
			# Convert to RGB
			rgb = colorsys.hsv_to_rgb(hsv[0],hsv[1],hsv[2])
			# Define the new colour
			colors.append(rgb)

	else:

		#print "\nColouring objects using PyMOL defined colours\n"

		# List of available colours
		colours = []
		# for PyMol only:
		coloursPyMol = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
				   'salmon', 'lime', 'pink', 'slate', 'magenta', 'orange', 'marine', \
				   'olive', 'purple', 'teal', 'forest', 'firebrick', 'chocolate',    \
				   'wheat']
		# for Jmol only:
		coloursJmol = ['red', 'green', 'blue', 'yellow', 'violet', 'cyan',    \
				   'salmon', 'lime', 'pink', 'mediumslateblue', 'magenta', 'orange', 'dodgerblue', \
				   'olive', 'purple', 'teal', 'forestgreen', 'firebrick', 'chocolate',    \
				   'wheat']
		coloursJmolHex = ['FF0000', '008000', '0000FF', 'FFFF00', 'EE82EE', '00FFFF',    \
				   'FA8072', '00FF00', 'FFC0CB', '7B68EE', 'FF00FF', 'FFA500', '1E90FF', \
				   '808000', '800080', '008080', '228B22', 'B22222', 'D2691E',    \
				   'F5DEB3']
				   
		if jmol == 1:
			colours = coloursJmol
		elif jmol == 2:
			colours = coloursJmolHex
		else:
			colours = coloursPyMol

		ncolours = len(colours)

		# Loop over objects
		j=0
		for i in range(nobj):
			colors.append(colours[j])
			j = j+1
			if(j == ncolours):
				j = 0
	return colors



