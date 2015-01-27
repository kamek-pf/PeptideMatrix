#!/usr/bin/python

from peptide_matrix import PeptideMatrix

# Read 
dataMatrix = PeptideMatrix('input')
dataMatrix.init()

# Remap the matrix to get a better size
dataMatrix.remap()

print('')

print('Data matrix :')
dataMatrix.output()

print('')

print('Matched matrix data')
dataMatrix.match(7750,7822)
