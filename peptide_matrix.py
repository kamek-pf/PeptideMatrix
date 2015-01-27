import os.path
import math

class PeptideMatrix:

	def __init__(self, input = ''):
		super(PeptideMatrix, self).__init__()
		self.hasFile = False
		self.N = 0
		self.M = 0
		self.peptides = 0
		self.setFile(input)

	# Initialize the matrix
	# If a size is passed, the matrix is initialized to that size and filled with the singleton 'None'
	# Else, if a file path is set, the matrix is initialized with the file's content
	def init(self, N = 0, M = 0):
		if N > 0 and M > 0:
			self.N = N
			self.M = M
			self.matrix = [[None for i in range(self.M)] for i in range(self.N)]
		elif self.hasFile:
			self.computeSize()
			self.matrix = [[None for i in range(self.M)] for i in range(self.N)]
			self.fillFromFile()
		else:
			print('You need to either set a size or pass a file to parse')

	# Initialize a matrix with the best possible size
	def autoInit(self, peptides):
		if peptides > 0:
			self.peptides = peptides
		else:
			self.peptides = 121 
		size = self.guessBestSize(self.peptides)
		self.init(size[0], size[1])

	# Set the input file to read from
	def setFile(self, input):
		if input != '':
			if not os.path.exists(os.path.abspath(input)):
				print('The file', input, 'does not exist')
			else:
				self.input = input
				self.hasFile = True

	# Read the input file to compute matrix size
	def computeSize(self):	
		if self.hasFile:
			self.N = 0
			self.M = 0
			with open(self.input) as f: # open file
				for line in f:	# for each line
					# number columns
					if self.M < len(line.split(',')):
						self.M = len(line.split(','))
					# number of lines
					self.N += 1
			f.close()
			print('Matrix read : ', self.N, 'x', self.M)
		else:
			print('No file to read from')

	# Fill the matrix with reference numbers
	# For a NxN matrix, 1st line will be 1, 2, 3, ..., N
	# 2nd line will be N+1, N+2, N+3, ..., N+N
	def fillWithReference(self, max):
		value = 1
		for i in range(self.N):
			for j in range(self.M):
				if value <= max:
					self.matrix[i][j] = value
					value += 1

	# Read the input file and fill the matrix with its values
	def fillFromFile(self):
		i = 0
		j = 0
		with open(self.input) as f: 
			for line in f:	
				for value in map(int, line.split(',')):
					self.matrix[i][j] = value
					self.peptides += 1
					j += 1
					if j == self.M:
						j = 0
						i += 1
		f.close()

	# Display the matrix
	def output(self):
		lineIndicator = self.M + 1
		columnIndicator = [(i + 1) for i in range(self.M)]
		print('\t', end = '')
		for i in range(self.M):
			print('PL', columnIndicator[i], '\t', end = '')
		print('')

		for i in range(self.N):
			print('PL', lineIndicator, '\t', end = '')
			lineIndicator += 1
			for j in range(self.M):
				if self.matrix[i][j] is not None:
					print(self.matrix[i][j], '\t', end = '')
			print('')

	# Display the matrix line by line and column by column
	def align(self):
		# Print columns :
		pool = 1
		for column in range(self.M):
			label = "POOL " + str(pool) + " : "
			print(label, end=''),
			for line in range(self.N):
				if self.matrix[line][column] is not None:
					print(self.matrix[line][column], ' ' ,end='')
			print('')
			pool += 1

		# Print lines
		for line in range(self.N):
			label = "POOL " + str(pool) + " : "
			print(label, end=''),
			for column in range(self.M):
				if self.matrix[line][column] is not None:
					print(self.matrix[line][column], ' ' ,end='')
			print('')
			pool += 1

	# Find a value in the matrix, returns the indices as a tuple
	def find(self, value):
		for i in range(self.N):
			for j in range(self.M):
				if self.matrix[i][j] == value:
					return (i, j)
		return (None, None)

	# Create a new matrix using the values within the specified range only, then
	# match these values with the reference matrix 
	def match(self, min, max):
		if min < max:
			refMatrix = PeptideMatrix()
			refMatrix.init(self.N, self.M)
			refMatrix.fillWithReference(self.peptides)

			# Get indices for each values within the specified range
			indices = []
			for i in range(self.N):
				for j in range(self.M):
					if self.matrix[i][j] in range(min, max + 1):
						res = self.find(self.matrix[i][j])
						indices.append(res)

			# Get the best matrix size
			subMatrixSize = self.guessBestSize(len(indices))

			# Create a new submatrix
			subMatrix = PeptideMatrix()
			subMatrix.init(subMatrixSize[0], subMatrixSize[1])

			# Match the indices with the reference matrix and fill the submatrix
			position = 0
			for i in range(subMatrixSize[0]):
				for j in range(subMatrixSize[1]):
					if position < len(indices):
						line = indices[position][0]
						col = indices[position][1]
						subMatrix.matrix[i][j] = refMatrix.matrix[line][col]
					position += 1

			subMatrix.output()

	# Guess the best matrix size, what we want to avoid is a pool with only 1 peptide
	def guessBestSize(self, total):
		# Get the size of the smallest square matrix and work from there
		size = math.ceil(math.sqrt(total))
		N = size
		M = size
		emptySlots = N * M - total

		# If the matrix is full, we're good
		if emptySlots == 0:
			return (N,M)

		# If it doesn't fit perfectly, try to find out what size will have
		# less empty slots
		tmpN = N
		tmpM = M
		size = tmpN * tmpM
		for x in range(10):
			# Decrease lines and increase columns until we find the size that is
			# the closest to our number of peptides
			tmpN -= 1
			tmpM += 1
			newSize = tmpN * tmpM
			# If a better size is found, swap values
			# else, we already found the closest match, so we return it
			if newSize >= total and newSize < size:
				N = tmpN
				M = tmpM
				size = newSize
			else:
				return (N,M)

	# Read a matrix and recreate it with a better size
	def remap(self):
		if self.hasFile:
			data = []

			# Keep all values in a list
			for i in range(self.N):
				for j in range(self.M):
					data.append(self.matrix[i][j])

			# Reinitialize the matrix to the best possible size
			self.autoInit(self.peptides)

			# Reinject data
			currentPeptide = 0
			for i in range(self.N):
				for j in range(self.M):
					self.matrix[i][j] = data[currentPeptide]
					currentPeptide += 1

			print('Remapped matrix : ', self.N, 'x', self.M)
