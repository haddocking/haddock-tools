from math import acos, pi, sqrt
import re
import string
import sys

class Dipole:

  def __init__(self, id, D, segid, Dxx, Dyy, Dzz):

    Dmin = min(min(Dxx, Dyy), Dzz)
    Dmax = max(max(Dxx, Dyy), Dzz)

    if (D < Dmin):

      print 'changing dipole id', id, segid, 'from', D, 'to', Dmin
      D = Dmin

    elif (D > Dmax):

      print 'changing dipole id', id, segid, 'from', D, 'to', Dmax
      D = Dmax

    self.id = id
    self.D = D
    self.segid = segid

    x = (D - Dzz) / (Dyy - Dzz) # yz plane
    y = (D - Dxx) / (Dzz - Dxx) # xz plane
    z = (D - Dyy) / (Dxx - Dyy) # xy plane

    #print 'x, y, z =', x, y, z

    if ((x < 0) or (x > 1)):

      s1 = sqrt(y)
      c1 = sqrt(1-y)
      s2 = sqrt(z)
      c2 = sqrt(1-z)

      self.excluded = 0
      self.circle1 = [ (c1, 0, s1), (c1, 0, -s1), (s2, c2, 0), (s2, -c2, 0) ] 
      self.circle2 = [ (-c1, 0, s1), (-c1, 0, -s1), (-s2, c2, 0), (-s2, -c2, 0) ] 
      
    elif ((y < 0) or (y > 1)):

      s1 = sqrt(z)
      c1 = sqrt(1-z)
      s2 = sqrt(x)
      c2 = sqrt(1-x)

      self.excluded = 1
      self.circle1 = [ (s1, c1, 0), (-s1, c1, 0), (0, s2, c2), (0, s2, -c2) ]
      self.circle2 = [ (s1, -c1, 0), (-s1, -c1, 0), (0, -s2, c2), (0, -s2, -c2) ]

    else: # z <= 0  or z >= 1

      s1 = sqrt(x)
      c1 = sqrt(1-x)
      s2 = sqrt(y)
      c2 = sqrt(1-y)

      self.excluded = 2
      self.circle1 = [ (0, s1, c1), (0, -s1, c1), (c2, 0, s2), (-c2, 0, s2) ]
      self.circle2 = [ (0, s1, -c1), (0, -s1, -c1), (c2, 0, -s2), (-c2, 0, -s2) ]

class DipolePair:

  def __init__(self, id1, id2, segid1, segid2, phimin1, phimax1, phimin2, phimax2):

    if (phimin2 < phimin1):

      (phimin1, phimin2) = (phimin2, phimin1)
      (phimax1, phimax2) = (phimax2, phimax1)

    if (phimin2 < phimax1): # one range

      phimax1 = 90
      phimin2 = 90

    included = phimax2 - phimin2 + phimax1 - phimin1

    self.id1 = id1
    self.id2 = id2
    self.segid1 = segid1
    self.segid2 = segid2
    self.phimin1 = phimin1
    self.phimax1 = phimax1
    self.phimin2 = phimin2
    self.phimax2 = phimax2
    self.phimid1 = 0.5 * (self.phimin1 + self.phimax1)
    self.phiwidth1 = 0.5 * (self.phimax1 - self.phimin1)
    self.phimid2 = 0.5 * (self.phimin2 + self.phimax2)
    self.phiwidth2 = 0.5 * (self.phimax2 - self.phimin2)

    self.excluded = 1 - included/180

  def write(self, file):

    file.write('assign (resid ')
    file.write(self.id1)
    file.write(' and name N and segid ')
    file.write(self.segid1)
    file.write(' ) (resid ')
    file.write(self.id1)
    file.write(' and name HN and segid ')
    file.write(self.segid1)
    file.write(' ) (resid ')
    file.write(self.id2)
    file.write(' and name N and segid ')
    file.write(self.segid2)
    file.write(' ) (resid ')
    file.write(self.id2)
    file.write(' and name HN and segid ')
    file.write(self.segid2)
    file.write(' ) ')
    file.write('%2.1f' % self.phimid1)
    file.write(' ')
    file.write('%2.1f' % self.phiwidth1)
    file.write(' ')
    file.write('%2.1f' % self.phimid2)
    file.write(' ')
    file.write('%2.1f' % self.phiwidth2)
    file.write(' ! excluded ')
    file.write('%4.3f' % self.excluded)
    file.write('\n')

def compare_dipoles(d1, d2):

  dmin = 2
  dmax = - 2

  for point1 in d1.circle1:

    (x1, y1, z1) = point1
    #print '1:', x1, y1, z1

    for point2 in d2.circle1:

      (x2, y2, z2) = point2
      #print '  2:', x2, y2, z2

      d = x1*x2 + y1*y2 + z1*z2
    
      dmin = min(d, dmin)
      dmax = max(d, dmax)
      #print '  3:', d, dmin, dmax

  phimin1 = acos(dmax) * 180 / pi
  phimax1 = acos(dmin) * 180 / pi

  dmin = 2
  dmax = - 2

  for point1 in d1.circle2:

    (x1, y1, z1) = point1
    #print '4:', x1, y1, z1

    for point2 in d2.circle1:

      (x2, y2, z2) = point2
      #print '  5:', x2, y2, z2

      d = x1*x2 + y1*y2 + z1*z2
    
      dmin = min(d, dmin)
      dmax = max(d, dmax)
      #print '  6:', d, dmin, dmax

  phimin2 = acos(dmax) * 180 / pi
  phimax2 = acos(dmin) * 180 / pi

  dipole_pair = DipolePair(d1.id, d2.id, d1.segid, d2.segid, phimin1, phimax1, phimin2, phimax2)

  return dipole_pair

def write_dipole_pairs(file, dipole_pairs):

  for dipole_pair in dipole_pairs:

    dipole_pair.write(file)

def read_dipoles(input_file, Dxx, Dyy, Dzz):

  file = open(input_file, 'r')

  for i in range(1): # header information
    line = file.readline()

  lines = file.readlines()

  file.close()

  dipoles = []
  #count = [ 0, 0, 0 ]

  for line in lines:

    line = string.strip(line)

    if (line):

      #[ id, junk1, junk2, value ] = re.split('\t', line)
      [ id, value, segid ] = re.split('\t', line)

      dipole = Dipole(id, string.atof(value), segid, Dxx, Dyy, Dzz)
      dipoles.append(dipole)

      #count[dipole.excluded] = count[dipole.excluded] + 1

  print 'found', len(dipoles), 'dipoles'
  #print 'exclusion cases:', count[0], count[1], count[2]

  return dipoles

def compare_exclusions(dipole_pair1, dipole_pair2):

  if (dipole_pair1.excluded > dipole_pair2.excluded):

    return -1

  elif (dipole_pair1.excluded < dipole_pair2.excluded):

    return 1

  else:

    return 0

def process_dipoles(input_file, output_file, Dxx, Dyy, Dzz,
                    correct_trace = 1, min_excluded = 0, max_pairs = 0):

  if (correct_trace):

    t = - (Dxx + Dyy + Dzz) / 3

    if (t != 0):

      print 'adding', t, 'to Dxx, Dyy, Dzz'

      Dxx = Dxx + t
      Dyy = Dyy + t
      Dzz = Dzz + t

  dipoles = read_dipoles(input_file, Dxx, Dyy, Dzz)
  file = open(output_file, 'w')

  dipole_pairs = []

  for n1 in range(len(dipoles)-1):

    print 'working on dipole id', dipoles[n1].id

    for n2 in range(n1+1, len(dipoles)):

      dipole_pair = compare_dipoles(dipoles[n1], dipoles[n2])

      if (dipole_pair.excluded >= min_excluded):

        dipole_pairs.append(dipole_pair)

  dipole_pairs.sort(compare_exclusions)

  if (max_pairs):

    dipole_pairs = dipole_pairs[:max_pairs]

  write_dipole_pairs(file, dipole_pairs)

  file.close()

  print 'output', len(dipole_pairs), 'dipole pairs'

if (__name__ == '__main__'):

  if (len(sys.argv) != 6):

    print 'args required: input file, output file, Dxx, Dyy, Dzz'

    sys.exit(1);

  input_file = sys.argv[1]
  output_file = sys.argv[2]
  Dxx = string.atof(sys.argv[3])
  Dyy = string.atof(sys.argv[4])
  Dzz = string.atof(sys.argv[5])

  process_dipoles(input_file, output_file, Dxx, Dyy, Dzz)
