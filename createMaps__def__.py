# createMaps__def__.py
import pyfits
from Nicola import *
from math import pi

#Retrieve dictionary
DropboxDirectory = os.getcwd().split('Dropbox')[0]
lib_path = os.path.abspath(DropboxDirectory+'Dropbox/PhD_Analysis/Library') 
sys.path.append(lib_path)
from galaxyParametersDictionary_v9 import *
# from astroquery.vizier import Vizier



# Vizier.ROW_LIMIT = -1
# cat_list = Vizier.get_catalogs('J/MNRAS/414/888/tableb1')

# cat = cat_list[0]

# cat.colnames
# # Galaxies = np.array(cat['Gal'])[np.where((np.array(cat['Rmax']) <= 2.0) & (np.array(cat['Rmax']) > 1.7))]
# Galaxies = ['NGC2699']

# b_a = {}
# for GalName in Galaxies:
#   b_a[GalName] = 1 - np.array(cat['epse'])[np.where(np.array(cat['Gal']) == GalName)]

# cat_list = Vizier.get_catalogs('J/MNRAS/413/813/atlas3d')
# cat = cat_list[0]
# cat.colnames

# Reff = {}
# vel0 = {}
# for GalName in Galaxies:
#   Reff[GalName] = 10**np.array(cat['log_Re_'])[np.where(np.array(cat['Gal']) == GalName)]
#   vel0[GalName] = np.array(cat['Vhel'])[np.where(np.array(cat['Gal']) == GalName)]

# exclude = {}
# Theta_Kriging_ATLAS = {}
# Theta_Kriging = {}
# for GalName in Galaxies:
#   exclude[GalName] = 'None'
#   Theta_Kriging_ATLAS[GalName] = 5
#   Theta_Kriging[GalName] = 15

# PA0 = {}

# x = np.loadtxt(DropboxDirectory+'Dropbox/PhD_Analysis/Data/ATLAS3D/Krajnovic11_TableD1.txt', comments = '#', unpack = True, usecols = [0], dtype = 'str')
# y = np.loadtxt(DropboxDirectory+'Dropbox/PhD_Analysis/Data/ATLAS3D/Krajnovic11_TableD1.txt', comments = '#', unpack = True, usecols = [1])
# for GalName in Galaxies:
#   PA0[GalName] = y[np.where(x == GalName)]

def findDell(RA, Dec, PA0, b_a, unfolded=False):
  angleRot = (numpy.pi/180.)*(PA0-90.)
  xrot, yrot = (RA *numpy.cos(angleRot) - Dec * numpy.sin(angleRot), 
                RA *numpy.sin(angleRot) + Dec * numpy.cos(angleRot))
  # 
  if unfolded:
    sign = (xrot/numpy.abs(xrot))
  else:
    sign = 1.
  Rell = (b_a*(xrot**2)+(yrot**2)/b_a)**(0.5)
  #
  return sign*Rell

def retrieve_krigPos(path, galname, excludeitem, offset=[0.,0.], 
					 keyword='VEL', clean=True):
  inputFits = pyfits.open(path)[1].data
  #
  RAgal = convAngCoord(CentreCoordinates[galname][0])[4]*15.
  print RAgal
  Decgal = convAngCoord(CentreCoordinates[galname][1])[4] 
  print Decgal
  Rell = findDell(inputFits['RA']+offset[0]/3600.-RAgal, 
  					inputFits['Dec']+offset[1]/3600.-Decgal, 
  					PA0[galname], b_a[galname])
  #
  # RA, RAg = np.radians(inputFits['RA']), np.radians(RAgal)
  # Dec, Decg = np.radians(inputFits['DEC']), np.radians(Decgal)
  # x = (np.sin(RA-RAg)*math.cos(Decg))*3600.+offset[0] 
  # y = (inputFits['Dec']-Decgal)*3600.+offset[1]

  # Not assuming small angle to identify whether it produces the same result (spoiler: it does)
  RA, RAg = (inputFits['RA']-RAgal)*np.pi/180, np.radians(RAgal)
  Dec, Decg = np.radians(inputFits['DEC']), np.radians(Decgal)
  x = ((np.sin(RA)*math.cos(Decg)))*180/np.pi*3600+offset[0] 
  y = ((np.sin(Dec)*np.cos(Decg)-np.cos(RA)*np.cos(Dec)*np.sin(Decg)))*180/np.pi*3600+offset[1]
  

  # Amendment made by Sabine to allow a kriging map of SN to be made 
  if keyword == 'SN':
    z, ez = inputFits[keyword], inputFits[keyword]
  else:
    z, ez = inputFits[keyword], inputFits['ERR'+keyword]

  #
  if clean:
    # additional condition implemented by Sabine so that a single point can be manually excluded from the kriging map. 
    selected = numpy.nonzero((inputFits['SKIMSFLAG'] == 'A') & (inputFits['NAME'] != excludeitem)) 
    x, y, z, ez = x[selected], y[selected], z[selected], ez[selected]
    # print inputFits['FILE'][selected]

  return [x, y, z, ez]
# Minor axis mask needs conversion of positions because namely the slits have 
# the same coordinates as the Major axis mask

# I'm applying the offset after the rotation
def retrieve_krigPos_Minor(path, galname, offset=[0.,0.], 
					 keyword='VEL'):
  inputFits = pyfits.open(path)[1].data
  #
  # Fix slit position in minor axis SuperSKiMS
  skypa1 = -276.70499351
  gal_PA0 = PA0[galname]
  mask_C_RA = '02:40:22.25'
  mask_C_Dec = '+39:02:48.2'
  #
  RAgal = convAngCoord(mask_C_RA)[4]*15.
  Decgal = convAngCoord(mask_C_Dec)[4]
  #
  deltaPA = mod(90.+skypa1-gal_PA0, 360.) #Difference between mask alignment and galaxy PA0
  angleNE = mod(90.-gal_PA0, 360.) #angle between galaxy major axis and East axis
  maskPA = -mod(angleNE+deltaPA, 360) #angle between the East-West axis and the mask alignment
  #
  distRA = inputFits['RA']-RA_c    #Distance from mask CentreCoordinates (in degrees)
  distDEC = inputFits['Dec']-Dec_c #Distance
  #
  angrot = maskPA*numpy.pi/180.
  realRA = (distRA*numpy.cos(angrot)-distDEC*numpy.sin(angrot))+RA_c+offset[0]/3600.   #Coordinates given the rotation of the mask
  realDEC = (distRA*numpy.sin(angrot)+distDEC*numpy.cos(angrot))+Dec_c+offset[1]/3600. #Coordinates given the rotation of the mask
  #
  Rell = findDell(realRA-RAgal, realDEC-Decgal, PA0[galname], b_a[galname])#,unfolded=True)
  #
  x, y = (realRA-RAgal)*3600., (realDEC-Decgal)*3600.
  z, ez = SS_Minor_input[keyword], SS_Minor_input['ERR'+keyword]
  #
  return [x, y, z, ez]




def drawLayout(galname='NGC1023'):
  fig = figure(0); clf(); ax = subplot(111); grid(True)
  ax.set_xlim([40.16,40.04]); ax.set_ylim([39.03,39.11])
  ax.set_aspect('equal')
  #Galaxy centre
  RAgal = convAngCoord(CentreCoordinates[galname][0])[4]*15.
  Decgal = convAngCoord(CentreCoordinates[galname][1])[4] 
  ax.scatter(RAgal, Decgal, c='r', marker='x', s=50)
  #
  #Galaxy axes
  major = ax.plot([numpy.sin(math.radians(PA0[galname]))*2+RAgal, numpy.sin(math.radians(PA0[galname]))*-2+RAgal], 
     [numpy.cos(math.radians(PA0[galname]))*2+Decgal, numpy.cos(math.radians(PA0[galname]))*-2+Decgal], 
     'r-')
  minor = ax.plot([numpy.cos(math.radians(PA0[galname]))*2+RAgal, numpy.cos(math.radians(PA0[galname]))*-2+RAgal], 
     [numpy.sin(math.radians(PA0[galname]))*-2+Decgal, numpy.sin(math.radians(PA0[galname]))*2+Decgal],
     'r-')
  #
  # Measured offset
  ax.scatter(RAgal-6./3600., Decgal+9./3600., marker='o', c='b')
  return ax



def drawMask_design(ax, angleMask=0., offsetX=0., offsetY=0):
  #fig = figure(0)
  path = '../JacobOutputs/MajorAxis/DEIMOS_pPXF_NGC1023_v0.fits'
  inputFits = pyfits.open(path)[1].data
  #
  angrot = math.radians(angleMask)
  distRA = inputFits['RA']
  distDec = inputFits['Dec']
  RA = (distRA*numpy.cos(angrot)-distDec*numpy.sin(angrot))+offsetX/3600.
  DEC = (distRA*numpy.sin(angrot)+distDec*numpy.cos(angrot))+offsetY/3600.
  #
  ax.scatter(RA, DEC, marker=(2,0,88.3), c='k')
  plt.draw()
  return True

