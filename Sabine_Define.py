# import pyfits, 
import math, os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

DropboxDirectory = os.getcwd().split('Dropbox')[0]
lib_path = os.path.abspath(DropboxDirectory+'Dropbox/PhD_Analysis/Library') 
sys.path.append(lib_path)
from galaxyParametersDictionary_v9 import *
# from JAM_Sabine_def import *
from JAM_Dictionary import *


def parameterExtractor(inputDict, name):
        '''
        This function extracts the value and associated uncertainties from a specified parameter in a 
        chainconsumer dictionary. 
        '''
        value = inputDict[name][1]
        # try:
        lower = inputDict[name][1] - inputDict[name][0]
        upper = inputDict[name][2] - inputDict[name][1]
        # except:
        #         lower, upper = 0.0, 0.0
        # print 'Parameter extraction successful for', name
        return value, lower, upper


# extract a column of the data input as a separate array from a fits file. 
def column(columnnumber, array):
	Output=[]
	for i in range(len(array)):
		Output.append(array[i][columnnumber])
	return (Output)

def column_txt(columnnumber, filename):
	Output = np.loadtxt(filename, usecols = [columnnumber])
	return (Output)

def convert(RA, DEC):
	return(RA*np.pi/180, DEC*np.pi/180)

def radian(value):
	return(value*np.pi/180)

def degree(value):
	return(value*180/np.pi)

def radius(RA, Dec, GalRA, GalDec, phi, ellipticity):
	Q=1-ellipticity
	# We want the RA and Dec in radians.
	RA_r, Dec_r=convert(RA, Dec)
	# now we calculate the radius
	GalRA_r, GalDec_r=convert(GalRA, GalDec)
	rag=(RA_r-GalRA_r)*np.cos(GalDec_r)
	decg=(Dec_r-GalDec_r)
	PA=math.atan2(rag, decg)
	# Calculating the projected RA and Dec
	DEC_p=decg*np.cos(phi)+rag*np.sin(phi)
	RA_p=rag*np.cos(phi)-decg*np.sin(phi)
	radius=degree(np.sqrt((DEC_p**2)/Q+(Q*RA_p**2))) # whether or not to multiply by 3600 is given in this line
	return (radius)

def radiusCircular(RA, Dec, GalRA, GalDec, phi, ellipticity):
	Q=1-ellipticity
	# We want the RA and Dec in radians.
	RA_r, Dec_r=convert(RA, Dec)
	# now we calculate the radius
	GalRA_r, GalDec_r=convert(GalRA, GalDec)
	rag=(RA_r-GalRA_r)*np.cos(GalDec_r)
	decg=(Dec_r-GalDec_r)
	radius=degree(np.sqrt((decg**2)+(rag**2))) # whether or not to multiply by 3600 is given in this line
	return (radius)

def radiusArray(RA, Dec, phi, ellipticity):
    Q=1-ellipticity
    # We want the RA and Dec in radians.
    RA_r, Dec_r=convert(RA, Dec)
    DEC_p=Dec_r*np.cos(phi)+RA_r*np.sin(phi)
    RA_p=RA_r*np.cos(phi)-Dec_r*np.sin(phi)
    radius=degree(np.sqrt((DEC_p**2)/Q+(Q*RA_p**2))) # whether or not to multiply by 3600 is given in this line
    return (radius)

def coordinateRotation(RA, Dec, GalRA, GalDec, phi):
	# We want the RA and Dec in radians.
	RA_r, Dec_r=convert(RA, Dec)
	GalRA_r, GalDec_r=convert(GalRA, GalDec)
	rag=(RA_r-GalRA_r)*np.cos(GalDec_r)
	decg=(Dec_r-GalDec_r)
	# Calculating the projected RA and Dec
	DEC_p=decg*np.cos(phi)+rag*np.sin(phi)
	RA_p=rag*np.cos(phi)-decg*np.sin(phi)
	return (degree(RA_p), degree(DEC_p))

def positionAngle(RA, Dec, GalRA, GalDec):
	x=radian(GalRA-RA)*np.cos(radian(GalDec))
	y=radian(GalDec-Dec)
	theta=np.arctan2(-x, -y)
        if len(RA) > 1: # this function works on an array and on an individual value
                sel = np.where(theta<0)
                theta[sel] = theta[sel] + 2*np.pi
        elif theta < 0:
                theta = theta + 2*np.pi 
	return (degree(theta))

def maximumRadiusCoverageAnnulus(error, Radius, areaPixel, ellipticity):
	Radius = np.array(Radius)
	maximumRadius = 0
	Flag = True
	for radius in range(int(min(Radius))+1, int(max(Radius)), 2):
		Annulus_Range = np.where( (Radius >= (radius-error)) & (Radius <= (radius+error)) )[0]
		if (len(Annulus_Range) > 0) and Flag:

			OuterRadius = (radius+error)
			InnerRadius = (radius-error)
			areaTotal = np.pi * ((OuterRadius)**2 - (InnerRadius)**2) # area of annulus on Sky
			areaAnnulus = (len(Annulus_Range)*areaPixel)# area of pixels in annulus
			
			coverage = (areaAnnulus/areaTotal) * 100. 
			if coverage < 85.:
				Flag = False
			else:
				maximumRadius = radius
		# else:
		# 	print 'no pixels at ', radius
	return maximumRadius

# def maximumRadiusCoverageAnnulusAzimuth(error, Radius, areaPixel, ellipticity, X, Y):
# 	# doing this a different way. Trying it by assessing the degree of azimuthal coverage
# 	Radius = np.array(Radius)
# 	maximumRadius = 0
# 	Azimuth = np.arctan2(Y,  X) * 180/np.pi
# 	Flag = True
# 	MinimumGap = 0.15*360
# 	for radius in range(int(min(Radius))+1, int(max(Radius)), 2):
# 		AzimuthalGap = 0
# 		for ii in np.arange(0, 360, 2):
# 			Selection = np.where( (Azimuth <= ii+1) & (Azimuth >= ii-1) & (Radius >= (radius-error)) & (Radius <= (radius+error)) )
# 			if len(Selection[0]) == 0:
# 				AzimuthalGap += 2
# 			# else:
# 			# 	print Azimuth[Selection]
# 		if AzimuthalGap >= MinimumGap:
# 			Flag = False
# 		else:
# 			maximumRadius = radius

# 	return maximumRadius



def maximumRadiusSigmaResolution(error, Radius, Sigma, SigmaResolution):
	Sigma = np.array(Sigma)
	Radius = np.array(Radius)
	# # first I need to identify the whole possible sample at each annulus. 
	InaccurateRadius = []
	for radius in range(int(min(Radius))+1, int(max(Radius)), 2):
		Annulus_Range = []
		for i in range(len(Radius)):
			if Radius[i] > (radius-error) and Radius[i] < (radius+error):
				Annulus_Range.append(i)
		try:
			Annulus_Range = np.array(Annulus_Range)
			Sigma_Annulus = Sigma[Annulus_Range]
			Radius_Annulus = Radius[Annulus_Range]
			
			# calculate the fraction of points below the sigma resolution
			# print len(np.where(Sigma_Annulus < SigmaResolution)[0])
			# print len(Sigma_Annulus)
			Sigma_Fraction = len(np.where(Sigma_Annulus < SigmaResolution)[0]) / len(Sigma_Annulus)
			if Sigma_Fraction > 0.2:
				InaccurateRadius.append(radius)
		except:
			print 'error calculating coverage at radius = ', radius
	if len(InaccurateRadius) == 0:
		maximumRadius = Radius[-1]
	else:
		maximumRadius = InaccurateRadius[0]
	return maximumRadius

def maximumRadiusCoverageEllipse(error, Radius, Vel, Sigma, areaPixel):
	Radius = np.array(Radius)
	# print Radius[0]
	# # first I need to identify the whole possible sample at each annulus. 
	for radius in range(int(min(Radius))+1, int(max(Radius)), 2):
		try:
			Selection = np.where(Radius <= radius)
			
			# we only wish to continue the analysis if the bin has a good radial coverage
			areaTotal = np.pi * (radius)**2
			coverage = len(Selection)/(areaTotal/areaPixel) * 100.
			if coverage >= 85.:
				maximumRadius = radius
		except:
			print 'error calculating coverage at radius = ', radius
	return maximumRadius

def parameterBootstrap(radius, error, Radius, Vel, VelErr, VelDisp, VelDispErr):
	# radius - the desired radius at which the parameter is measured
	# error - the spread in radius over which SKiMS points are used to measure the lambda_R parameter
	# Radius, Vel, Sigma - radius, velocity and velocity dispersion of points. 

	# # first I need to identify the whole possible sample at each annulus. 
	Annulus_Range = []
	for i in range(len(Radius)):
		if Radius[i] > (radius-error) and Radius[i] < (radius+error):
			Annulus_Range.append(i)
	Annulus_Range = np.array(Annulus_Range)
	Vel_Annulus = Vel[Annulus_Range]
	VelErr_Annulus = VelErr[Annulus_Range]
	VelDisp_Annulus = VelDisp[Annulus_Range]
	VelDispErr_Annulus = VelDispErr[Annulus_Range]
	Radius_Annulus = Radius[Annulus_Range]

	# now I want to randomly sample 50 points from this array. Let repetitions occur
	# Note that this may cause smaller errors for annuli in which there are less data points, because they 
	# will be doubly selected more frequently...
	Selection = []
	Lambda = []
	# print len(Annulus_Range)
	for jj in range(100):
		for ii in range(50):
			Selection.append(np.random.randint(0, len(Annulus_Range)-1))
	
		Vel_Selection = Vel_Annulus[Selection]
		VelErr_Selection = VelErr_Annulus[Selection]
		VelDisp_Selection = VelDisp_Annulus[Selection]
		VelDispErr_Selection = VelDispErr_Annulus[Selection]
		Radius_Selection = Radius_Annulus[Selection]
	
		Lambda.append(np.sum(Radius_Selection*np.fabs(Vel_Selection))/np.sum(Radius_Selection*np.sqrt(Vel_Selection**2 + VelDisp_Selection**2)))

	# X = np.sqrt(Vel_Annulus**2 + VelDisp_Annulus**2)
	# # X_Err = np.sqrt((Vel_Annulus**2 * VelErr_Annulus**2 + VelDisp_Annulus**2 * VelDispErr_Annulus**2) / (Vel_Annulus**2 + VelDisp_Annulus**2))
	# X_Err = ((VelErr_Annulus/Vel_Annulus) + (VelDispErr_Annulus/VelDisp_Annulus)) * X
	# B = np.sum(Radius_Annulus * X)
	# B_Err = np.sqrt(np.sum((Radius_Annulus * X_Err)**2))
	# A = np.sum(Radius_Annulus * Vel_Annulus)
	# A_Err = np.sqrt(np.sum((Radius_Annulus * VelErr_Annulus)**2))
	# Total_Lambda_Error = np.array(Lambda).mean() * np.sqrt((A_Err/A)**2 + (B_Err/B)**2)

	return (np.array(Lambda).mean(), np.array(Lambda).std())#Total_Lambda_Error)

def parameterBootstrapFlux(radius, error, Radius, Vel, Sigma, Flux):
	# radius - the desired radius at which the parameter is measured
	# error - the spread in radius over which SKiMS points are used to measure the lambda_R parameter
	# Radius, Vel, Sigma - radius, velocity and velocity dispersion of points. 

	Vel = np.array(Vel)
	Sigma = np.array(Sigma)
	Radius = np.array(Radius)
	Flux = np.array(Flux)
	# # first I need to identify the whole possible sample at each annulus. 
	Annulus_Range = []
	for i in range(len(Radius)):
		if Radius[i] > (radius-error) and Radius[i] < (radius+error):
			Annulus_Range.append(i)
	Annulus_Range = np.array(Annulus_Range)
	Vel_Annulus = Vel[Annulus_Range]
	Sigma_Annulus = Sigma[Annulus_Range]
	Radius_Annulus = Radius[Annulus_Range]
	Flux_Annulus = Flux[Annulus_Range]

	# now I want to randomly sample 50 points from this array. Let repetitions occur
	# Note that this may cause smaller errors for annuli in which there are less data points, because they 
	# will be doubly selected more frequently...
	Selection = []
	Lambda = []
	for jj in range(100):
		for ii in range(50):
			# print len(Annulus_Range)
			Selection.append(np.random.randint(0, len(Annulus_Range)-1))
	
		Vel_Selection = Vel_Annulus[Selection]
		Sigma_Selection = Sigma_Annulus[Selection]
		Radius_Selection = Radius_Annulus[Selection]
		Flux_Selection = Flux_Annulus[Selection]
	
		Lambda.append(np.sum(Flux_Selection*Radius_Selection*np.fabs(Vel_Selection))/np.sum(Flux_Selection*Radius_Selection*np.sqrt(Vel_Selection**2 + Sigma_Selection**2)))

	return (np.array(Lambda).mean(), (np.array(Lambda).max()-np.array(Lambda).min())/2)


def FitsToTxt(fitsfile):
# converting a fits file into a txt file. 
	data=pyfits.getdata(fitsfile, 0)
	file=open(fitsfile.replace('.fits', '.txt'), 'w')
	for ii in range(len(data)):
		for i in range(37):
			file.write(str(column(i, data)[ii])+'\t')
		file.write('\n')

def onclick(event):
  global coords
  try: length = len(coords)
  except: coords = []
  if event.button==1:
    print 'Point selected to remove: xdata=%f, ydata=%f'%(event.xdata, event.ydata)
    plt.scatter(event.xdata, event.ydata, s=60, marker = 's', c='r')
    plt.draw()
    coords.append((event.xdata, event.ydata))


def PlotPointSelector(x, y, z, Sel, namegal):
  # making an interactive plot that enables points to be selected and then removed. 
  print 'Select points to be removed by clicking on them. '
  print 'Close the plot when selection is complete. '

  fig=plt.figure()
  ax=fig.add_subplot(111)
  cax=ax.scatter(numpy.array(x)[Sel], numpy.array(y)[Sel], c=numpy.array(z)[Sel], s=50, marker='o', \
    norm=mpl.colors.Normalize(vmin=np.min(numpy.array(z)[Sel]), vmax=np.max(numpy.array(z)[Sel])))
  ax.set_xlabel('x')
  ax.set_ylabel('y')
  ax.set_xlim(max(numpy.array(x)[Sel])+5, min(numpy.array(x)[Sel])-5)
  cb1=fig.colorbar(cax)
  cb1.set_label('z')
  
  plt.title(str(namegal))
  cid = fig.canvas.mpl_connect('button_press_event', onclick)
  plt.show()
  plt.close()

  # now I need to identify which array indices these selected coordinates relate to, 
  # by identifying the closest in x/y space
  try:
    print coords
    Selection = True
  except:
    Selection = False

  if Selection:
    Rejection = []
    for kk in range(len(coords)):
      x_DIFF, y_DIFF = 1000, 1000
      for jj in range(len(x)):
        if (abs(x[jj]-coords[kk][0]) < x_DIFF) and (abs(y[jj]-coords[kk][1]) < y_DIFF):
          Rejection_element = jj
          x_DIFF = abs(x[jj]-coords[kk][0])
          y_DIFF = abs(y[jj]-coords[kk][1])
      Rejection.append(Rejection_element)
    Rejection = numpy.array(Rejection)
  
    Sel = numpy.array(Sel[0])
    Sel = numpy.array(Sel[numpy.where(-(numpy.in1d(Sel, Rejection)-1))])
  
    fig=plt.figure()
    ax=fig.add_subplot(111)
    cax=ax.scatter(numpy.array(x)[Sel], numpy.array(y)[Sel], c=numpy.array(z)[Sel], s=50, marker='o', \
      norm=mpl.colors.Normalize(vmin=np.min(numpy.array(z)[Sel]), vmax=np.max(numpy.array(z)[Sel])))
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim(max(numpy.array(x)[Sel])+5, min(numpy.array(x)[Sel])-5)
    cb1=fig.colorbar(cax)
    cb1.set_label('z')
    
    plt.title(str(namegal))
    plt.show()
    plt.close()
  return Sel

def LocalMinimaSelector(x, Vel, namegal):
  # making an interactive plot that enables points to be selected and then removed. 
  print 'Select the local maximum and minimum by clicking on them. '
  print 'Close the plot when selection is complete. '

  fig=plt.figure()
  ax=fig.add_subplot(111)
  ax.scatter(x, Vel, s=2, marker='o')
  ax.set_xlabel('x')
  ax.set_ylabel('Vel')

  plt.title(str(namegal))
  cid = fig.canvas.mpl_connect('button_press_event', onclick)
  plt.show()
  plt.close()

  # print coords[0][0], coords[1][0]
  if coords[-2][1] < coords[-1][1]:
  	Minimum = coords[-2][1] 
  	Maximum = coords[-1][1]
  else:
  	Minimum = coords[-1][1] 
  	Maximum = coords[-2][1]
  
  return (Minimum, Maximum)

def CoordinateSelection(SLUGGS_RA, SLUGGS_Dec, ATLAS_RA, ATLAS_Dec):
	Indices = []
	for ii in range(len(SLUGGS_RA)):
		Duplicate = False
		for jj in range(len(ATLAS_RA)):
			if (abs(SLUGGS_RA[ii] - ATLAS_RA[jj]) < 1) and (abs(SLUGGS_Dec[ii] - ATLAS_Dec[jj]) < 1) and (Duplicate == False):
				Indices.append(ii)
				Duplicate = True
	return Indices

def pixelArea(RA, Dec, Radius):
	# identify the central point as being that which has the smallest radius
	CentralPixelIndex = np.argsort(Radius)[0]
	centralRA, centralDec = RA[CentralPixelIndex], Dec[CentralPixelIndex]
	# now I need to identify the pixels that are the neighbouring pixels. 
	# first identify the closes Dec (but not the same)
	Diff = 100
	for ii in range(len(Dec)):
		if (abs(Dec[ii] - centralDec) < Diff) and (Dec[ii] != centralDec):
			MinimumDec = Dec[ii]
			Diff = abs(Dec[ii] - centralDec)
	Diff = 100
	for ii in range(len(RA)):
		if (abs(RA[ii] - centralRA) < Diff) and (RA[ii] != centralRA):
			MinimumRA = RA[ii]
			Diff = abs(RA[ii] - centralRA)

	lengthPixel = abs(MinimumRA - centralRA)
	heightPixel = abs(MinimumDec - centralDec)
	
	areaPixel = lengthPixel * heightPixel
	print lengthPixel, heightPixel
	return areaPixel


def TotalAngularMomentum(Radius, Vel, VelDisp, GalRA, GalDec, phi, Ellipticity, maximumRadius):
	Total_Lambda = 0
	# print len(Radius)
	Sel = np.where(Radius < maximumRadius)
	Total_Lambda = Total_Lambda + np.sum(Radius[Sel]*np.fabs(Vel[Sel]))/np.sum(Radius[Sel]*np.sqrt(Vel[Sel]**2 + VelDisp[Sel]**2))
	return Total_Lambda

def TotalAngularMomentumError(GalName, Radius, Vel, VelErr, VelDisp, VelDispErr, maximumRadius):
	Sel = np.where(Radius < maximumRadius)
	X = np.sqrt(Vel[Sel]**2 + VelDisp[Sel]**2)
	X_Err = np.sqrt((Vel[Sel]**2 * VelErr[Sel]**2 + VelDisp[Sel]**2 * VelDispErr[Sel]**2) / (Vel[Sel]**2 + VelDisp[Sel]**2))
	B = np.sum(Radius[Sel] * X)
	B_Err = np.sqrt(np.sum((Radius[Sel] * X_Err)**2))
	A = np.sum(Radius[Sel] * Vel[Sel])
	A_Err = np.sqrt(np.sum((Radius[Sel] * VelErr[Sel])**2))
	Total_Lambda_Error = np.sqrt((A_Err/A)**2 + (B_Err/B)**2)
	return Total_Lambda_Error

def AngularMomentumRe(Radius, Vel, VelDisp, Reff, extent):
	Sel = np.where((Radius/Reff < extent+0.05) & (Radius/Reff > extent-0.05))
	Inner_Lambda = np.sum(Radius[Sel]*np.fabs(Vel[Sel]))/np.sum(Radius[Sel]*np.sqrt(Vel[Sel]**2 + VelDisp[Sel]**2))
	return Inner_Lambda

def AngularMomentumReError(Radius, Vel, VelErr, VelDisp, VelDispErr, Reff, extent):
	Sel = np.where((Radius/Reff < extent+0.05) & (Radius/Reff > extent-0.05))
	X = np.sqrt(Vel[Sel]**2 + VelDisp[Sel]**2)
	X_Err = np.sqrt((Vel[Sel]**2 * VelErr[Sel]**2 + VelDisp[Sel]**2 * VelDispErr[Sel]**2) / (Vel[Sel]**2 + VelDisp[Sel]**2))
	B = np.sum(Radius[Sel] * X)
	B_Err = np.sqrt(np.sum((Radius[Sel] * X_Err)**2))
	A = np.sum(Radius[Sel] * Vel[Sel])
	A_Err = np.sqrt(np.sum((Radius[Sel] * VelErr[Sel])**2))
	Inner_Lambda_Error = np.sqrt((A_Err/A)**2 + (B_Err/B)**2)
	return Inner_Lambda_Error

def Flux(GalName, Radius):
	if len(GalName) < 7:
		GalName = 'NGC0'+GalName.split('C')[1]
	try:
		MGE = MGE_Source[GalName]
	except:
		MGE = 'S13'
	surf, sigma, qobs = MGEfinder(MGE, GalName)


	Gaussian = np.zeros(len(Radius))

	for ii in range(len(surf)):
		Gaussian += surf[ii]*np.exp((-Radius**2)/(2*sigma[ii]**2))

	return Gaussian

def MGETotalLuminosity(GalName):
	if len(GalName) < 7:
		GalName = 'NGC0'+GalName.split('C')[1]
	MGE = MGE_Source[GalName]
	surf, sigma, qobs = MGEfinder(MGE, GalName)
	TotalLuminosity = np.sum(2*np.pi*surf*sigma**2*qobs)
	return TotalLuminosity

def Flux_Spitzer(GalName, Radius):
	# given a radial point, I want to calculate the associated luminosity. 
	Radius_Flux, Luminosity, Error = SpitzerProfileFinder(GalName)

	Output = []
	for jj in range(len(Radius)):
		if Radius[jj] < Radius_Flux[0]:
			# print 'lower'
			value = Luminosity[0]
		else:
			match = False
			for ii in range(len(Radius_Flux)):
				if (Radius_Flux[ii] < Radius[jj]) & (match == False):
					lower = ii
				elif (match == False) & (Radius_Flux[ii] > Radius[jj]):
					upper = ii
					match = True	
			
			# elif Radius[jj] > Radius_Flux[]
			value = Luminosity[lower] + ((Radius[jj] - Radius_Flux[lower]) / (Radius_Flux[upper] - Radius_Flux[lower])) * (Luminosity[upper] - Luminosity[lower])
		Output.append(value)
	Output = np.array(Output)
	Output[np.where(np.isfinite(Output) == False)] = 0
	return Output

def SpitzerProfileFinder(GalName):
	DropboxDirectory = os.getcwd().split('Dropbox')[0]
	Spitzer_path = DropboxDirectory+'Dropbox/PhD_Analysis/Data/Spitzer_Luminosity'
	filename = glob.glob(Spitzer_path+'/ngc'+GalName.split('C')[1]+'*logscale.ell')[0]
	# print filename
	file=open(filename, 'r')
	lines=file.readlines()
	file.close()

	Radius, Luminosity, Error = [], [], []
	for line in lines:
		if line[0] != '#':
			try:
				Radius.append(1.22*float(line.split()[1])) # to convert pixel scale to arcsecond scale
				AbsMag = float(line.split()[18]) + 5 - 5*np.log10(distMpc[GalName] * 10**6) # need to convert magnitude into absolute magnitude to convert to luminosity. 
				Luminosity.append(10**((3.24-AbsMag)/2.5))# the solar magnitude in the Spitzer 3.6 micron band is 3.24
				Error.append(10**((3.24-(float(line.split()[18]) + float(line.split()[19])))/2.5) - 10**((3.24-float(line.split()[18]))/2.5))
			except:
				# print 'empty line'
				Line = True

	Radius = np.array(Radius)
	Luminosity = np.array(Luminosity)
	Error = np.array(Error)
	return Radius, Luminosity, Error

def SpitzerMagProfileFinder(GalName, Spitzer_path):
    DropboxDirectory = os.getcwd().split('Dropbox')[0]
    filename = glob.glob(Spitzer_path+'/ngc'+GalName.split('C')[1]+'*logscale.ell')[0]


    Arcseconds, Mag, LowerErr, UpperErr = np.loadtxt(filename, usecols = [1, 18, 19, 20], skiprows = 5, unpack = True)
    Radius = 1.22*Arcseconds # to convert pixel scale to arcsecond scale
    Error = np.maximum(LowerErr, UpperErr) # taking the largest uncertainty out of the
                                           # upper and the lower uncertainty. 
    return Radius, Mag, Error

def Ellipticity_Spitzer(GalName, Radius):
	# given a radial point, I want to calculate the associated luminosity. 
	Radius_Ellipticity, Ellipticity, Error = SpitzerEllipticityFinder(GalName)

	Output = []
	for jj in range(len(Radius)):
		if Radius[jj] < Radius_Ellipticity[0]:
			# print 'lower'
			value = Ellipticity[0]
		elif Radius[jj] > Radius_Ellipticity[-1]:
			# print 'lower'
			value = Ellipticity[-1]
		else:
			match = False
			for ii in range(len(Radius_Ellipticity)):
				if (Radius_Ellipticity[ii] < Radius[jj]) & (match == False):
					lower = ii
				elif (match == False) & (Radius_Ellipticity[ii] > Radius[jj]):
					upper = ii
					match = True	
			
			value = Ellipticity[lower] + ((Radius[jj] - Radius_Ellipticity[lower]) / (Radius_Ellipticity[upper] - Radius_Ellipticity[lower]))\
			 * (Ellipticity[upper] - Ellipticity[lower])
		Output.append(value)
	Output = np.array(Output)
	return Output

def SpitzerEllipticityFinder(GalName):
	DropboxDirectory = os.getcwd().split('Dropbox')[0]
	Spitzer_path = DropboxDirectory+'Dropbox/PhD_Analysis/Data/Spitzer_Luminosity'
	filename = glob.glob(Spitzer_path+'/ngc'+GalName.split('C')[1]+'*logscale.ell')[0]

	Radius, Ellipticity, Error = np.loadtxt(filename, unpack = True, usecols = [1, 6, 7], comments = '#', dtype = 'str')
	Selection = np.where(Error != 'INDEF')
	Radius = Radius[Selection].astype(float)
	Ellipticity = Ellipticity[Selection].astype(float)
	Error = Error[Selection].astype(float)
	Radius = 1.22*Radius

	return Radius, Ellipticity, Error

def AngularMomentumReFlux(GalName, Radius, Vel, VelDisp, GalRA, GalDec, phi, Ellipticity, Reff, extent, FluxSource = 'MGE'):
	Sel = np.where(Radius/Reff < extent)
	if FluxSource == 'MGE':
		Weighting = Flux(GalName, Radius[Sel])
	elif FluxSource == 'Spitzer':
		Weighting = Flux_Spitzer(GalName, Radius[Sel])
	Inner_Lambda = np.sum(Weighting*Radius[Sel]*np.fabs(Vel[Sel]))/np.sum(Weighting*Radius[Sel]*np.sqrt(Vel[Sel]**2 + VelDisp[Sel]**2))
	return Inner_Lambda

def AngularMomentumFlux(GalName, Radius, Vel, VelDisp, GalRA, GalDec, phi, Ellipticity, maximumRadius, FluxSource = 'MGE'):
	Sel = np.where(Radius < maximumRadius)
	# print Sel
	if FluxSource == 'MGE':
		Weighting = Flux(GalName, Radius[Sel])
	elif FluxSource == 'Spitzer':
		Weighting = Flux_Spitzer(GalName, Radius[Sel])
	Total_Lambda = np.sum(Weighting*Radius[Sel]*np.fabs(Vel[Sel]))/np.sum(Weighting*Radius[Sel]*np.sqrt(Vel[Sel]**2 + VelDisp[Sel]**2))
	return Total_Lambda

def AngularMomentumFluxError(GalName, Radius, Vel, VelErr, VelDisp, VelDispErr, maximumRadius, FluxSource = 'MGE'):
	Sel = np.where(Radius < maximumRadius)
	if FluxSource == 'MGE':
		Weighting = Flux(GalName, Radius[Sel])
	elif FluxSource == 'Spitzer':
		Weighting = Flux_Spitzer(GalName, Radius[Sel])
	X = np.sqrt(Vel[Sel]**2 + VelDisp[Sel]**2)
	X_Err = np.sqrt((Vel[Sel]**2 * VelErr[Sel]**2 + VelDisp[Sel]**2 * VelDispErr[Sel]**2) / (Vel[Sel]**2 + VelDisp[Sel]**2))
	B = np.sum(Weighting * Radius[Sel] * X)
	B_Err = np.sqrt(np.sum((Weighting * Radius[Sel] * X_Err)**2))
	A = np.sum(Weighting * Radius[Sel] * Vel[Sel])
	A_Err = np.sqrt(np.sum((Weighting * Radius[Sel] * VelErr[Sel])**2))
	Total_Lambda_Error = np.sqrt((A_Err/A)**2 + (B_Err/B)**2)
	return Total_Lambda_Error

def AngularMomentumFluxAnnulus(GalName, Radius, Vel, VelDisp, GalRA, GalDec, phi, Ellipticity, minimumRadius, maximumRadius):
	Sel = np.where((Radius < maximumRadius) & (Radius > minimumRadius))
	Weighting = Flux(GalName, Radius[Sel])
	Total_Lambda = np.sum(Weighting*Radius[Sel]*np.fabs(Vel[Sel]))/np.sum(Weighting*Radius[Sel]*np.sqrt(Vel[Sel]**2 + VelDisp[Sel]**2))
	return Total_Lambda

def krigingFileRead(Kriging_path, GalName):
	dtype1 = np.dtype([('ra', 'float'), ('dec', 'float'), ('vel', '|S17'), ('velerr', '|S17')])
	Array = np.loadtxt(Kriging_path+str(GalName)+'/'+str(GalName)+'_gridKrig_vel.txt', dtype = dtype1)
	
	Valid = np.where((Array['vel'] != 'NA') & (Array['vel'] != 'nan'))
	RA = Array['ra'][Valid]
	Dec = Array['dec'][Valid]
	Vel = np.asfarray(Array['vel'][Valid], dtype='float')
	VelErr = np.asfarray(Array['velerr'][Valid], dtype='float')
	
	dtype1 = np.dtype([('ra', 'float'), ('dec', 'float'), ('veldisp', '|S17'), ('veldisperr', '|S17')])
	Array = np.loadtxt(Kriging_path+str(GalName)+'/'+str(GalName)+'_gridKrig_sigma.txt', dtype = dtype1)
	
	Valid = np.where((Array['veldisp'] != 'NA') & (Array['veldisp'] != 'nan'))
	VelDisp = np.asfarray(Array['veldisp'][Valid], dtype='float')
	VelDispErr = np.asfarray(Array['veldisperr'][Valid], dtype='float')
	return (RA, Dec, Vel, VelErr, VelDisp, VelDispErr)


def krigingFileReadAll(Kriging_path, GalName):
        dtype1 = np.dtype([('ra', 'float'), ('dec', 'float'), ('vel', '|S17'), ('velerr', '|S17')])
        Array = np.loadtxt(Kriging_path+str(GalName)+'/'+str(GalName)+'_gridKrig_vel.txt', dtype = dtype1)
        
        RA = Array['ra']
        Dec = Array['dec']
        Vel = np.asfarray(Array['vel'], dtype='float')
        VelErr = np.asfarray(Array['velerr'], dtype='float')
        
        dtype1 = np.dtype([('ra', 'float'), ('dec', 'float'), ('veldisp', '|S17'), ('veldisperr', '|S17')])
        Array = np.loadtxt(Kriging_path+str(GalName)+'/'+str(GalName)+'_gridKrig_sigma.txt', dtype = dtype1)
        
        VelDisp = np.asfarray(Array['veldisp'], dtype='float')
        VelDispErr = np.asfarray(Array['veldisperr'], dtype='float')
        return (RA, Dec, Vel, VelErr, VelDisp, VelDispErr)

def krigingFileReadAtlas(Kriging_path, GalName):
	dtype1 = np.dtype([('ra', 'float'), ('dec', 'float'), ('vel', '|S17'), ('velerr', '|S17')])
	Array = np.loadtxt(Kriging_path+str(GalName)+'/'+str(GalName)+'_SAURON_gridKrig_vel.txt', dtype = dtype1)
	
	Valid = np.where((Array['vel'] != 'NA') & (Array['vel'] != 'nan'))
	RA = Array['ra'][Valid]
	Dec = Array['dec'][Valid]
	Vel = np.asfarray(Array['vel'][Valid], dtype='float')
	VelErr = np.asfarray(Array['velerr'][Valid], dtype='float')
	
	dtype1 = np.dtype([('ra', 'float'), ('dec', 'float'), ('veldisp', '|S17'), ('veldisperr', '|S17')])
	Array = np.loadtxt(Kriging_path+str(GalName)+'/'+str(GalName)+'_SAURON_gridKrig_sigma.txt', dtype = dtype1)
	
	Valid = np.where((Array['veldisp'] != 'NA') & (Array['veldisp'] != 'nan'))
	VelDisp = np.asfarray(Array['veldisp'][Valid], dtype='float')
	VelDispErr = np.asfarray(Array['veldisperr'][Valid], dtype='float')
	return (RA, Dec, Vel, VelErr, VelDisp, VelDispErr)



def rotationParameterBootstrap(radius, error, Radius, Vel, VelErr, VelDisp, VelDispErr):
	# radius - the desired radius at which the parameter is measured
	# error - the spread in radius over which SKiMS points are used to measure the lambda_R parameter
	# Radius, Vel, Sigma - radius, velocity and velocity dispersion of points. 

	# # first I need to identify the whole possible sample at each annulus. 
	Annulus_Range = []
	for i in range(len(Radius)):
		if Radius[i] > (radius-error) and Radius[i] < (radius+error):
			Annulus_Range.append(i)
	Annulus_Range = np.array(Annulus_Range)
	Vel_Annulus = Vel[Annulus_Range]
	VelErr_Annulus = VelErr[Annulus_Range]
	VelDisp_Annulus = VelDisp[Annulus_Range]
	VelDispErr_Annulus = VelDispErr[Annulus_Range]
	Radius_Annulus = Radius[Annulus_Range]

	'''
	# implementing the SAURON V/sigma definition:
	'''
	Total_RotationParameter = np.sqrt( (np.sum(Vel_Annulus**2)) / (np.sum(VelDisp_Annulus**2)) )
	Total_RotationParameter_Error = 0.0

	return Total_RotationParameter, Total_RotationParameter_Error




def Rot_Parm_Flux(GalName, Vel, VelDisp, Radius, UpperRadius, FluxSource = 'MGE'):
	Sel = np.where(Radius <= UpperRadius)
	if FluxSource == 'MGE':
		Weighting = Flux(GalName, Radius[Sel])
	elif FluxSource == 'Spitzer':
		Weighting = Flux_Spitzer(GalName, Radius[Sel])
	Rotation_Parameter = np.sqrt(np.sum(Weighting * Vel[Sel]**2) / np.sum(Weighting * VelDisp[Sel]**2))
	return Rotation_Parameter

def Rot_Parm_Flux_Error(GalName, Vel, VelErr, VelDisp, VelDispErr, Radius, UpperRadius, FluxSource = 'MGE'):
	Sel = np.where(Radius <= UpperRadius)
	if FluxSource == 'MGE':
		Weighting = Flux(GalName, Radius[Sel])
	elif FluxSource == 'Spitzer':
		Weighting = Flux_Spitzer(GalName, Radius[Sel])
	Rotation_Parameter = np.sqrt(np.sum(Weighting * Vel[Sel]**2) / np.sum(Weighting * VelDisp[Sel]**2))
	A = np.sum(Weighting * Vel[Sel]**2)
	A_err = np.sum( Weighting * 2 * VelErr[Sel] * np.absolute(Vel[Sel]))
	B = np.sum(Weighting * VelDisp[Sel]**2)
	B_err = np.sum( Weighting * 2 * VelDispErr[Sel] * VelDisp[Sel]) 
	Rotation_Parameter_Err = Rotation_Parameter * 0.5 * np.sqrt((A_err/A) + (B_err/B))
	return Rotation_Parameter_Err


def KMagToMass(Magnitude):
	Mag_Solar = 3.28
	Mass = 10**((Mag_Solar-Magnitude)/2.5)
	return Mass

def arcsecToKpc(arcsecondArray, distanceKpc):
	kpcArray = arcsecondArray * distanceKpc / 206265
	return kpcArray

def MGEfinder(MGE, GalName, units = 'original'):
  DropboxDirectory = os.getcwd().split('Dropbox')[0]
  MGE_path = DropboxDirectory+'Dropbox/PhD_Analysis/Analysis/JAM Modelling/mge_parameters/'
 
  if GalName == 'NGC0821':
      GalName = 'NGC821'
 
  def ArrayExtractorArcsec(filename):
    file=open(filename, 'r')
    lines=file.readlines()
    file.close()
    # need to calculate a conversion factor between L/pc^2 and L/"^2. 
    try:
     Conversion = ((distMpc[GalName]*10**6) / (206265))**2
    except:
      Conversion = ((distMpc_ATLAS[GalName]*10**6) / (206265))**2
  
    surf, sigma, qobs = [], [], []
      
    i=0
    for line in lines:
        if i != 0 and line[0] != '#':
            surf.append(10**float(line.split()[0]) * Conversion)
            sigma.append(10**float(line.split()[1]))
            qobs.append(float(line.split()[2]))
            i += 1
        else:
            i += 1
    return np.array(surf), np.array(sigma), np.array(qobs)
 
  def ArrayExtractorKpc(filename):
    file=open(filename, 'r')
    lines=file.readlines()
    file.close()
    # need to calculate a conversion factor between L/pc^2 and L/kpc^2. 
    Conversion = 10**6
 
    surf, sigma, qobs = [], [], []
     
    i=0
    for line in lines:
        if i != 0 and line[0] != '#':
            surf.append(float(line.split()[0]) * Conversion)# converting this value from L/pc^2 to L/kpc^2. 
            sigma.append(arcsecToKpc(float(line.split()[1]), distMpc[GalName]*1000)) # converting from arcseconds to kpc
            qobs.append(float(line.split()[2]))
            i += 1
        else:
            i += 1
    return np.array(surf), np.array(sigma), np.array(qobs)
 
  def ArrayExtractorOriginal(filename):
    file=open(filename, 'r')
    lines=file.readlines()
    file.close()
    surf, sigma, qobs = [], [], []
     
    i=0
    for line in lines:
        if i != 0 and line[0] != '#':
            surf.append(float(line.split()[0]))
            sigma.append(float(line.split()[1]))
            qobs.append(float(line.split()[2]))
            i += 1
        else:
            i += 1
    return np.array(surf), np.array(sigma), np.array(qobs)
 
  if MGE == 'S13':
    filename = MGE_path+'/mge_parameters_Scott13/mge_'+GalName+'.txt'
 
  elif MGE == 'S09':
    filename = MGE_path+'/mge_parameters_Scott09/mge_'+GalName+'_Scott09.txt'
 
  elif MGE == 'C06':
    filename = MGE_path+'/mge_parameters_Cappellari06/mge_'+GalName+'_C06.txt'
 
  elif MGE == 'E99':
    filename = MGE_path+'/mge_parameters_Emsellem99/mge_'+GalName+'.txt'
 
  elif MGE == 'Sabine':
    filename = MGE_path+'/mge_parameters_Sabine/mge_'+GalName+'.txt'
 
  if units == 'arcsec':
    print 'MGE used with units of arcsec'
    return ArrayExtractorArcsec(filename)
  elif units == 'kpc':
    print 'MGE used with units of kpc'
    return ArrayExtractorKpc(filename)
  elif units == 'original':
    # print 'MGE used with original units'
    return ArrayExtractorOriginal(filename)


def MovingAverage(x, y, PointsPerBin = 15):
    SortedXArray = np.argsort(x)
    Y, Y_upper, Y_lower = [], [], []
    Bins = int(len(x) / PointsPerBin) +1
    X = []
    for ii in range(Bins):
        Sel = SortedXArray[ii*PointsPerBin+1:ii*PointsPerBin+PointsPerBin]
        Y.append(np.median(y[Sel]))
        Y_upper.append(np.percentile(y[Sel], 84))
        Y_lower.append(np.percentile(y[Sel], 16))
        X.append(np.mean(x[Sel]))
    
    Y = np.array(Y)
    Y_upper = np.array(Y_upper)
    Y_lower = np.array(Y_lower)
    X = np.array(X)
    return X, Y, Y_upper, Y_lower

def localEllipticity(X, Y, Density, RadialInterval = 2.0, annulusWidth = 0.3, minimumRadius = 1.0):
	phi=(90*np.pi/180-np.pi/2.0)
	maximumRadius = np.max(X)
	Radius_EllipProfile = np.arange(minimumRadius, maximumRadius, RadialInterval)
	Ellipticity_EllipProfile, Uncertainty_EllipProfile = [], []
	for majorAxisRadius in Radius_EllipProfile: # sampling the local ellipticity along the major axis in kpc.
		# first calculating the mean density on the major axis at this radius
		MeanDensity = np.mean(Density[np.where( (np.abs(X) > majorAxisRadius - annulusWidth) & \
			(np.abs(X) < majorAxisRadius + annulusWidth) & (np.abs(Y) < annulusWidth))]) 
		# now at each radial increment, test various ellipticity values to check which has the 
		# lowest scatter in the density distribution along that ellipse. 
		try:

			# coarse text first
			ellipticityRange = np.arange(0.05, 0.95, 0.05)
			densityDistribution = []
			meanDensityDeviation = []
			for coarseEllipticity in ellipticityRange:
				# recalculate the radius of each pixel according to the test ellipticity value
				Radius = radiusArray(X, Y, phi, coarseEllipticity)
	
				# now selecting the pixels that correspond to the major axis value
				circularisedRadius = np.sqrt(1 - coarseEllipticity) * majorAxisRadius
				annulusSelection = np.where( (Radius > (circularisedRadius - annulusWidth)) & (Radius < (circularisedRadius + annulusWidth)) )
				densityDistribution.append(np.std(Density[annulusSelection]))
				meanDensityDeviation.append(abs(MeanDensity - np.mean(Density[annulusSelection])))

	
			# identifying the ellipticity that is the closest match from the coarse test:
			# firstGuessEllipticity = np.mean(ellipticityRange[np.where(densityDistribution == np.min(densityDistribution))])
			firstGuessEllipticity = np.mean(ellipticityRange[np.where(meanDensityDeviation == np.min(meanDensityDeviation))])
			
			# fine test second
			ellipticityRange = np.arange(firstGuessEllipticity-0.05, firstGuessEllipticity+0.05, 0.01)
			densityDistribution = []
			meanDensityDeviation = []

			for fineEllipticity in ellipticityRange:
				if fineEllipticity != 1.0:
					# recalculate the radius of each pixel according to the test ellipticity value
					Radius = radiusArray(X, Y, phi, fineEllipticity)
		
					# now selecting the pixels that correspond to the major axis value
					circularisedRadius = np.sqrt(1 - fineEllipticity) * majorAxisRadius
					annulusSelection = np.where( (Radius > (circularisedRadius - annulusWidth)) & (Radius < (circularisedRadius + annulusWidth)) )
					densityDistribution.append(np.std(Density[annulusSelection]))
					meanDensityDeviation.append(abs(MeanDensity - np.mean(Density[annulusSelection])))

			# FinalArray = ellipticityRange[np.where(densityDistribution == np.min(densityDistribution))]
			FinalArray = ellipticityRange[np.where(meanDensityDeviation == np.min(meanDensityDeviation))]
			if np.isfinite(np.mean(FinalArray)):
				Ellipticity_EllipProfile.append(np.mean(FinalArray))
				# now calculating the uncertainty in the ellipticity. Here, the major source of uncertainty is simply the resolution. 
				Err = (np.max(FinalArray) - np.min(FinalArray)) / 2
				Uncertainty_EllipProfile.append(Err)
			else:
				Ellipticity_EllipProfile.append(0.5)
				Uncertainty_EllipProfile.append(0.5)
		except:
			Ellipticity_EllipProfile.append(0.5)
			Uncertainty_EllipProfile.append(0.5)
	Ellipticity_EllipProfile = np.array(Ellipticity_EllipProfile)
	Uncertainty_EllipProfile = np.array(Uncertainty_EllipProfile)

	return 	Radius_EllipProfile, Ellipticity_EllipProfile, Uncertainty_EllipProfile

def localEllipticityFitting(X, Y, Density, Radius_EllipProfile, annulusWidth = 0.3):
	phi=(90*np.pi/180-np.pi/2.0)
	# maximumRadius = np.max(X)
	# Radius_EllipProfile = np.arange(minimumRadius, maximumRadius, RadialInterval)
	Ellipticity_EllipProfile, Uncertainty_EllipProfile = [], []
	for majorAxisRadius in Radius_EllipProfile: # sampling the local ellipticity along the major axis in kpc.
		# first calculating the mean density on the major axis at this radius
		MeanDensity = np.mean(Density[np.where( (np.abs(X) > majorAxisRadius - annulusWidth) & \
			(np.abs(X) < majorAxisRadius + annulusWidth) & (np.abs(Y) < annulusWidth))]) 
		# now at each radial increment, test various ellipticity values to check which has the 
		# lowest scatter in the density distribution along that ellipse. 
		try:

			# coarse text first
			ellipticityRange = np.arange(0.2, 0.8, 0.2)
			densityDistribution = []
			meanDensityDeviation = []
			for coarseEllipticity in ellipticityRange:
				# recalculate the radius of each pixel according to the test ellipticity value
				Radius = radiusArray(X, Y, phi, coarseEllipticity)
	
				# now selecting the pixels that correspond to the major axis value
				circularisedRadius = np.sqrt(1 - coarseEllipticity) * majorAxisRadius
				annulusSelection = np.where( (Radius > (circularisedRadius - annulusWidth)) & (Radius < (circularisedRadius + annulusWidth)) )
				densityDistribution.append(np.std(Density[annulusSelection]))
				meanDensityDeviation.append(abs(MeanDensity - np.mean(Density[annulusSelection])))

	
			# identifying the ellipticity that is the closest match from the coarse test:
			# firstGuessEllipticity = np.mean(ellipticityRange[np.where(densityDistribution == np.min(densityDistribution))])
			firstGuessEllipticity = np.mean(ellipticityRange[np.where(meanDensityDeviation == np.min(meanDensityDeviation))])
			
			# fine test second
			ellipticityRange = np.arange(firstGuessEllipticity-0.2, firstGuessEllipticity+0.2, 0.05)
			densityDistribution = []
			meanDensityDeviation = []

			for fineEllipticity in ellipticityRange:
				if fineEllipticity != 1.0:
					# recalculate the radius of each pixel according to the test ellipticity value
					Radius = radiusArray(X, Y, phi, fineEllipticity)
		
					# now selecting the pixels that correspond to the major axis value
					circularisedRadius = np.sqrt(1 - fineEllipticity) * majorAxisRadius
					annulusSelection = np.where( (Radius > (circularisedRadius - annulusWidth)) & (Radius < (circularisedRadius + annulusWidth)) )
					densityDistribution.append(np.std(Density[annulusSelection]))
					meanDensityDeviation.append(abs(MeanDensity - np.mean(Density[annulusSelection])))

			# FinalArray = ellipticityRange[np.where(densityDistribution == np.min(densityDistribution))]
			FinalArray = ellipticityRange[np.where(meanDensityDeviation == np.min(meanDensityDeviation))]
			Ellipticity_EllipProfile.append(np.mean(FinalArray))
			# now calculating the uncertainty in the ellipticity. Here, the major source of uncertainty is simply the resolution. 
			Err = (np.max(FinalArray) - np.min(FinalArray)) / 2
			Uncertainty_EllipProfile.append(Err)
		except:
			Ellipticity_EllipProfile.append(0.5)
			Uncertainty_EllipProfile.append(0.5)
	Ellipticity_EllipProfile = np.array(Ellipticity_EllipProfile)
	Uncertainty_EllipProfile = np.array(Uncertainty_EllipProfile)

	return 	Radius_EllipProfile, Ellipticity_EllipProfile, Uncertainty_EllipProfile

def DensitytoLuminosity(X, Y, Density):
		PixelScale = abs(X[0] - X[1])
		SurfaceBrightness_lumperpc = (1/((1000*PixelScale)**2)) * Density
		
		return SurfaceBrightness_lumperpc
