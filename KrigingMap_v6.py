# Create velocity map excluding the elements with velocities << overlapping ATLAS3d values

######################################################################################################
# v.3 applies the spatial offset (-6, +9 arcsec) to the points 
# v.4 - Edited by Sabine - SAURON data can now be included in certain percentages (to replicate the 
#                         distribution of the SKiMS points) 
#                        - map resolution can be easily determined at the start.
# v.5 - produced monte carlo iterations of the maps to calculate errors in future.  
# v.6 - gives the option of bi-symmetrising the input maps and storing them in a separate location. 
######################################################################################################

#Retrieve dictionaries
import os, sys, random
DropboxDirectory = os.getcwd().split('Dropbox')[0]
lib_path = os.path.abspath(DropboxDirectory+'Dropbox/PhD_Analysis/Library') 
sys.path.append(lib_path)
from galaxyParametersDictionary_v9 import *
from Sabine_Define import *
from Nicola import *
from KrigingMapping_def_v6 import *
from createMaps__def__ import *
import glob


############ All input parameters to be determined here ################################################
Galaxies = ['NGC3115']
# Gal_Name = ['NGC5846']
# Gal_Name = ['NGC2549', 'NGC4459', 'NGC4474', 'NGC7457'] # name of galaxies to be mapped:
# Gal_Name = Galaxies
To_Plot = ['vel', 'sigma'] # property of galaxy to be mapped:
                   # options: vel; sigma; h3; h4; SN; Z (metallicity only for SAURON maps)
SLUGGSplusATLAS = False                  
SAURON = False
Norris = False
GC = False
addPoints = True # add points to the kriging map?
addContours = True # add contours to the Kriging map?
pixelNumber = 160 # number of pixels in the final map - higher number means larger map resolution. 
symmetrize = False # symmetrising the points to make a symmetric map. 
######################################################################################################

if SAURON:
  # if SAURON/ATLAS3D data is used, download the galaxy parameters from online. 
  from astroquery.vizier import Vizier

  Vizier.ROW_LIMIT = -1
  cat_list = Vizier.get_catalogs('J/MNRAS/414/888/tableb1')
  
  cat = cat_list[0]
  
  cat.colnames
  Galaxies = np.array(cat['Gal'])#[np.where((np.array(cat['Rmax']) <= 2.0) & (np.array(cat['Rmax']) > 1.7))]
  
  b_a = {}
  for GalName in Galaxies:
    b_a[GalName] = 1 - np.array(cat['epse'])[np.where(np.array(cat['Gal']) == GalName)]
  
  cat_list = Vizier.get_catalogs('J/MNRAS/413/813/atlas3d')
  cat = cat_list[0]
  cat.colnames
  
  Reff = {}
  vel0 = {}
  for GalName in Galaxies:
    Reff[GalName] = 10**np.array(cat['log_Re_'])[np.where(np.array(cat['Gal']) == GalName)]
    vel0[GalName] = np.array(cat['Vhel'])[np.where(np.array(cat['Gal']) == GalName)]
  
  exclude = {}
  Theta_Kriging_ATLAS = {}
  Theta_Kriging = {}
  for GalName in Galaxies:
    exclude[GalName] = 'None'
    Theta_Kriging_ATLAS[GalName] = 5
    Theta_Kriging[GalName] = 15
  
  PA0 = {}
  
  x = np.loadtxt('Data/ATLAS3D/Krajnovic11_TableD1.txt', comments = '#', unpack = True, usecols = [0], dtype = 'str')
  y = np.loadtxt('Data/ATLAS3D/Krajnovic11_TableD1.txt', comments = '#', unpack = True, usecols = [1])
  for GalName in Galaxies:
    PA0[GalName] = y[np.where(x == GalName)]

for GalName in Galaxies:
  # try:
  print GalName
  for Property in To_Plot:
    print Property
    if not os.path.exists(str(GalName)): 
      os.mkdir(str(GalName))
    
    # Extract positions and kinematics from masks
    print vel0[GalName]
    vel_sys = float(vel0[GalName])
    if (SAURON == False):
      pathMajor = 'Data/SKiMS_latest/deimos_ppxf_'+str(GalName)+'.fits'
      
      Identifier = {'vel':'VEL', 'sigma':'VELDISP', 'h3':'H3', 'h4':'H4', 'SN':'SN'}
      if (Property=='vrms'):
        krig_positions_vel = retrieve_krigPos(pathMajor, GalName, exclude[GalName], keyword='VEL')
        krig_positions_sig = retrieve_krigPos(pathMajor, GalName, exclude[GalName], keyword='VELDISP')
      else:
        krig_positions = retrieve_krigPos(pathMajor, GalName, exclude[GalName], keyword=Identifier[Property])
      
      if Property == 'SN':
        x_SS, y_SS, z_SS, ez_SS = krig_positions
        ez_SS=[]
        for i in range(len(z_SS)):
          ez_SS.append(1.0)
      elif  Property == 'vrms':
        x_SS, y_SS, z_SS_vel, ez_SS_vel = krig_positions_vel
        x_SS, y_SS, z_SS_sig, ez_SS_sig = krig_positions_sig
        for ii in range(len(z_SS_vel)): # subtracting the systemic velocity. 
          z_SS_vel[ii] = z_SS_vel[ii] - vel_sys
        z_SS = np.sqrt((z_SS_vel)**2 + z_SS_sig**2)
        ez_SS = np.sqrt((ez_SS_vel)**2 + ez_SS_sig**2)
      else:
        x_SS, y_SS, z_SS, ez_SS = krig_positions
    
      if Property=='vel':
        for ii in range(len(z_SS)):
          z_SS[ii] = z_SS[ii] - vel_sys
    
      if symmetrize:
        # here we wish to symmetrise the points along both axis - taking into account whether the 
        # parameter being mapped is antisymmetric or not. 
        if Property == 'SN':
          print 'Will not symmetrise SN maps or MC generated error maps'
        elif Property == 'vel' or Property == 'h3':
          # these maps are antisymmetric
          theta = np.radians(PA0[GalName] - 90.)
          xNew = x_SS*np.cos(theta) - y_SS*np.sin(theta)
          yNew = x_SS*np.sin(theta) + y_SS*np.cos(theta)
    
          print z_SS
          for i in range(len(z_SS)): # flipping along the x axis - needs to be opposite sign
            xNew=np.append(xNew, -xNew[i])
            yNew=np.append(yNew, yNew[i])
            z_SS=np.append(z_SS, -z_SS[i])
            ez_SS=np.append(ez_SS, ez_SS[i])
          for i in range(len(z_SS)): # flipping along y axis - sign remains same
            xNew=np.append(xNew, xNew[i])
            yNew=np.append(yNew, -yNew[i])
            z_SS=np.append(z_SS, z_SS[i])
            ez_SS=np.append(ez_SS, ez_SS[i])
          for i in range(len(z_SS)): # flipping in both axis - sign is swapped
            xNew=np.append(xNew, -xNew[i])
            yNew=np.append(yNew, -yNew[i])
            z_SS=np.append(z_SS, -z_SS[i])
            ez_SS=np.append(ez_SS, ez_SS[i])
    
          theta = np.radians(360.+90.-PA0[GalName])
          x_SS = xNew*np.cos(theta) - yNew*np.sin(theta)
          y_SS = xNew*np.sin(theta) + yNew*np.cos(theta)
        else:
          # these maps are symmetric
          theta = np.radians(PA0[GalName] - 90.)
          xNew = x_SS*np.cos(theta) - y_SS*np.sin(theta)
          yNew = x_SS*np.sin(theta) + y_SS*np.cos(theta)
    
          for i in range(len(z_SS)): # flipping along the x axis - sign remains same
            xNew=np.append(xNew, -xNew[i])
            yNew=np.append(yNew, yNew[i])
            z_SS=np.append(z_SS, z_SS[i])
            ez_SS=np.append(ez_SS, ez_SS[i])
          for i in range(len(z_SS)): # flipping along y axis - sign remains same
            xNew=np.append(xNew, xNew[i])
            yNew=np.append(yNew, -yNew[i])
            z_SS=np.append(z_SS, z_SS[i])
            ez_SS=np.append(ez_SS, ez_SS[i])
          for i in range(len(z_SS)): # flipping in both axis - sign remains same
            xNew=np.append(xNew, -xNew[i])
            yNew=np.append(yNew, -yNew[i])
            z_SS=np.append(z_SS, z_SS[i])
            ez_SS=np.append(ez_SS, ez_SS[i])
    
          theta = np.radians(360.+90.-PA0[GalName])
          x_SS = xNew*np.cos(theta) - yNew*np.sin(theta)
          y_SS = xNew*np.sin(theta) + yNew*np.cos(theta)
    
    ##  extract velocities from ATLAS3d/SAURON
    if SAURON or SLUGGSplusATLAS:
      ATLAS_path = 'Data/ATLAS3D/'
      try:
        # print glob.glob(ATLAS_path+'atlas3d_stellar_kinematics/*/PXF_*'+str(GalName)+'*.fits')
        data=pyfits.getdata(glob.glob(ATLAS_path+'atlas3d_stellar_kinematics/*/PXF_*'+str(GalName)+'*.fits')[0], 0)
        print 'ATLAS3d source file detected: '
        print glob.glob(ATLAS_path+'atlas3d_stellar_kinematics/*/PXF_*'+str(GalName)+'*.fits')[0]
      except:
        data=pyfits.getdata(glob.glob(ATLAS_path+'ATLAS3d/PXF_*'+str(GalName)+'*.fits')[0], 0)
        print 'ATLAS3d source file detected: '
        print glob.glob(ATLAS_path+'ATLAS3d/PXF_*'+str(GalName)+'*.fits')[0]
      
      X, Y=data['xs'], data['ys']
      # SAURON data has flipped the RA axis - we swap this to make it consistent with SKiMS data
      for i in range(len(X)):
        X[i]=-X[i]
      Vel, VelErr=data['vpxf']-vel_sys, data['evpxf']
      if SLUGGSplusATLAS:
        Sig, SigErr=data['spxf']-ATLAS_SigmaOffset[GalName], data['espxf']
      else:
        Sig, SigErr=data['spxf'], data['espxf']
      h3, h3Err=data['h3pxf'], data['eh3pxf']
      h4, h4Err=data['h4pxf'], data['eh4pxf']
      
      x_atlas, y_atlas = numpy.array(X), numpy.array(Y)
      if Property=='vel':
        z_atlas, ez_atlas = numpy.array(Vel), numpy.array(VelErr)
      if Property=='sigma':
        z_atlas, ez_atlas = numpy.array(Sig), numpy.array(SigErr)
      if Property=='h3':
        z_atlas, ez_atlas = numpy.array(h3), numpy.array(h3Err)
      if Property=='h4':
        z_atlas, ez_atlas = numpy.array(h4), numpy.array(h4Err)
      
      if GalName == 'NGC4459': 
        # determining which points correspond to the foreground star, so that they can be masked
        index=[]
        X, Y, Vel=numpy.array(X), numpy.array(Y), numpy.array(Vel)
        for ii in range(len(X)):
         if X[ii] < -17 and X[ii] > -20 and Y[ii] > 5 and Y[ii] < 7.5 and Vel[ii] < 1290:
          index.append(ii)
        # now deleting the points corresponding to the star
        x_atlas, y_atlas = numpy.delete(x_atlas, index, None), numpy.delete(y_atlas, index, None)
        z_atlas, ez_atlas = numpy.delete(z_atlas, index, None), numpy.delete(ez_atlas, index, None)

    if Norris:
      data=pyfits.getdata('Data/Norris06/Norris06_NGC3115.fits', 0)
      
      X, Y=data['ra'], data['dec']
      # SAURON data has flipped the RA axis - we swap this to make it consistent with SKiMS data
      # for i in range(len(X)):
      #   X[i]=-X[i]
      Vel, VelErr=data['vel'], data['errvel']
      Sig, SigErr=data['veldisp'], data['errveldisp']
      h3, h3Err=data['h3'], data['errh3']
      h4, h4Err=data['h4'], data['errh4']
      
      x_norris, y_norris = numpy.array(X), numpy.array(Y)
      if Property=='vel':
        z_norris, ez_norris = numpy.array(Vel), numpy.array(VelErr)
      if Property=='sigma':
        z_norris, ez_norris = numpy.array(Sig), numpy.array(SigErr)
      if Property=='h3':
        z_norris, ez_norris = numpy.array(h3), numpy.array(h3Err)
      if Property=='h4':
        z_norris, ez_norris = numpy.array(h4), numpy.array(h4Err)

    if GC:
      GC_filename = 'Data/SLUGGS_GC_RVs/'+GalName+'_GC_BinnedKinematics.txt'
      # if not os.path.isfile(GC_filename):
      #   GC = False
      #   print 'No GC data'
      # else:
      if Property == 'vel':
        X_gc, Y_gc, Z_gc, eZ_gc = np.loadtxt(GC_filename, unpack = True, usecols = [0, 1, 2, 3])
      elif Property == 'sigma':
        X_gc, Y_gc, Z_gc, eZ_gc = np.loadtxt(GC_filename, unpack = True, usecols = [0, 1, 4, 5])

    # Combining all SAURON and SKiMS data points into arrays 
    
    X, Y, Z, eZ, Vel, eVel, Sig, eSig = [], [], [], [], [], [], [], []
    Check = []
    if SAURON:
      for ii in numpy.arange(len(x_atlas)):
        X.append(x_atlas[ii])
        Y.append(y_atlas[ii])
        Z.append(z_atlas[ii])
        eZ.append(ez_atlas[ii])

    if Norris:
      for ii in numpy.arange(len(x_norris)):
        X.append(x_norris[ii])
        Y.append(y_norris[ii])
        Z.append(z_norris[ii])
        eZ.append(ez_norris[ii])

    if GC:
      for ii in numpy.arange(len(X_gc)):
        X.append(X_gc[ii])
        Y.append(Y_gc[ii])
        Z.append(Z_gc[ii])
        eZ.append(eZ_gc[ii])

    if SLUGGSplusATLAS:
      # I want to only add a selection of ATLAS points
      Sample = xrange(0,len(x_atlas))
      ATLASSelection = numpy.random.choice(Sample, 140, replace = False)
      for ii in ATLASSelection:
        X.append(x_atlas[ii])
        Y.append(y_atlas[ii])
        Z.append(z_atlas[ii])
        eZ.append(ez_atlas[ii])
      for ii in numpy.arange(len(x_SS)):
        X.append(x_SS[ii])
        Y.append(y_SS[ii])
        Z.append(z_SS[ii])
        eZ.append(ez_SS[ii])

    
    if SAURON == False:
      for ii in numpy.arange(len(x_SS)):
         X.append(x_SS[ii])
         Y.append(y_SS[ii])
         Z.append(z_SS[ii])
         eZ.append(ez_SS[ii])
         Check.append(1)
    
    if symmetrize:
      if not os.path.exists(str(GalName)+'/Symmetrized/'):
        os.mkdir(str(GalName)+'/Symmetrized/')
      Directory = str(GalName)+'/Symmetrized'
    else:
      Directory = str(GalName)

    if SAURON:
      SAURONimplementation = 'SAURON'
    elif Norris:
      SAURONimplementation = 'Norris06'
    elif SLUGGSplusATLAS:
      SAURONimplementation = 'SLUGGSandATLAS'
    elif GC:
      SAURONimplementation = 'GC'
    else:
      SAURONimplementation = ''
    
 
    genTable = transpose(numpy.array([numpy.array(X), numpy.array(Y), numpy.array(Z), numpy.array(eZ)]))
    print len(genTable)
    #Saving new files
    fileout = open(Directory+'/'+str(GalName)+'_listElements'+str(SAURONimplementation)+'_'+str(Property)+'.txt', 'wb')
    
    numpy.savetxt(fileout, genTable, delimiter='\t', header='x\ty\tz\terrz')
    fileout.close()
    if SAURON:
      theta_sigma = Theta_Kriging_ATLAS[GalName]
    else:
      theta_sigma = Theta_Kriging[GalName]
    
    # the filename of the gridKrig file generated by R is given by pathOutput+'gridKrig_'+label+'.txt'
    # Therefore to save these as more specific, the pathOutput parameter needs to change. 
    # the inputpath for the KrigingMapPython function then also needs to be modified. 
    if SAURON:
      prefix = Directory+'/'+str(GalName)+'_SAURON_'
    elif Norris:
      prefix = Directory+'/'+str(GalName)+'_Norris_'
    elif SLUGGSplusATLAS:
      prefix = Directory+'/'+str(GalName)+'_SLUGGSandATLAS_'
    elif GC:
      prefix = Directory+'/'+str(GalName)+'_GC_'
    else:
      prefix = Directory+'/'+str(GalName)+'_'

    if Property=='vel':
      dummy = KrigingR(Directory+'/'+str(GalName)+'_listElements'+str(SAURONimplementation)+'_vel.txt', visualize=False, 
              theta_r = theta_sigma, coeff_r = 3, savePdf = False, 
              pathOutput = prefix, label='Vel', sizePixelMap=pixelNumber) 
      dummy = KrigingMapPython(prefix, GalName, genTable, label='Vel',
                                 limits = [-300, 300], sizePixelMap=pixelNumber, 
                                 points=addPoints, contours=addContours) #For the visualization
    elif Property=='sigma':
      dummy = KrigingR(Directory+'/'+str(GalName)+'_listElements'+str(SAURONimplementation)+'_sigma.txt', visualize=False, 
              theta_r = theta_sigma, coeff_r = 3, savePdf = False, 
              pathOutput = prefix, label='sigma', sizePixelMap=pixelNumber) 
      dummy = KrigingMapPython(prefix, GalName, genTable, label='sigma',
                                 limits = [0, 400], sizePixelMap=pixelNumber, 
                                 points=addPoints, contours=addContours) #For the visualization
    elif Property=='h3':
      dummy = KrigingR(Directory+'/'+str(GalName)+'_listElements'+str(SAURONimplementation)+'_h3.txt', visualize=False, 
              theta_r = theta_sigma, coeff_r = 3, savePdf = False, 
              pathOutput = prefix, label='h3', sizePixelMap=pixelNumber) 
      dummy = KrigingMapPython(prefix, GalName, genTable, label='h3',
                                 limits = [-0.4, 0.4], sizePixelMap=pixelNumber, 
                                 points=addPoints, contours=addContours) #For the visualization
    elif Property=='h4':
      dummy = KrigingR(Directory+'/'+str(GalName)+'_listElements'+str(SAURONimplementation)+'_h4.txt', visualize=False, 
              theta_r = theta_sigma, coeff_r = 3, savePdf = False, 
              pathOutput = prefix, label='h4', sizePixelMap=pixelNumber) 
      dummy = KrigingMapPython(prefix, GalName, genTable, label='h4',
                                 limits = [-0.4, 0.4], sizePixelMap=pixelNumber, 
                                 points=addPoints, contours=addContours) #For the visualization
    elif Property=='SN':
      dummy = KrigingR(Directory+'/'+str(GalName)+'_listElements'+str(SAURONimplementation)+'_SN.txt', visualize=False, 
              theta_r = theta_sigma, coeff_r = 3, savePdf = False, 
              pathOutput = prefix, label='SN', sizePixelMap=pixelNumber) 
      dummy = KrigingMapPython(prefix, GalName, genTable, label='SN',
                                 limits = [0, 200], sizePixelMap=pixelNumber, 
                                 points=addPoints, contours=addContours) #For the visualization
    elif Property=='Z':
      dummy = KrigingR(Directory+'/'+str(GalName)+'_listElements'+str(SAURONimplementation)+'_Z.txt', visualize=False, 
              theta_r = theta_sigma, coeff_r = 3, savePdf = False, 
              pathOutput = prefix, label='Z', sizePixelMap=pixelNumber) 
      dummy = KrigingMapPython(prefix, GalName, genTable, label='Z',
                                 limits = [-0.4, 0.4], sizePixelMap=pixelNumber, 
                                 points=addPoints, contours=addContours) #For the visualization
    elif Property=='vrms':
      dummy = KrigingR(Directory+'/'+str(GalName)+'_listElements'+str(SAURONimplementation)+'_vrms.txt', visualize=False, 
              theta_r = theta_sigma, coeff_r = 3, savePdf = False, 
              pathOutput = prefix, label='vrms', sizePixelMap=pixelNumber) 
      dummy = KrigingMapPython(prefix, GalName, genTable, label='vrms',
                                 limits = [0, 300], sizePixelMap=pixelNumber, 
                                 points=addPoints, contours=addContours) #For the visualization
    
    
  
    print 'DONE'   
  # except:
  #   print 'error with', GalName 