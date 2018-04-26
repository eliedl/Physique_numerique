#!/usr/bin/env python
# -*- coding: utf-8 -*-
# TP reconstruction TDM (CT)
# Prof: Philippe Després
# programme: Dmitri Matenine (dmitri.matenine.1@ulaval.ca)


# libs
import numpy as np
import time
import matplotlib.pyplot as plt

# local files
import geometry as geo
import util as util
import CTfiltre as CTfilter
#fuckoff
def distance_point_droite(x, y, x0, y0, theta):
    dx =  (x0-x)*np.cos(theta) - (y0-y)*np.sin(theta)
    return dx

## créer l'ensemble de données d'entrée à partir des fichiers
def readInput():
    # lire les angles
    [nbprj, angles] = util.readAngles(geo.dataDir+geo.anglesFile)

    print("nbprj:",nbprj)
    print("angles min and max (rad):")
    print("["+str(np.min(angles))+", "+str(np.max(angles))+"]")

    # lire le sinogramme
    [nbprj2, nbpix2, sinogram] = util.readSinogram(geo.dataDir+geo.sinogramFile)

    if nbprj != nbprj2:
        print("angles file and sinogram file conflict, aborting!")
        exit(0)

    if geo.nbpix != nbpix2:
        print("geo description and sinogram file conflict, aborting!")
        exit(0)

    return [nbprj, angles, sinogram]


## reconstruire une image TDM en mode retroprojection
def laminogram():
    
    [nbprj, angles, sinogram] = readInput()

    # initialiser une image reconstruite
    image = np.zeros((geo.nbvox, geo.nbvox))
    # "etaler" les projections sur l'image
    # ceci sera fait de façon "voxel-driven"
    # pour chaque voxel, trouver la contribution du signal reçu
    x0 = y0 = geo.nbvox/2
    half_sino = sinogram.shape[1]/2
    for j in range(geo.nbvox): # colonnes de l'image
        #print("working on image column: "+str(j+1)+"/"+str(geo.nbvox))
        for i in range(geo.nbvox): # lignes de l'image
                #rotation et translation de l'espace
                dist = np.array(distance_point_droite(i, j, x0, y0, angles))
                dist = np.sqrt(1.9) * dist * half_sino / (geo.nbvox)
                index_inf = (np.floor(half_sino + dist)).astype(int)
                index_sup = (np.ceil(half_sino + dist)).astype(int)
                residue = dist%1
                x = np.arange(0, len(sinogram))
                inf = sinogram[(x, index_inf)]
                sup = sinogram[(x, index_sup)]
                image[j,i] += np.sum(residue*(sup-inf) + inf)

    #util.saveImage(image, "lam")


## reconstruire une image TDM en mode retroprojection
def backproject():
    
    [nbprj, angles, sinogram] = readInput()
    
    # initialiser une image reconstruite
    image = np.zeros((geo.nbvox, geo.nbvox))
    
    ### option filtrer ###
    #CTfilter.filterSinogram(sinogram)
    ######
    
    # "etaler" les projections sur l'image
    # ceci sera fait de façon "voxel-driven"
    # pour chaque voxel, trouver la contribution du signal reçu
    x0 = y0 = geo.nbvox / 2
    half_sino = sinogram.shape[1] / 2
    for j in range(geo.nbvox): # colonnes de l'image
        #print("working on image column: "+str(j+1)+"/"+str(geo.nbvox))
        for i in range(geo.nbvox): # lignes de l'image
            dist = np.array(distance_point_droite(i, j, x0, y0, angles))
            dist = np.sqrt(1.9) * dist * half_sino / (geo.nbvox)
            index_inf = (half_sino + dist).astype(int)
            x = np.arange(0,len(sinogram))
            image[j, i] += np.sum(sinogram[(x,index_inf)])

    #util.saveImage(image, "fbp")


def rotate_line(x, x0, y0, theta):
    nx = -(x - x0)*np.sin(theta) + x0
    ny = (x - x0)*np.cos(theta) + y0
    return nx, ny

## reconstruire une image TDM en mode retroprojection
def reconFourierSlice():
    sinogram = np.loadtxt('data/sino_fftfiltre.txt')

    angles = np.loadtxt('data/angles.txt')
    [nbprj, angles, sinogram] = readInput()

    # m = int(np.round(0.5*(np.sqrt(2)-1)*sinogram.shape[1]))
    # pad = np.zeros((720, m))
    # sinogram = np.hstack((pad, sinogram, pad))

    # Initialisation de l'image et de la TdF du sinogramme
    image = np.zeros((geo.nbvox, geo.nbvox), 'complex')
    fft_sino = np.zeros(sinogram.shape, 'complex')
    print(sinogram.shape)
    for i in range(len(sinogram)):
        # TdF de la slice
        line = sinogram[i, :]
        fftline = np.fft.fft(line)
        fftshift = np.fft.fftshift(fftline)

        # Filtration des basses fréquences
        freq = np.fft.fftfreq(len(line), 1)
        freqshift = np.fft.fftshift(freq)
        fftline = np.conj(np.abs(freqshift) * fftshift)

        fft_sino[i, :] = fftline
    sinogram = fft_sino
    x0 = y0 = len(sinogram[0]) / 2
    half_sino = sinogram.shape[1] / 2
    indexes = np.arange(len(sinogram[0]))
    factor = geo.nbvox / (len(sinogram[0])*np.sqrt(1.9))
    count = 0
    dilat = geo.nbvox/np.sqrt(1.9)
    linefunc = lambda x, a, b: a*x+b
    for j in range(geo.nbvox):  # colonnes de l'image
        print("working on image column: " + str(j + 1) + "/" + str(geo.nbvox))
        for i in range(geo.nbvox):
            x, y  = (j - geo.nbvox/2),(i - geo.nbvox/2)
            if x == 0  or y == 0:
                continue
            theta = np.arctan(y/x)%(2*np.pi)
            r = np.linalg.norm([x,y])#/factor
            index = np.argmin(np.abs(angles - theta))
            rindex = linefunc(r, -168/68, 168)
            rinf = (np.floor(rindex)).astype(int)
            rsup = (np.ceil(rindex)).astype(int)
            reste = rindex%1
            line = sinogram[index]
            print(index)
            inf = line[rinf]
            sup = line[rsup]
            image[j, i] += reste * (sup - inf) + inf
        #plt.show()
    plt.imshow(image.imag, cmap = "Greys")
    plt.show()



    '''
    for line, angle in zip(sinogram, angles):

        xvalues, yvalues = rotate_line(indexes, x0, y0, angle)
        xindex, yindex = (factor *xvalues).astype(int), (factor*yvalues).astype(int)
        diff = np.sqrt(((xindex - xvalues*factor)**2 + (xindex - xvalues*factor)**2))
        for i,x,y, val in zip(diff, xindex, yindex, line):
            if diffs[x,y] > i:
                diffs[x,y] = i
                image[x,y] = val
        #image[(xvalues,yvalues)] += line
#
        if count == 0:
            positions = np.array( [xvalues, yvalues, line])
        else:
            positions = np.hstack((positions,np.array( [xvalues, yvalues, line])))

    # lignes de l'image

        #plt.plot(xvalues,yvalues,"o")
    #image[int(factor*x0),int(factor*y0)] = 0 + 0*1j
    print(image[48])
    plt.imshow(np.sqrt(np.abs(image.real)))
    plt.show()
    #for j in range(geo.nbvox):  # colonnes de l'image
    #    print("working on image column: " + str(j + 1) + "/" + str(geo.nbvox))
    #    for i in range(geo.nbvox):  # lignes de l'image
    #        pass
    
    for j in range(geo.nbvox):  # colonnes de l'image
        print("working on image column: " + str(j + 1) + "/" + str(geo.nbvox))
        for i in range(geo.nbvox):  # lignes de l'image
            # rotation et translation de l'espace
            dist = np.array(distance_point_droite(i, j, x0, y0, angles))
            dist = np.sqrt(1.9) * dist * half_sino / (geo.nbvox)
            index_inf = (np.floor(half_sino + dist)).astype(int)
            index_sup = (np.ceil(half_sino + dist)).astype(int)
            residue = dist % 1
            x = np.arange(0, len(sinogram))
            inf = sinogram[(x, index_inf)]
            sup = sinogram[(x, index_sup)]
            image[j, i] += np.sum(residue * (sup - inf) + inf)
    '''

    # for i in range(geo.nbvox):
    #    slice = image[i, :]
    #    image_ifft[i, :] = np.fft.ifft(slice)

    image_ifft = np.fft.ifft2(((image)))

    #plt.imshow(abs(np.fft.ifftshift(image_ifft)))
    plt.imshow(image_ifft.real, cmap = "Greys")
    plt.show()
    # util.saveImage(abs(image_ifft), "toto")


## main ##
start_time = time.time()
#laminogram()
#backproject()
reconFourierSlice()
#reconFourierSlice()

print("--- %s seconds ---" % (time.time() - start_time))

