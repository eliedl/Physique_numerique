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
        print("working on image column: "+str(j+1)+"/"+str(geo.nbvox))
        for i in range(geo.nbvox): # lignes de l'image
                #rotation et translation de l'espace
                dist = np.array(distance_point_droite(i, j, x0, y0, angles))
                dist = np.sqrt(1.9) * dist * half_sino / (geo.nbvox)
                index_inf = (half_sino + dist).astype(int)
                x = np.arange(0, len(sinogram))
                image[j, i] += np.sum(sinogram[(x, index_inf)])


    util.saveImage(image, "lam")


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
        print("working on image column: "+str(j+1)+"/"+str(geo.nbvox))
        for i in range(geo.nbvox): # lignes de l'image
            dist = np.array(distance_point_droite(i, j, x0, y0, angles))
            dist = np.sqrt(1.9) * dist * half_sino / (geo.nbvox)
            index_inf = (half_sino + dist).astype(int)
            x = np.arange(0,len(sinogram))
            image[j, i] += np.sum(sinogram[(x,index_inf)])

    util.saveImage(image, "fbp")


## reconstruire une image TDM en mode retroprojection
def reconFourierSlice():

    sinogram = np.loadtxt('data/sinogram-patient.txt')

    angles = np.loadtxt('data/angles.txt')

    nbprj = len(angles)

    # initialiser une image reconstruite, complexe
    # pour qu'elle puisse contenir sa version FFT d'abord
    image = np.zeros((geo.nbvox, geo.nbvox), 'complex')

    #image reconstruite
    image = np.zeros((geo.nbvox, geo.nbvox), 'complex')

    fft_sino = np.zeros(sinogram.shape, 'complex')
    print(sinogram.shape)
    for i in range(len(sinogram)):
        line = sinogram[i, :]
        fftline = np.fft.fft(line)
        fftshift = np.fft.fftshift(fftline)
        freq = np.fft.fftfreq(len(line), 1)
        freqshift = np.fft.fftshift(freq)
        fftline = np.conj(np.abs(freqshift) * fftshift)
        fft_sino[i, :] = fftline

    #plt.imshow(sinogram)
    #plt.show()

    plt.imshow(abs(fft_sino))
    plt.show()

    x0 = y0 = geo.nbvox / 2
    half_sino = sinogram.shape[1] / 2
    for j in range(geo.nbvox):  # colonnes de l'image
        print("working on image column: " + str(j + 1) + "/" + str(geo.nbvox))
        for i in range(geo.nbvox):  # lignes de l'image
            # rotation et translation de l'espace
            dist = np.array(distance_point_droite(i, j, x0, y0, angles))
            #dist = np.sqrt(1.9) * dist * half_sino / (geo.nbvox)
            index_inf = (half_sino + dist).astype(int)
            x = np.arange(0, len(sinogram))
            image[j, i] += np.sum(fft_sino[(x, index_inf)])

    plt.imshow(abs(image))
    plt.show()
    image_ifft = np.fft.fft2(image)
    plt.imshow(abs(image_ifft))
    plt.show()
    #util.saveImage(abs(image_ifft), "toto")


## main ##
start_time = time.time()
#laminogram()
reconFourierSlice()
#reconFourierSlice()
print("--- %s seconds ---" % (time.time() - start_time))

