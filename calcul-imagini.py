from unittest.main import main
from PIL import Image
import numpy as np
import scipy.ndimage.morphology as morpho
import matplotlib.pyplot as plt
import pandas as pd

def intensitate_medie_membrana (path, x):

    img1 = np.asarray(Image.open(path))

    k = 0
    pixum = 0

    for i in img1:
        for j in i:
            if j[1] >= x:
                pixum += j[1]
                k += 1
    return float(pixum)/k


def rgb2gray(rgb):
    return np.dot(rgb[...,:3], [0.2989, 0.5870, 0.1140])



def contrast(img,a,b,Ta,Tb):
  s=img.shape
  img2=np.zeros((s[0],s[1]))
  img=img.astype('float')
  if len(s)==3 and s[2]==3:
    img_gri=(0.299*img[:,:,0]+0.587*img[:,:,1]+0.114*img[:,:,2])*255
    img_gri=np.clip(img_gri,0,255)
  else:
    img_gri=img
  for i in range(s[0]):
    for j in range(s[1]):
      if img_gri[i,j]<a:
        img2[i,j]=(Ta/a)*img_gri[i,j]
      elif img_gri[i,j]>a and img_gri[i,j]<b:
        img2[i,j]=Ta+((Tb-Ta)/(b-a))*(img_gri[i,j]-a)
      else:
        img2[i,j]=Tb+((255-Tb)/(255-b))*(img_gri[i,j]-b)
  img2=np.clip(img2,0,255)
  img2=img2.astype('uint8')
  return img2


def intensitate_medie_nucleu (path, x, y):

    img1 = np.asarray(Image.open(path))

    k = 0
    pixum = 0

    se = np.ones((1,5))

    img2 = np.zeros(img1.shape)

    for i in range (len(img1)):
        for j in range (len(img1[0])):
            if img1[i][j][2] > x and img1[i][j][1] < y:
                img2[i][j][2] = img1[i][j][2]
    img3 = rgb2gray(img2)
    er_se=morpho.binary_dilation(img3,se)
    contur = er_se - img3 #contur extern
    for i in range(len(contur)):
      for j in range(len(contur[0])):
        if contur[i][j] < 0:
          contur[i][j] = abs(contur[i][j])
        if contur[i][j] > 0:
          k += 1
          pixum += contur[i][j]
    return float(pixum)/k


def intensitate_max_membrana(path, x):

  img1 = np.asarray(Image.open(path))

  maxim = 0
  img = np.zeros(img1.shape)

  for i in range(len(img1)):
    for j in range(len(img1[0])):
      if img1[i][j][1] > x:
        img[i][j] = img1[i][j]
  for i in range(len(img)):
    for j in range(len(img[0])):
      if img[i][j][1] > maxim:
        maxim = img[i][j][1]
  return maxim


def intensitate_max_nucleu(path, x, y):

  img1 = np.asarray(Image.open(path))
  img2 = np.zeros(img1.shape)
  se = np.ones((1,5))
  maxim = 0

  for i in range (len(img1)):
      for j in range (len(img1[0])):
           if img1[i][j][2] > x and img1[i][j][1] < y:
             img2[i][j][2] = img1[i][j][2]
  img3 = rgb2gray(img2)
  er_se=morpho.binary_dilation(img3,se)
  contur = er_se - img3 #contur extern
  for i in range(len(contur)):
    for j in range(len(contur[0])):
      if abs(contur[i][j]) > maxim:
        maxim = abs(contur[i][j])
  return maxim


def rezultate(path, x, y, z):
  d = {"intensitate medie membrana": intensitate_medie_membrana(path, z), 
      "intensitate maxima membrana": intensitate_max_membrana(path, z), 
      "intensitate medie zona perinuceara": intensitate_medie_nucleu(path, x, y), 
      "intensitate maxima zona preinucleara": intensitate_max_nucleu(path, x, y)}
  return pd.Series(d)



if __name__ == "__main__":

    print(rezultate('/home/iustin/ROI5.jpg', 35, 15, 15))