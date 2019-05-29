q# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 17:31:35 2019

@author: Maple
"""

import numpy as np
import matplotlib.pyplot as plt
import os 
from PIL import Image as im

# Set up Data

d = 'CroppedYale'   # Main folder Name
students = os.listdir(d)    # Lists directories in folder d

# Initialize matrix that will store all student images as columns
pic = im.open(os.listdir(d + '\\' + students[0])[0])    
pnum = len(students)*len(os.listdir(d + '\\' + students[0]))    # Number of pictures total for all students
fw = pic.size[0]        # picture width
fh = pic.size[1]        # picture height
colength = fw*fh

studentall_fft = np.zeros([colength,pnum])      # FFT Filtered Images
studentall = np.zeros([colength,pnum])          # Unmodified Images

# Create filter for FFT Filter
x = np.arange(0, fw)
y = np.arange(0, fh)

X, Y = np.meshgrid(x,y)
fil = 1 - np.exp(-0.01*((X - fw / 2) ** 2 + (Y - fh / 2) ** 2))

k = 0;
for s in students:
    sf = d + '\\' + s      # Name of specific student folder
    files = os.listdir(sf)  # Obtains all files in student folder    
    for f in files:
        # Reshapes image into column and add to matrix
        I = im.open(sf+'\\'+f)
        col = np.reshape(I, [colength,1])
        studentall[:, k] = col[:,0]
        
        # Applies FFT Filter onto Images and stores in matrix
        It = np.fft.fft2(I)
        Its = np.fft.fftshift(It)
        Itsf = Its*fil
        Itf = np.fft.ifftshift(Itsf)
        If = np.fft.ifft2(Itf) 
        col = np.reshape(abs(If), [colength,1])
        studentall_fft[:, k] = col[:,0]
        
        k = k+1
# %% Preview Images as necessary
# plt.imshow(np.reshape(studentall[:,200], [fh, fw]), cmap = 'gray')
# plt.imshow(np.reshape(studentall_fft[:,200], [fh, fw]), cmap = 'gray')        
     
#%% Receives a data set and returns the mean centered version of the data
studentall_ori = studentall
def center_mean(data):
    means = np.mean(data, axis = 0)
    data = data-means
    data = data*(1/(len(data)-1))
    return data

studentall = center_mean(studentall)
studentall_fft = center_mean(studentall_fft)

# %% Apply SVD - note that in numpy, v is returned transposed
u, s, v = np.linalg.svd(studentall, full_matrices=False, compute_uv=True)

# %% Apply SVD to FFT Filtered data
uf, sf, vf = np.linalg.svd(studentall_fft, full_matrices=False, compute_uv=True)

# %% Function receives a sigma vector of singular values and 
#       plots the variances in each PC. Then it titles the graph
#       with the inputed title and saves it as the inputed file name
def plot_pcvars(sigma, file_name, title):
    plt.plot(np.arange(1,pnum+1), (sigma**2)/np.sum(sigma**2), 'ro')
    plt.title(title)
    plt.ylim([0,1])
    plt.xlabel('Principal Component Number')
    plt.ylabel('Variance')
    plt.savefig(file_name)

plot_pcvars(s, 'pc_var.png', 'Variance in Principal Components')
# %%
plot_pcvars(sf, 'pc_varz_fft.png', 'Variance in Principal Components with FFT')

# %% Receives the components of the SVD decomposition and returns
#       how many principal components contain 90% of the variance
def get_rank(u,s,v):
    var_tot = 0; 
    rank = 0; 
    while var_tot < 0.90:
        var_tot = var_tot + (s[rank]**2) / np.sum(s**2)
        rank = rank+1
    return rank

rank = get_rank(u,s,v)
rank_fft = get_rank(uf,sf,vf)

# %% Receives u,s,v, the components of the SVD decomposition and returns the
#        r-rank version of the original data
def recon_r_rank(rank,u,s,v): 
    rank = rank
    u_rec = u[:, 0:rank]
    s_rec = np.diag(s)[0:rank, 0:rank]
    v_rec = v[0:rank,:]
    # Reconstruct data with r-rank versions of U,S,V
    x_rec = np.matmul(np.matmul(u_rec, s_rec), v_rec)
    return x_rec

x_rec = recon_r_rank(rank,u,s,v)
x_rec_fft = recon_r_rank(rank_fft,uf,sf,vf)

# %% Create all Visualizations
faces = np.arange(0,9)*64
pf = np.concatenate((studentall_ori[:,faces], x_rec[:,faces], x_rec_fft[:,faces]),axis=1)
fig,ax = plt.subplots(3,9, figsize=(12,6))
xl = np.shape(ax)[0]
yl = np.shape(ax)[1]
i = 0
j = 0
k = 0
ax[2][1].set_xlabel('Original Faces')
ax[2][4].set_xlabel('Reconstructed Faces')
ax[2][7].set_xlabel('Reconstructed FFT Filtered Faces')
plt.suptitle('Visualization of Different Face Experiments', ha='center', fontsize=16)
for p in np.arange(0,xl*yl):
    pic = np.reshape(pf[:,p], [fh, fw])
    ax[i][j].imshow(pic, cmap='gray')
    ax[i][j].set_title('Face' + str(k+1))
    ax[i][j].set_xticklabels([])
    ax[i][j].set_yticklabels([])
    i = i+1
    k=k+1
    if i == xl:
        i = 0
        j = j+1
    if k == yl:
        k = 0
#plt.savefig('all_faces.png')
plt.show()

#%% 
fig, ax = plt.subplots(2,6, figsize=(8,4))
i = 0
j = 0
pcs = u[:,0:6]
plt.suptitle('Visualization of Dominate Features', ha='center', fontsize=16)
ax[1][1].set_xlabel('Unmodified Faces')
ax[1][4].set_xlabel('FFT Filtered Faces')
for p in np.arange(0,6):
    pcim = np.reshape(pcs[:,p], [fh, fw])
    ax[i][j].imshow(pcim, cmap='gray')
    ax[i][j].set_title('Feature' + str(p+1))
    ax[i][j].set_xticklabels([])
    ax[i][j].set_yticklabels([])
    j = j+1
    if j == 3:
        j = 0
        i = i+1

i = 0
j = 3
pcs = uf[:,0:6]
for p in np.arange(0,6):
    pcim = np.reshape(pcs[:,p], [fh, fw])
    ax[i][j].imshow(pcim, cmap='gray')
    ax[i][j].set_title('Feature' + str(p+1))
    ax[i][j].set_xticklabels([])
    ax[i][j].set_yticklabels([])
    j = j+1
    if j == 6:
        j = 3
        i = i+1
plt.savefig('dom_faces.png')    
plt.show()

