# -*- coding: utf-8 -*-
"""
Created on Sun Mar 10 19:54:22 2019

@author: Maple
"""

import imageio
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
from skimage import color

# %% Load data
# Necessary File Names and Variables for Different Trials
# stuff = ['penguin.mp4', 'varp.png', 'evalsp.png', 'penguin1.png', 'framediffp.png', 60]
stuff = ['spidey2.wmv', 'vars.png', 'evalss.png', 'spidey1.png', 'framediffs.png', 120]
file = stuff[0]
vid = imageio.get_reader(file, 'ffmpeg')
dt = 1/60       # 1/Framerate
k = 0
for n, frame in enumerate(vid):
    k = k+1
hei = np.shape(frame)[0] 
wid = np.shape(frame)[1]
croph = stuff[5]
cropw = 0
hnew=(hei-2*croph)
wnew = (wid-2*cropw)
hei = hnew
wid = wnew
X = np.zeros([hei*wid, k])

# Reshape and Store into Matrix
for n, frame in enumerate(vid):
    f = color.rgb2gray(frame)
    fc = f[croph:-croph,:]
    X[:,n] = np.reshape(fc, [hei*wid, 1])[:,0]

plt.imshow(fc, cmap='gray')
plt.show()

#%% Preview Frames as necessary
#plt.figure(figsize=(5,4))
#plt.imshow(color.rgb2gray(frame), cmap='gray')
#plt.axis('off')
#plt.savefig('uncropped.png')
#
#
##%%
#plt.figure(figsize=(5,4))
#plt.imshow(np.reshape(X[:,6], [hei,wid]),cmap='gray')
#plt.axis('off')
#plt.savefig('cropped.png')


#%%
endp = np.shape(X)[1]
X1 = X[:, np.arange(0, endp - 1)]
X2 = X[:, np.arange(1, endp)]

#%%
u,s,vh = np.linalg.svd(X1, full_matrices=False)
v = np.transpose(vh)


#%% Plot Singular Values
sl = np.shape(s)[0]
plt.scatter(np.arange(0,sl), s**2 / sum(s**2), color = 'red')
plt.xlabel('Principal Component')
plt.ylabel('Variance')
plt.title('Amount of Variance in each Principal Component')
plt.ylim([0, 1])
plt.savefig(stuff[1])
plt.show()


# %%
r=1
u_r = u[:, 0:r]
s_r = np.diag(s[0:r])
v_r = v[:, 0:r]
Atilde = np.transpose(u_r)@X2@v_r@la.inv(s_r)

# %%
w, d = la.eig(Atilde)
Phi = X2@v_r@la.inv(s_r)@d

# %% Plot Eigenvalues
omega = np.log(w)/dt
x = np.real(omega)
y = np.imag(omega)
plt.figure(figsize=(8,5))
plt.axhline(color='black')
plt.axvline(color='black')  
plt.title('Plot of Eigenvalues of $\widetilde{A}$ for ' + stuff[0][0].upper() + stuff[0][1:-4] + ' video')
plt.xlabel('Real')
plt.ylabel('Im')
plt.scatter(x,y, zorder=3)
plt.xlim([-3, 3])
plt.grid()
plt.savefig(stuff[2])

#%% Solve for b 
x0 = X1[:,0:1]
b = la.lstsq(Phi, x0, rcond=None)[0]

# %% Create DMD Modes
mm1 = np.shape(X1)[1]
time_dyn = np.zeros([r, mm1], dtype=complex)
t = np.arange(0, mm1)*dt
for tt in np.arange(len(t)):
    time_dyn[:, tt] = (b[:,0]*np.exp(omega*t[tt]))


Xdmd = Phi@time_dyn

#%%
res = Xdmd - abs(Xdmd)
X_sparse = X1 - abs(Xdmd)


#%% Plot Forground and Background Frames against Original
n = 50
frp_act = np.reshape(X[:, n], [hei,wid])
frp_sparse = np.reshape(Xdmd[:, n], [hei,wid])
frp_dmd = np.reshape(X_sparse[:, n], [hei,wid])
frp_res = np.reshape(res[:, n], [hei,wid])

fig,ax = plt.subplots(3,1, figsize=(6,8))
ax[0].imshow(frp_act, cmap='gray')
ax[0].axis('off')
ax[0].set_title('Original Frame')
ax[2].imshow(abs(frp_sparse - frp_res), cmap='gray')
ax[2].axis('off')
ax[2].set_title('Low Rank DMD of Frame (Background)')
ax[1].imshow(255-abs(frp_dmd), cmap='gray')
ax[1].axis('off')
ax[1].set_title('Sparse DMD of Frame (Forground)')
plt.savefig(stuff[3])
plt.show()

#%% Plot Different Frames
n = 40
frp_beg = np.reshape(X[:,4], [hei,wid])
frp_end = np.reshape(X[:,n], [hei, wid])

fig,axs = plt.subplots(2,1)
axs[0].imshow(frp_beg, cmap='gray')
axs[0].axis('off')
axs[0].set_title('Frame 4')
axs[1].imshow(frp_end, cmap='gray')
axs[1].axis('off')
axs[1].set_title('Frame ' + str(n))
plt.savefig(stuff[4])









