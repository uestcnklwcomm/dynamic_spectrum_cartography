from tsvd import  *
import  numpy as np
import  time
import  matplotlib.pyplot as plt
import h5py
from scipy.io import  loadmat
from scipy.io import  savemat

#loading spectrum data, size I * J * K * T
RTFiles = 'RayTracingMaps.mat'
SamplingMask = 'SamplingMask.mat'

spectrum_data = loadmat(RTFiles)
sampling_mask = loadmat(SamplingMask)



Spectrum_tens = np.array(spectrum_data['X4DT'])
Spectrum_tens = Spectrum_tens.transpose([0,3,1,2]) # Transpose to I*T*J*K

mask = np.array(sampling_mask['Wtens'])
mask = mask.transpose([0,2,1])


n1, n2, n3, n4 = Spectrum_tens.shape

Xcube = Spectrum_tens[:, :, :, 6]

X = Tensor(Xcube)
Xfrob = tfrobnorm(X)


np.random.seed(0)

tube = False
#data missing rate, i.e., sampling ratio = 1 - rho







#Tensor SVD
################################ TOUCAN #######################################
rank = 1
toucan_err_curve = np.zeros(n2)

# p = 0.9
#
# if (tube is False):
#     mask = np.random.rand(n1, n2, n3)
#     mask[mask > p] = 1
#     mask[mask <= p] = 0
#     mask = mask.astype(int)
# else:
#     mask = np.random.rand(n1, n2)
#     mask[mask > p] = 1
#     mask[mask <= p] = 0
#     mask = mask.astype(int)
#     mask = np.repeat(mask[:, :, np.newaxis], n3, axis=2)

sig = 0

# fun = lambda L: [0, tfrobnorm(L - X) / Xfrob]
fun = lambda Q,k: [0, tfrobnorm_array(Q.array()[:,k,:] - X.array()[:,k,:]) / tfrobnorm_array(X.array()[:,k,:])]
# fun = lambda Q,k : [0,0]
Y_hat_assem = np.zeros([n1,n2,n3,n4])

for freq in range(n4):
    Xcube_freq = Spectrum_tens[:, :, :, freq]
    Y = Tensor(Xcube_freq * mask)
    Y_hat_toucan, U_toucan,stats_toucan, tElapsed_toucan = toucan(Y,mask,rank,tube=False,outer=1,cgtol=1e-6,mode='online',fun=fun,
                                                     randomOrder=False,verbose=False)
    Y_hat_assem[:,:,:,freq] = Y_hat_toucan.array()

Y_hat_assem = Y_hat_assem.transpose([0,2,3,1])
mdicXhat = {"X4DHatTOUCAN": Y_hat_assem}
savemat("X4DHatTOUCAN.mat",mdicXhat)

