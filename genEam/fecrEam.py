import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

'''
Interatomic Potential to Simulate
Radiation Damage in Fe-Cr Alloys
G. Bonny, R.C. Pasianot, D. Terentyev and L. Malerba
'''

plot = False


def fcut(X, Xcut):
    if X <= 1.0:
        return 1.0
    if X > 1 and X < Xcut:
        return 1.0 - np.power(X-1, 3)/np.power(Xcut-1, 3)
    if X > Xcut:
        return 0.0


def densityCrCr(r):
    a0 = 2.878
    r0 = a0 * 0.866
    beta = 5.0
    Xcut = 1.65 #* r0
    des0 = 0.0676504617
    X = r / r0


    des1 = np.exp(-beta * X) / X
    des2 = np.exp(-beta * Xcut) / Xcut

    des =  (des1 - des2) * fcut(X, Xcut) / des0

    return des

def VCrCr(r):
    rk_n = [4.112494835E+00,3.738631668E+00,3.364768501E+00,2.990905335E+00,\
        2.617042168E+00]

    ak_n = [-6.079639840E-02, -8.216224618E-01, 2.424000178E+00, -1.275324968E+00,\
        7.942650649E-01]

    rk_n = np.flip(np.asarray(rk_n), 0)
    ak_n = np.flip(np.asarray(ak_n), 0)

    a = ak_n * (rk_n -r) ** 3
    return  np.sum((rk_n > r) *a)

rk = [0.5,0.6905,0.8811,1.0716,1.2621,1.4526,\
1.6432,1.8337,2.0242,2.2147,2.4053,2.5958,\
2.7863,2.9768,3.1674,3.3579,3.5484,3.7389,\
3.9295,4.12]

Vcrcr_rk=[14.0403,10.8094,8.0938,5.8496,4.0328,\
2.5993,1.5053,0.7065,0.1592,-0.1809,\
-0.3576,-0.4149,-0.3932,-0.3108,-0.1859,\
-0.0715,-0.0166,-0.0032,-0.0004,0]


rho_rk=[17.5647,8.6776,4.64,2.6025,1.5071,0.8929,\
0.538,0.3284,0.2024,0.1257,0.0784,0.049,\
0.0304,0.0185,0.0108,0.0059,0.0029,0.0011,\
0.0002,0]

rho=[0,0.1579,0.3158,0.4737,0.6316,0.7895,\
0.9474,1.1053,1.2632,1.4211,1.5789,\
1.7368,1.8947,2.0526,2.2105,2.3684,\
2.5263,2.6842,2.8421,3]

F_rho=[-0.0004, -1.196, -1.3247, -1.4043, -1.4346,\
-1.4158, -1.3477, -1.2304, -1.064, -0.8483,\
-0.5834, -0.2693, 0.0939, 0.5064, 0.9681,\
1.479, 2.0391, 2.6484, 3.3069, 4.0146]

#xnew = np.linspace(0, rk[len(rk)-1], num=20, endpoint=True)
for i,r in enumerate(rk):
  print r, VCrCr(r),rho_rk[i]/(densityCrCr(r) + .0001)


F = interp1d(rho, F_rho,  kind='cubic') #embed
phi = interp1d(rk, rho_rk,  kind='cubic') # potUU
Vcrcr = interp1d(rk, Vcrcr_rk,  kind='cubic') # potUU

if plot:
    plt.plot(rk, Vcrcr_rk) 
    plt.show()

    xnew = np.linspace(rk[0], rk[len(rk)-1], num=41, endpoint=True)
    plt.plot(xnew, phi(xnew), '--')
    plt.plot(rk, rho_rk) 
    plt.show()

    xnew = np.linspace(rho[0], rho[len(rho)-1], num=41, endpoint=True)
    plt.plot(rho, F_rho, 'X')
    plt.plot(xnew, F(xnew), '--')
    plt.grid(True)
    plt.show()
