import numpy as np
from scipy.interpolate import interp1d

def resolution_from_collisionality(configuration):

    configuration = configuration.rsplit("/",1)[-1]
    if configuration=="lhd113208t4640.bc":

        # from manual
        # nuPrime_data = np.array([0.001,0.01,0.1,0.3,1.,10.,100.])
        # Ntheta_data = np.array([85,21,15,15,15,15,15])
        # Nzeta_data = np.array([41,25,13,13,13,13,13])
        # Nxi_data = np.array([103,70,37,34,13,13,13])
        # Nx_data = np.array([5,5,5,6,8,8,8])
        # NL_data = np.array([4,4,4,4,4,4,4])

        # from my anus
        nuPrime_data = np.array([0.001,100.])
        Ntheta_data = np.array([91,21])
        Nzeta_data = np.array([59,21])
        Nxi_data = np.array([200,20])
        Nx_data = np.array([6,8])
        NL_data = np.array([4,4])
        

        # nuPrime_low = np.array([0.001])
        # nuPrime_hi = np.array([0.05])
        
        # Ntheta = np.array([61])
        # Nzeta = np.array([59])
        # Nxi = np.array([200])
        # Nx = np.array([6])
        # NL = np.array([4])
    elif configuration=="w7x-m111-pb2.bc":

        # from manual
        # nuPrime_data = np.array([0.001,0.01,0.1,0.3,1.,10.,100.])
        # Ntheta_data = np.array([29,11,11,11,11,11,11])
        # Nzeta_data = np.array([83,64,37,37,37,37,37])
        # Nxi_data = np.array([180,100,37,30,24,13,13])
        # Nx_data = np.array([5,5,5,5,6,7,8])
        # NL_data = np.array([4,4,4,4,4,4,4])
        # from my anus
        nuPrime_data = np.array([0.001,100.])
        Ntheta_data = np.array([41,21])
        Nzeta_data = np.array([111,21])
        Nxi_data = np.array([200,20])
        Nx_data = np.array([9,12])
        NL_data = np.array([8,8])

        # nuPrime_low = np.array([0.02])
        # nuPrime_hi = np.array([2.0])
        
        # Ntheta = np.array([27])
        # Nzeta = np.array([111])
        # Nxi = np.array([160])
        # Nx = np.array([10])
        # NL = np.array([8])
    else:
        raise ValueError("Unrecognized configuration: '" + configuration + "'")

    Ntheta_interp = interp1d(np.log(nuPrime_data),Ntheta_data,fill_value = (Ntheta_data[0],Ntheta_data[-1]),bounds_error=False)
    Nzeta_interp = interp1d(np.log(nuPrime_data),Nzeta_data,fill_value = (Nzeta_data[0],Nzeta_data[-1]),bounds_error=False)
    Nxi_interp = interp1d(np.log(nuPrime_data),Nxi_data,fill_value = (Nxi_data[0],Nxi_data[-1]),bounds_error=False)
    Nx_interp = interp1d(np.log(nuPrime_data),Nx_data,fill_value = (Nx_data[0],Nx_data[-1]),bounds_error=False)
    NL_interp = interp1d(np.log(nuPrime_data),NL_data,fill_value = (NL_data[0],NL_data[-1]),bounds_error=False)

    return Ntheta_interp,Nzeta_interp,Nxi_interp,Nx_interp,NL_interp
