import numpy as np
from scipy.interpolate import interp1d
from numpy import array

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

        # from my simulations
        nuPrime_data = array([8.24036162e-04, 8.24036162e-03, 8.24036162e-02, 8.24036162e-01, 8.24036162e+00, 8.24036162e+01, 8.24036162e+02, 8.24036162e+03])
        Ntheta_data = array([91, 91, 83, 61, 39, 25, 23, 21])
        Nzeta_data = array([59, 59, 57, 51, 45, 37, 27, 21])
        Nxi_data = array([200, 200, 179, 129,  79,  44,  29,  20])
        Nx_data = array([ 6,  6,  7,  7,  9, 12, 16, 16])
        NL_data = array([4, 4, 4, 4, 4, 4, 4, 4])
        
        precond_x_threshold = 5.0
        
    elif configuration=="w7x-m111-pb2.bc":

        # from manual
        # nuPrime_data = np.array([0.001,0.01,0.1,0.3,1.,10.,100.])
        # Ntheta_data = np.array([29,11,11,11,11,11,11])
        # Nzeta_data = np.array([83,64,37,37,37,37,37])
        # Nxi_data = np.array([180,100,37,30,24,13,13])
        # Nx_data = np.array([5,5,5,5,6,7,8])
        # NL_data = np.array([4,4,4,4,4,4,4])
        # from Albert's PoP scan without electrons
        # lots of constant resolution with sudden jumps
        
        nuPrime_data = array([1e-03, 0.99e+00, 1e+00, 1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([35, 35, 35, 35, 35, 35])

        Nzeta_data = array([201, 201, 83, 83, 49, 49])
        Nxi_data = array([220, 220, 90, 90, 36, 36])

        Nx_data = array([6, 6, 8, 8, 10, 10])

        NL_data = array([8, 8, 8, 8, 8, 8])
        
        precond_x_threshold = None
        
        # nuPrime_low = np.array([0.02])
        # nuPrime_hi = np.array([2.0])
        
        # Ntheta = np.array([27])
        # Nzeta = np.array([111])
        # Nxi = np.array([160])
        # Nx = np.array([10])
        # NL = np.array([8])

    elif configuration=="filtered_ref66.bc":
        # Very hi and low collisionality datapoints: 
        #   from Albert's PoP scan without electrons 
        # Other points:
        #   Updated based on experimentation in ~/sfincsSimuls/1/
        nuPrime_data = array([1e-03,3e-02, 0.99e+00, 1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([31,  19, 19, 19, 19, 31, 35, 35, 35])

        Nzeta_data  = array([201, 65, 65, 65, 41, 41, 41, 41, 41])
        Nxi_data    = array([220, 80, 80, 80, 30, 30, 30, 30, 30])
        Nx_data     = array([6,   6,  6,  8,  10, 14, 14, 14, 14])
        NL_data     = array([6,   6,  6,  6,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0

    
    elif configuration=="scaled+filtered_ref66.bc":
        # Very hi and low collisionality datapoints: 
        #   from Albert's PoP scan without electrons 
        # Other points:
        #   Updated based on experimentation in ~/sfincsSimuls/4/
        # NOTE: B field is not that different from un-scaled version
        nuPrime_data = array([1e-03,3e-02, 7e-2,0.99e+00, 1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([31,  19, 19, 19, 19, 19, 31, 35, 35, 35])

        Nzeta_data  = array([201, 65, 45, 45, 45, 31, 31, 31, 31, 31])
        Nxi_data    = array([220, 80, 80, 80, 80, 27, 27, 27, 27, 27])
        Nx_data     = array([6,   6,  11, 6,  8,  10, 14, 14, 14, 14])
        NL_data     = array([6,   6,  4,  6,  6,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0

    elif configuration=="boz_EIM_20190220182940.bc":
        # hi collisionality datapoints, nuPrime >= 1 
        #   from from ref66
        # very hi col: from Albert's PoP col scan without e-
        # nuPrime<1:
        # From resolution scan for this geometry
        nuPrime_data = array([1.17e-3, 3e-03, 1e-02, 6e-02, 6.5e-1,    0.99e+00, 1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([19,  19,  19,  25,  25,   19, 19, 19, 31, 35, 35, 35])
        Nzeta_data  = array([105, 105, 105, 105, 45,   45, 45, 31, 31, 31, 31, 31])
        Nxi_data    = array([140, 140, 123, 80,  60,   60, 60, 27, 27, 27, 27, 27])
        Nx_data     = array([10,  8,   8,   6,   6,    6,  8,  10, 14, 14, 14, 14])
        NL_data     = array([4,   4,   4,   4,   4,    4,  4,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0


    
    elif configuration=="boz_EIM_20190220182950.bc":
        # hi collisionality datapoints, nuPrime >= 1 
        #   from from ref66
        # very hi col: from Albert's PoP col scan without e-
        # nuPrime<1:
        # From resolution scan for boz_EIM_20190220182940.bc and verification
        nuPrime_data = array([1.17e-3, 3e-03, 1e-02, 6e-02, 6.5e-1, 0.95e+00,     1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([19,  19,  19,  25,  25,   21,     19, 19, 31, 35, 35, 35])
        Nzeta_data  = array([105, 105, 105, 105, 45,   45,     45, 31, 31, 31, 31, 31])
        Nxi_data    = array([140, 140, 123, 80,  60,   60,     60, 27, 27, 27, 27, 27])
        Nx_data     = array([6,   6,   6,   6,   6,    6,      8,  10, 14, 14, 14, 14])
        NL_data     = array([4,   4,   4,   4,   4,    4,      4,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0

    elif configuration=="boz_EIM_20190220182951.bc":
        # hi collisionality datapoints, nuPrime >= 1 
        #   from from ref66
        # very hi col: from Albert's PoP col scan without e-
        # nuPrime<1:
        # From resolution scan for boz_EIM_20190220182940.bc
        # and verified by a resolution scan.
        nuPrime_data = array([1.17e-3, 3e-03, 0.02814753, 0.03920992, 0.09330015, 6e-02, 0.30997011, 6.5e-1, 0.95e+00,     1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([19,       19,    25,         25,          25,        25,    25,         25,   21,     19, 19, 31, 35, 35, 35])
        Nzeta_data  = array([105,      105,   101,        101,         95,        101,   69,         45,   45,     45, 31, 31, 31, 31, 31])
        Nxi_data    = array([140,      140,   124,        124,         77,        80,    68,         60,   60,     60, 27, 27, 27, 27, 27])
        Nx_data     = array([8,        8,     8,          8,           6,         6,     6,          8,    8,      8,  10, 14, 14, 14, 14])
        NL_data     = array([4,        4,     4,          4,           4,         4,     4,          4,    4,      4,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0


    elif configuration=="boz_EIM_20190220182952.bc":
        # hi collisionality datapoints, nuPrime >= 1 
        #   from from ref66
        # very hi col: from Albert's PoP col scan without e-
        # nuPrime<1:
        # From resolution scan for boz_EIM_20190220182951.bc
        # NEEDS VERIFICATION
        nuPrime_data = array([1.17e-3, 3e-03, 0.02814753, 0.03920992, 0.09330015, 6e-02, 0.30997011, 6.5e-1, 0.95e+00,     1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([19,       19,    25,         25,          25,        25,    25,         25,   21,     19, 19, 31, 35, 35, 35])
        Nzeta_data  = array([105,      105,   101,        101,         95,        101,   69,         45,   45,     45, 31, 31, 31, 31, 31])
        Nxi_data    = array([140,      140,   124,        124,         77,        80,    68,         60,   60,     60, 27, 27, 27, 27, 27])
        Nx_data     = array([8,        8,     8,          8,           6,         6,     6,          8,    8,      8,  10, 14, 14, 14, 14])
        NL_data     = array([4,        4,     4,          4,           4,         4,     4,          4,    4,      4,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0

    elif configuration=="boz_EIM_20190220182954.bc":
        # hi collisionality datapoints, nuPrime >= 1 
        #   from from ref66
        # very hi col: from Albert's PoP col scan without e-
        # nuPrime<1:
        # From resolution scan for boz_EIM_20190220182951.bc
        # NEEDS VERIFICATION
        nuPrime_data = array([1.17e-3, 3e-03, 0.02814753, 0.03920992, 0.09330015, 6e-02, 0.30997011, 6.5e-1, 0.95e+00,     1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([19,       19,    25,         25,          25,        25,    25,         25,   21,     19, 19, 31, 35, 35, 35])
        Nzeta_data  = array([105,      105,   101,        101,         95,        101,   69,         45,   45,     45, 31, 31, 31, 31, 31])
        Nxi_data    = array([140,      140,   124,        124,         77,        80,    68,         60,   60,     60, 27, 27, 27, 27, 27])
        Nx_data     = array([8,        8,     8,          8,           6,         6,     6,          8,    8,      8,  10, 14, 14, 14, 14])
        NL_data     = array([4,        4,     4,          4,           4,         4,     4,          4,    4,      4,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0

    elif configuration=="boz_EIM_20190220182932.bc":
        # hi collisionality datapoints, nuPrime >= 1 
        #   from from ref66
        # very hi col: from Albert's PoP col scan without e-
        # nuPrime<1:
        # From resolution scan for boz_EIM_20190220182951.bc
        # NEEDS VERIFICATION
        nuPrime_data = array([1.17e-3, 3e-03, 0.02814753, 0.03920992, 0.09330015, 6e-02, 0.30997011, 6.5e-1, 0.95e+00,     1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([19,       19,    25,         25,          25,        25,    25,         25,   21,     19, 19, 31, 35, 35, 35])
        Nzeta_data  = array([105,      105,   101,        101,         95,        101,   69,         45,   45,     45, 31, 31, 31, 31, 31])
        Nxi_data    = array([140,      140,   124,        124,         77,        80,    68,         60,   60,     60, 27, 27, 27, 27, 27])
        Nx_data     = array([8,        8,     8,          8,           6,         6,     6,          8,    8,      8,  10, 14, 14, 14, 14])
        NL_data     = array([4,        4,     4,          4,           4,         4,     4,          4,    4,      4,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0

    elif configuration=="boz_EIM_20190220182941.bc":
        # hi collisionality datapoints, nuPrime >= 1 
        #   from from ref66
        # very hi col: from Albert's PoP col scan without e-
        # nuPrime<1:
        # From resolution scan for boz_EIM_20190220182951.bc
        # NEEDS VERIFICATION
        nuPrime_data = array([1.17e-3, 3e-03, 0.02814753, 0.03920992, 0.09330015, 6e-02, 0.30997011, 6.5e-1, 0.95e+00,     1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([19,       19,    25,         25,          25,        25,    25,         25,   21,     19, 19, 31, 35, 35, 35])
        Nzeta_data  = array([105,      105,   101,        101,         95,        101,   69,         45,   45,     45, 31, 31, 31, 31, 31])
        Nxi_data    = array([140,      140,   124,        124,         77,        80,    68,         60,   60,     60, 27, 27, 27, 27, 27])
        Nx_data     = array([8,        8,     8,          8,           6,         6,     6,          8,    8,      8,  10, 14, 14, 14, 14])
        NL_data     = array([4,        4,     4,          4,           4,         4,     4,          4,    4,      4,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0
    
    elif configuration=="boz_EIM_20190220182939.bc":
        # hi collisionality datapoints, nuPrime >= 1 
        #   from from ref66
        # very hi col: from Albert's PoP col scan without e-
        # nuPrime<1:
        # From resolution scan for boz_EIM_20190220182940.bc and boz_EIM_20190220182950.bc
        nuPrime_data = array([1.17e-3, 3e-03, 1e-02, 6e-02, 6.5e-1, 0.95e+00, 1.49+00,     2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([19,  19,  19,  25,  25,   21, 19,    19, 31, 35, 35, 35])
        Nzeta_data  = array([105, 105, 105, 105, 45,   45, 45,    31, 31, 31, 31, 31])
        Nxi_data    = array([140, 140, 123, 80,  60,   60, 60,    27, 27, 27, 27, 27])
        Nx_data     = array([10,  8,   6,   6,   6,    6,  8,     10, 14, 14, 14, 14])
        NL_data     = array([4,   4,   4,   4,   4,    4,  4,     4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0
    
    elif configuration=="scaled_ref172.bc":
        # hi collisionality datapoints, nuPrime >= 1 
        #   from from ref66
        # very hi col: from Albert's PoP col scan without e-
        # nuPrime<1:
        # based on "boz_EIM_20190220182951.bc":
        # NOT VERIFIED
        nuPrime_data = array([1.17e-3, 3e-03, 0.02814753, 0.03920992, 0.09330015, 6e-02, 0.30997011, 6.5e-1, 0.95e+00,     1e+00, 2.5e+00, 80 ,1.49e+02, 1.5e+02, 1.00000000e+04])
        
        Ntheta_data = array([19,       19,    25,         25,          25,        25,    25,         25,   21,     19, 19, 31, 35, 35, 35])
        Nzeta_data  = array([105,      105,   101,        101,         95,        101,   69,         45,   45,     45, 31, 31, 31, 31, 31])
        Nxi_data    = array([140,      140,   124,        124,         77,        80,    68,         60,   60,     60, 27, 27, 27, 27, 27])
        Nx_data     = array([8,        8,     8,          8,           6,         6,     6,          8,    8,      8,  10, 14, 14, 14, 14])
        NL_data     = array([4,        4,     4,          4,           4,         4,     4,          4,    4,      4,  4,  4,  4,  4,  4])
        
        precond_x_threshold = 10.0




    elif configuration=="real_tj20_jmc.bc":
        nuPrime_data = np.array([2.27e-3, 2.27e-2, 2.27e-1, 2.27e+0, 2.27e+1, 2.27e+2, 2.27e+3])
        Ntheta_data = np.array([61, 61, 51, 41, 23, 21, 21])
        Nzeta_data = array([111, 111,  71,  51,  23,  21, 21])
        Nxi_data = array([160, 160,  120,  100,  40,  29, 20])
        Nx_data = array([ 6,  6,  8,  8, 12, 16, 16])
        NL_data = array([4, 4, 4, 4, 4, 4, 4])        
        precond_x_threshold = 5.0

        
    else:
        raise ValueError("Unrecognized configuration: '" + configuration + "'")

    Ntheta_interp = interp1d(np.log10(nuPrime_data),Ntheta_data,fill_value = (Ntheta_data[0],Ntheta_data[-1]),bounds_error=False)
    Nzeta_interp = interp1d(np.log10(nuPrime_data),Nzeta_data,fill_value = (Nzeta_data[0],Nzeta_data[-1]),bounds_error=False)
    Nxi_interp = interp1d(np.log10(nuPrime_data),Nxi_data,fill_value = (Nxi_data[0],Nxi_data[-1]),bounds_error=False)
    Nx_interp = interp1d(np.log10(nuPrime_data),Nx_data,fill_value = (Nx_data[0],Nx_data[-1]),bounds_error=False)
    NL_interp = interp1d(np.log10(nuPrime_data),NL_data,fill_value = (NL_data[0],NL_data[-1]),bounds_error=False)

    def f(log10NuPrime):
        if log10NuPrime > np.log10(precond_x_threshold):
            return 0
        else:
            return 1
    precond_x = np.vectorize(f,otypes=[int])
    
    return Ntheta_interp,Nzeta_interp,Nxi_interp,Nx_interp,NL_interp,precond_x
