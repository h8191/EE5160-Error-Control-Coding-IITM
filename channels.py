import random
from scipy.special import erfc

def qfunc_custom(x):
    return 0.5 * erfc(x/pow(2,0.5))

'''channels'''
class BSC(object):
    def __init__(self,p):
        self.p = p

    def transmit(self,x,x_len):
        error, mask = 0, 1
        for i in range(x_len):
            if random.random() < self.p:
                error ^= mask
            mask <<= 1
        return x ^ error 

'''modulation schemes'''
class BPSK(object):
    def __init__(self, Eb):
        self.Eb = Eb
        self.d  = pow(Eb, 0.5)

    def modulate(self, bit_list):
        return [self.d * (1 - 2*bit) for bit in bit_list]
    
    def demodulate(self, val_list):
        return [0 if val >= 0 else 1 for val in val_list]

