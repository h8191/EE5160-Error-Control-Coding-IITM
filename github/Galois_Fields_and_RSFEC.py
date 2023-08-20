import numpy as np
import matplotlib.pyplot as plt
import time
import math

time_stack = []
tic = lambda: time_stack.append(time.time())

def toc():
    if len(time_stack) > 0:
        print('Time elapsed := ',time.time() - time_stack.pop())

def hamming_distance(x,y):
    return sum([i!=j for i,j in zip(x,y)])

class GF(object):
    def __init__(self,m=3,prim_poly=0b1011):
        self.m = m
        self.prim_poly = prim_poly
        self.size = 1<<m
        self.sizem1 = self.size - 1#size minus 1
        self.gf_exp = [0 for i in range(2*self.size)]#input i output alpha pow i 
        self.gf_log = [0 for i in range(self.size)]#input alpha pow i output i 

        tmplog = 1
        for tmpexp in range(self.size-1):# iterate over vals of exponent
            self.gf_exp[tmpexp] = tmplog
            self.gf_exp[tmpexp + self.sizem1] = tmplog
            self.gf_log[tmplog] = tmpexp
            tmplog = tmplog << 1
            if tmplog & self.size:
                tmplog = tmplog ^ self.prim_poly #xor to keep it within GF

    
    def add(self,x,y):
        return x^y #addition is xor since characteristic = p is 2

    def mult(self,x,y):
        if x==0 or y==0:
            return 0
        #modulo op to keep the exponent less than size
        #return self.gf_exp[(self.gf_log[x]+self.gf_log[y])%self.sizem1]
        #or use twice as much memory to avoid modulo opoeration
        return self.gf_exp[self.gf_log[x]+self.gf_log[y]]

    def pow(self,x,p):
        if x==0:
            return 0
        return self.gf_exp[(self.gf_log[x]*p)%self.sizem1]
    
    def poly_eval(self,poly,x):
        #left side of poly list is MSB and right side is LSB
        val  = 0
        xpow = 1#for mutlipying with constant of poly
        for ind in range(len(poly)-1,-1,-1):
            tmp = self.mult(poly[ind], xpow)#multiply coeff and power of x
            val = self.add(val, tmp)  #add it to previous terms in polynomial
            xpow = self.mult(xpow, x) #raise the power of x
        return val
    
    def log_poly_roots(self, poly):
        #find i for non zero roots(alpha power i)
        #if poly[-1] == 0:
            #if LSB of poly is zero then zero is a root
        roots = [i for i in range(self.sizem1) if self.poly_eval(poly,self.gf_exp[i])==0]
        return roots

    def poly_roots(self, poly):
        roots = [val for val in range(self.sizem1) if self.poly_eval(poly,val)==0]
        return roots

    def print(self):
        print("poly \t power \t power \t poly")
        print("0 \t -inf \t -inf \t 0")
        for i in range(1,self.size):
            print("{} \t {} \t {} \t {}".format(i,self.gf_log[i], \
                                                i - 1,self.gf_exp[i - 1]) )


class reed_solomon(object):
    def __init__(self, n,d,prim_poly = 0b1011,m=None):
        self.n = n
        self.d = d
        self.k = n - d + 1
        if m is None:
            self.m = int(math.ceil(math.log2(n + 1)))
        else:
            self.m = m

        self.two_power_m = 1 << self.m
        self.GF = GF(self.m,prim_poly)
        self.gen_poly = self.calc_generator_poly()
    
    def calc_generator_poly(self):
        tmp = [1]
        for i in range(1,self.d):
            tmp = self.poly_mult(tmp, [1,self.GF.gf_exp[i]])
        return tmp

    def poly_mult(self,poly1,poly2):
        len_poly1, len_poly2 = len(poly1), len(poly2)
        len_poly3 = len(poly1) + len(poly2) - 1
        poly3 = [0]*(len_poly3)
        for ind1 in range(len(poly1)):
            for ind2 in range(len(poly2)):
                poly3[ind1 + ind2] ^= self.GF.mult(poly1[ind1], poly2[ind2])
        return poly3

    def encode(self,m):
        return self.poly_mult(self.gen_poly, m)

    def decode(self,r):
        sydrome = [self.GF.poly_eval(r,self.GF.gf_exp[i]) for i in range(1,self.d)]
        #print(sydrome)

    def gen_all_codewords(self):
        message = [0]*self.k
        codewords = []
        for i in range(self.GF.size ** self.k):
            codewords.append(self.poly_mult(message, self.gen_poly))
            for j in range(self.k-1,-1,-1):
                if message[j] == self.GF.sizem1:
                    message[j] = 0
                else:
                    message[j]+=1
                    break
        return codewords


if __name__ == '__main__':

    RS0  = reed_solomon(360,35,0x409,m=10)
    print("n bits per symbol",RS0.m)
    print("n message symb per frame",RS0.k)
    print("n parity  symb per frame",RS0.n - RS0.k)

    message_frame = np.random.randint(low=0,high=(1<<RS0.m),size=RS0.k)
    encoded_frame = RS0.encode(message_frame)

    for n_errors in [0,10,16,17,18,22]:
        received_frames = deepcopy(encoded_frame)
        uniq_error_locations = np.random.sample(list(range(RS0.n),  n_errors))
        uniq_error_magnitudes= np.random.sample(list(range(1,1024), n_errors))
        
