'''
This file contains functions for some of the frequently used operation
when doing efficient simulations for error control codes
'''
def euclidean_distance(a, b):
    assert len(a) == len(b)
    return sum([pow(i-j, 2) for i,j in zip(a,b)])

def int_to_bin_list(x,size = None):
    if size == None: return [int(i) for i in bin(x)[2:]]
    return [(x>>i) & 1 for i in range(size-1,-1,-1)]

def bin_list_to_int(l):
    tmp, mask = 0, 1
    for i in range(len(l)-1,-1,-1):
        if l[i]: tmp |= mask
        mask <<= 1
    return tmp

def bin_mat_to_list(mat):
    return [bin_list_to_int(row) for row in mat]
    
def list_to_bin_mat(l, size = None):
    return [int_to_bin_list(integer, size) for integer in l]

def transpose_2D_list(l):
    r, c = len(l), len(l[0])
    t = [[0]*r for i in range(c)]
    for i in range(r):
        for j in range(c):
            t[j][i] = l[i][j]
    return t

def transpose_int_vector(A, n_cols_A):
    '''Transpose (a vector of integers representation of) binary matrix'''
    n_cols_B = n_rows_A = len(A)
    n_rows_B = n_cols_A
    B = [0] * n_rows_B
    for j, jc in zip(range(n_cols_A), range(n_cols_A - 1, -1, -1)):
        mask_col_A = 1 << jc
        for i, ic in zip(range(n_rows_A), range(n_rows_A - 1, -1, -1)):
            # B[indB][indA] = A[indA][indB] (or) B[j][i] = A[i][j]
            B[j] ^= ( (A[i] & mask_col_A) >> jc) << ic 
    return B

def print_int_vector(A, n_cols_A):
    for i in A:
        print(bin(i)[2:].zfill(n_cols_A))

def eye_mat(n):
    return [[1 if j == i else 0 for j in range(n)] for i in range(n)]

def eye_int_vector(n):#return a identity matrix of dimension n
    return [1<<i for i in range(n-1,-1,-1)]

def sum_of_bits(x):
    ans = 0
    while x > 0:
        ans += x & 1
        x >>= 1
    return ans
    #return sum([i!='0' for i in bin(c)])

def xor_of_bits(x):
    return sum_of_bits(x) & 1

def hamming_distance(c1, c2):
    return sum_of_bits(c1 ^ c2)
    #return sum([i!='0' for i in bin(c1 ^ c2)])

class sum_of_bits_class():
    def __init__(self, nbits = 8):
        '''nbits >= 16 or 24 takes too much pre compute time'''
        self.nbits  = nbits
        self.mem    = [0] * (1 << nbits)
        self.mem[1] = 1
        for j in range(1, nbits):
            l = width = 1 << j
            r = 1 << (j + 1)
            for ind0, ind1 in zip(range(l,r), range(l)):
                self.mem[ind0] = self.mem[ind1] + 1
    
    def __call__(self,x):
        ans  = 0
        mask = (1 << self.nbits) - 1
        while x > 0:
            ans += self.mem[x & mask]
            x >>= self.nbits
        return ans

class linear_block_code(object):
    def __init__(self,G=None,H=None,n=None,k=None):
        self.G = G
        self.H = H
        self.n, self.k, self.p = n, k, n - k 
        self.H_t = transpose_int_vector(self.H, n)
        self.construct_syndrome_table()

    def encode(self, message):
        return self.matmul(self.G, message, self.k - 1) 

    def decode(self, received_msg):
        return self.correct_error(received_msg) >> self.p #discarding parity bits

    def correct_error(self,received_msg):
        syndrome = self.matmul(self.H_t, received_msg, self.n - 1)
        return received_msg ^ self.syndrome_table[syndrome]
    
    def matmul(self, A, x, x_len_minus_1):
        tmp = 0
        for ind, row in enumerate(A):
            if x & (1 << (x_len_minus_1 - ind)):
                tmp ^= row
        return tmp

    def construct_syndrome_table(self):
        self.syndrome_table = [ -1 ] * ( 1 << (self.n - self.k) )
        self.syndrome_table[0], n_syndromes_found = 0, 1
        frontier = [ [0,0] ] # (syndrome, error vector) pairs
        nxt = []
        while n_syndromes_found < len(self.syndrome_table) \
                and len(frontier) > 0: #to account for non full rank H
            for s,e in frontier:
                for ind in range(self.n):
                    mask = 1 << (self.n - 1 - ind)
                    if e & mask: continue
                    s1, e1 = s ^ self.H_t[ind], e ^ mask
                    if self.syndrome_table[s1] == -1:
                        self.syndrome_table[s1] = e1
                        n_syndromes_found += 1
                    nxt.append([s1, e1])
            frontier = []
            nxt = []

def hamming_code_matrices(r):
    '''(n,k,d) = (2**r - 1, 2**r - 1 - r, 3)'''
    n = (1<<r) - 1
    k = n - r
    p = n_minus_k = n - k 
    
    x_power_of_2 = 1 
    def is_power_of_2(x):
        nonlocal x_power_of_2
        if x == x_power_of_2:
            x_power_of_2 <<= 1  #get next power of two
            return 1
        return 0

    #P contains all non zero values which are not powers of 2
    P = [i for i in range(1, n + 1) if not is_power_of_2(i) ]
    P_transpose = transpose_int_vector(P,r)
    #P_mat = list_to_bin_mat(P, size = p)
    #P_mat_transpose = transpose_2D_list(P_mat)
    #P_transpose = bin_mat_to_list(P_mat_transpose)

    Ik = eye_int_vector(k)#identity matrix dimension k
    Ir = eye_int_vector(r)#identity matrix dimension r

    G = [(i << n_minus_k) + p for i,p in zip(Ik,P)]
    H = [(p << n_minus_k) + i for i,p in zip(Ir,P_transpose)]
    return G,H,n,k

