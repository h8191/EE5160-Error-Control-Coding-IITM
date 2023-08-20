"""
This file contains a python object which can do pretty much everything 
related to Non-recursive Convolutional Codes(some mods may be necessary)

encoding an input, trellis decoding(both soft ml and hard ml), 
state mapping from current stage to both next stage and previous stage. 

Author: Harshavardhan P
"""
from pprint import pprint
from functions import xor_of_bits, sum_of_bits, hamming_distance, \
bin_list_to_int, int_to_bin_list, euclidean_distance

class convolutional_codes(object):
    def __init__(self,G = [[0b1],[0b1]], memory = 2):
        self.G              = G             #Generator matrix each integer corresponds to an output
        self.memory         = memory        #number of D-Flip flops (bits per state)
        self.num_states     = 1 << memory
        self.num_outputs    = len(G)        #number of output ports in encoder (D-Flip flop)

        self.gen_adj_list()

        self.int_to_BPSK = [[-1 if (output >> i) & 1 else 1
                                for i in range(self.num_outputs-1,-1,-1)]
                                for output in range(1 << self.num_outputs)]

    def gen_adj_list(self):
        '''generates adjaceny lists to move to next state and previous state'''

        #mapping of states from current state to next state
        adj_next = [[None for j in range(2)] for i in range(self.num_states)]
        #mapping of states from current state to previous state
        adj_prev = [[None for j in range(2)] for i in range(self.num_states)]
        for state in range(self.num_states):
            for input_bit in (0,1):
                output, nxt_state = self.encoder(input_bit, state)
                adj_next[state][input_bit]      = [nxt_state, output]
                popped_bit = state & 1
                adj_prev[nxt_state][popped_bit]  = [state,     output]
        self.adj_next = adj_next
        self.adj_prev = adj_prev           

    def encoder(self, input_bit = 0, state = 0):
        '''given state(integer)(concatenation of bits) and input(bit) computes
        the outputs, (using Generator matrix) and next state'''
        #by default encoder starts in all zero state
        state = (input_bit << self.memory) | state  #append input bit to the left side of state
        nxt_state = state >> 1                      #pop the right most (least recent) bit of the state
         
        output = [xor_of_bits(g & state) for g in self.G]   #find the value of each output of encoder
        output = bin_list_to_int(output)               #stack the output bits into integer for faster computation
        return output, nxt_state

    def compute_section_of_trellis(self, r_val, trellis, current_stage_index,
                                    zero_termination = False):

        num_stages  = len(trellis)
        min_CM      = 0

        #precompute branch metrics for received vector and each possible output
        branch_metrics_list = [self.branch_metric(output, r_val) 
                                for output in range(1 << self.num_outputs)]
        
        next_stage_index = (current_stage_index + 1) % num_stages

        #reset parents for all states in next stage
        for next_state in range(self.num_states):
            trellis[next_stage_index][next_state][1] = -1
        
        for curr_state in range(self.num_states):

            curr_state_CM, activated = trellis[current_stage_index][curr_state]
            if activated == -1:
                continue #node(current state) does n't have a parent


            for next_state, output in self.adj_next[curr_state]:
                next_state_CM, parent = trellis[next_stage_index][next_state]
                
                tmp_CM = curr_state_CM + branch_metrics_list[output]
                #store which of the previous states gave the least cumulative metric
                if parent == -1 or tmp_CM < next_state_CM or zero_termination:
                    trellis[next_stage_index][next_state] = [tmp_CM, curr_state]            
                    min_CM                                = min(min_CM, tmp_CM)

                #inputs other than 0 are not used in the current section of trellis
                if zero_termination:
                    break
        print('curr next ind:', current_stage_index, next_stage_index)
        pprint(trellis[current_stage_index])
        pprint(trellis[next_stage_index])

        '''To prevent floating point overflow subtract the minimum CM from all CM'''
        #for state in range(self.num_states):
        #    trellis[ind + 1][state] -= min_CM

    def backtrack_trellis(self, delay, stage_index, trellis, len_output = 1):
        curr_state = 0
        backtrack_output = []
        num_stages = len(trellis)

        for i in range(delay):
            parent = trellis[stage_index][curr_state][1]

            if len_output != 1 or i == delay - 1:
                for input_bit in range(2):
                    state, output = self.adj_next[parent][input_bit] 
                    if state == curr_state:  break
                backtrack_output.append(input_bit)

            stage_index = (stage_index - 1) % num_stages
            curr_state  = parent

        return backtrack_output[0] if len_output == 1 else backtrack_output[::-1]

    def finite_delay_decoding(self, r_vec, soft = 0):
        self.branch_metric = self.BPSK_AWGN_branch_metric if soft else hamming_distance

        delay      = 5 * self.memory
        len_r_vec  = len(r_vec)
        num_stages = min(delay, len_r_vec) + 1

        trellis = [[ [-1, -1] for i in range(self.num_states) ]
                              for i in range(num_stages)]

        trellis[0][0][0] =  0
        trellis[0][0][1] =  0 #activate all zero state in first stage by setting its parent to non default value
        decoder_output   = [0] * (0 if delay >= len_r_vec else len_r_vec - delay)

        print('len r vec:', len_r_vec)
        current_stage_index = 0
        for r_index in range(min(delay, len_r_vec)):
            self.compute_section_of_trellis(r_vec[r_index], trellis, current_stage_index)
            current_stage_index += 1
        #by the end of this for loop, current_stage_index = delay = num_stages - 1 = len(trellis) - 1

        for r_index in range(delay, len_r_vec):
            zero_termination = len_r_vec - r_index <= self.memory

            decoder_output[r_index - delay] = \
            self.backtrack_trellis(delay, current_stage_index, trellis, len_output = 1)
             
            self.compute_section_of_trellis(r_vec[r_index], trellis, current_stage_index, zero_termination)
            current_stage_index = (current_stage_index + 1) % num_stages

        decoder_output.extend(
        self.backtrack_trellis(num_stages - 1, current_stage_index, trellis, len_output = -1))

        self.trellis = trellis
        return decoder_output
        

    def decode(self, r, soft=True, verbose = 0):
        branch_metric = self.BPSK_AWGN_branch_metric if soft else hamming_distance
        
        frontier = [0] #all zero state
        num_stages = len(r) + 1
        trellis  = [[[-1,-1]   for i in range(self.num_states)]
                                for i in range(num_stages)] #CM, prev state
        
        trellis[0][0][0] = 0 #setting Cumulative Metric of all zero state in 0th stage to 0
        all_states_activated = False
        
        branch_metrics_all_stages = []

        #computing cumulative metrics of nodes in trellis
        for ind, r_val in enumerate(r):
            if not all_states_activated: 
                iterator, frontier = frontier, set()#iterator will be frontier and frontier will be empty set
            else:
                iterator = range(self.num_states)

            branch_metrics_list = [branch_metric(output, r_val) 
                        for output in range(1 << self.num_outputs)]
            branch_metrics_all_stages.append(branch_metrics_list)
            
            for curr_state in iterator:
                curr_state_CM = trellis[ind][curr_state][0]

                for next_state, output in self.adj_next[curr_state]:
                    next_state_CM, parent = trellis[ind + 1][next_state]

                    tmp = curr_state_CM + branch_metrics_list[output]
                    #store which node in the previous stage gave it the current minimum cumulative metric
                    if parent is None or parent < 0 or tmp < next_state_CM:
                        trellis[ind + 1][next_state] = [tmp, curr_state]

                    if not all_states_activated:
                        frontier.add(next_state)

            all_states_activated |= len(frontier) == self.num_states
    
        #backtracing survivor path
        curr_state = 0 #all zero state
        path      = [[] for i in range(len(r))]
        for stage_index in range(num_stages - 1, 0, -1):
            parent = trellis[stage_index][curr_state][1]
            input_bit  = (curr_state >> (self.memory - 1)) & 1
            output     = self.adj_next[parent][input_bit]
            
            curr_state = parent
            path[stage_index - 1] = [input_bit, output]
        
        def myroundoff1(l,x=2):
            return [None if i is None else round(i,x) for i in l]

        def myroundoff2(l,x=2):
            return [myroundoff1(i,x) for i in l]

        if verbose:
            print('\n received vector')
            print(r)
            print('\ncumulative metrics')
            for stage in trellis:
                print(myroundoff2(stage))
            print('\nbranch metrics')
            for branch_metrics_list in branch_metrics_all_stages:
                print(myroundoff1(branch_metrics_list))
        
            print('\nadj next: (state, input)   -> (next state, output)\n', self.adj_next)
            print('\nadj prev: (state, pop bit) -> (prev state, output)\n', self.adj_prev)
                
            print('\nsurvivor path: input & output')
            print(path)

        transmitted_message_estimate  = [i[0] for i in path]
        transmitted_codeword_estimate = [i[1] for i in path]
        self.trellis  = trellis
        return transmitted_message_estimate, transmitted_codeword_estimate

    def encode(self, message, zero_termination = True):
        state = 0 #initial state of encoder
        output = []

        for input_bit in message:
            tmp, state = self.encoder(input_bit, state)
            output.append(tmp)

        if zero_termination:
            for i in range(self.memory):
                tmp, state = self.encoder(0, state)
                output.append(tmp)

        return output

    def BPSK_AWGN_branch_metric(self,bits,r):
        return euclidean_distance(self.int_to_BPSK[bits],r)

if __name__ == '__main__':
    G = [0b101, 0b110]
    code0 = convolutional_codes(G,memory = 2)
    #print(code0.adj_next)
    
    G = [0b111, 0b110]
    code1 = convolutional_codes(G, memory = 2)
    
    import numpy as np
    np.random.seed(42)
    num_bits = 20
    bits_tx = np.random.randint(low = 0, high = 2, size = num_bits)
    print(code1.encode(bits_tx))
    tmp     = np.array([j for i in code1.encode(bits_tx)
                          for j in int_to_bin_list(i, code1.num_outputs)    ])
    symb_tx = (1 - 2 * tmp).reshape(-1,2)
    symb_rx = symb_tx + 0.8 * np.random.randn(* symb_tx.shape) 
    bits_rx1 = code1.decode(symb_rx, 1, verbose = 0)[0]
    print('infinite')
    pprint(code1.trellis)
    bits_rx2 = code1.finite_delay_decoding(symb_rx, 1)
    print('finite')
    pprint(code1.trellis)
    print(bits_tx)
    print(np.array(bits_rx1).astype('uint8')) 
    print(np.array(bits_rx2).astype('uint8'))
