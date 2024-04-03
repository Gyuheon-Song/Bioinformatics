'''
생명공학과 2021163102 송규헌
https://github.com/Gyuheon-Song/Bioinformatics
2024-1 생명정보학 Problem Set 01
Problem no. 3
Viterbi Algorithm을 구현하여 amino acid sequence "DELIFFLIF" 에 대한 most likely state sequence 구하기
files needed : soluble_sequences_aa_freq.txt , transmembrane_sequences_aa_freq.txt , state_sequences_tr_freq.txt
usage : github에서 repo clone 후에 /Viterbi_algorithm.py  Ctrl+F5로 Run
result : console창에 viterbi matrix와 hidden state sequence출력
'''

'''
임의의 서열을 입력받는 경우에는 line19를 삭제하고,
긴 서열의 빠른 입출력을 위해 아래 코드로 진행

input = sys.stdin.readline
target_seq = input().rstrip()

'''

import sys
from collections import deque 
target_seq = "DELIFFLIF"   # 본 과제에서 decoding할 sequence
seq_matrix = [[['S', 0] for _ in range(len(target_seq))], [['T', 0] for _ in range(len(target_seq))]]
dp = [[['', 0] for _ in range(len(target_seq))] for __ in range(2)]
hidden_state_sequence = deque()

# 전역
sol_map    = {}      # soluble sequence의 frequency dictionary 
trans_map  = {}      # transmembrane sequence의 frequency dictionary
init_prob  = {}      # transition probabilities from the beginning state of the HMM(Starting probability)
trans_prob = {}      # transition probabilities

'''
Solution Sequence의 frequency에 대한 hash map 생성함수
key     : a.a
value   : frequency rate
'''
def Make_Solution_Sequence_Map():
    file = open('soluble_sequences_aa_freq.txt', 'r')
    lst = file.readlines()

    for str in lst :
        tmp = str.split()
        sol_map[tmp[2]] = float(tmp[4])

    file.close()
    return

'''
Transmembrane Sequence의 frequency에 대한 hash map 생성함수
key     : a.a
value   : frequency rate
'''
def Make_Transmembrane_Sequence_Map():
    file = open('transmembrane_sequences_aa_freq.txt', 'r')
    lst = file.readlines()

    for str in lst :
        tmp = str.split()
        trans_map[tmp[2]] = float(tmp[4])

    file.close()
    return

'''
Initializing Probability
시작 확률 설정
'''
def Init_Start_Prob() :
    file = open('state_sequences_tr_freq.txt', 'r')
    lst = file.readlines()
    for str in lst[:2] :
        tmp = str.split()
        init_prob[tmp[3]] = float(tmp[5])
    
    file.close()
    return

'''
Initializing Transition Probability
form 변화에 대한 확률 S->S, S->T, T->S, T->T
이중 dictionary 구조 {OuterKey1 : {InnerKey1 : v, InnerKey2 : vv, ...}, OuterKey2 : {InnerKey_A : vvv, ...}}
'''
def Init_Trans_Prob() :
    file = open('state_sequences_tr_freq.txt', 'r')
    lst = file.readlines()
    for str in lst[3:] :
        tmp = str.split()
        if tmp[3] not in trans_prob :
            trans_prob[tmp[3]] = {tmp[5] : float(tmp[7])}
        else :
            trans_prob[tmp[3]][tmp[5]] = float(tmp[7])

    file.close()
    return

'''
Viterbi Algorithm
'''
def Viterbi(sequence) :

    l = len(sequence)
    # initializing viterbi matrix
    for i in range(l) :
        seq_matrix[0][i][1] = sol_map[sequence[i]]
        seq_matrix[1][i][1] = trans_map[sequence[i]]

    # initial condition
    start_prob = [init_prob['S'], init_prob['T']]
    for j in range(2) :
        seq_matrix[j][0][1] *= start_prob[j]

    # viterbi algorithm using dp(3 dimensional array)
    for aaidx in range(1, l) :
        for j in range(2) :
            for k in range(2) :
                if aaidx == 1 :
                    t_prob = seq_matrix[k][aaidx-1][1] * seq_matrix[j][aaidx][1]
                    if t_prob * trans_prob[seq_matrix[k][aaidx-1][0]][seq_matrix[j][aaidx][0]] > dp[j][aaidx][1] :
                        dp[j][aaidx][0] = seq_matrix[k][aaidx-1][0]
                        dp[j][aaidx][1] = t_prob * trans_prob[seq_matrix[k][aaidx-1][0]][seq_matrix[j][aaidx][0]]
                    continue
                t_prob = dp[k][aaidx-1][1] * seq_matrix[j][aaidx][1]
                if t_prob * trans_prob[seq_matrix[k][aaidx-1][0]][seq_matrix[j][aaidx][0]] > dp[j][aaidx][1] :
                    dp[j][aaidx][0] = seq_matrix[k][aaidx-1][0]   # dp배열의 현재 원소를 [이전 state, 현재 아미노산까지의 1차 HMM 최대 누적곱]
                    dp[j][aaidx][1] = t_prob * trans_prob[seq_matrix[k][aaidx-1][0]][seq_matrix[j][aaidx][0]]

    return

'''
Backtracking으로 hidden state sequence 찾기
'''
def Find_Hidden_State() :
    next_state_idx = -1
    # 가장 마지막 아미노산의 hidden state로 시작
    if dp[0][len(target_seq)-1][1] > dp[1][len(target_seq)-1][1] :
        hidden_state_sequence.appendleft('S')
        next_state_idx = 0
    else :
        hidden_state_sequence.appendleft('T')
        next_state_idx = 1
    
    # 가장 확률이 높은 말단 state부터 역으로 추적
    for i in range(len(target_seq) - 1, 0, -1) :
        hidden_state_sequence.appendleft(dp[next_state_idx][i][0])
        if dp[next_state_idx][i][0] == 'S' :
            next_state_idx = 0
        else :
            next_state_idx = 1
    
    return


if __name__ == "__main__" :

    Make_Solution_Sequence_Map()
    Make_Transmembrane_Sequence_Map()
    Init_Start_Prob()
    Init_Trans_Prob()
    Viterbi(target_seq)
    Find_Hidden_State()

    hidden_state_sequence = ''.join(hidden_state_sequence)
    print("")
    print(f"Hidden state sequence of {target_seq} is", str(hidden_state_sequence))

    # 누적곱과 각 이전노드(경로) 보기 -> 위의 코드에서 배열명 'dp'에 저장
    '''
    for item in dp :
        print(dp)
    '''

