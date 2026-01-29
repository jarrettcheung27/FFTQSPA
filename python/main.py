import numpy as np
import fftqspa
from Inner_Code_DNA_Channel_Simulation import DNAChannel


def main():
    # code parameter file name
    parity_filename = "2080_320_16ary.dat"
    mapping_filename = "SignalSet_BPSK-4.txt"
    max_iteration = 50
    k_2 = 320  # inner code length of data bits
    # DNA channel parameters
    Pe = 0.1  # base error rate
    sequencingDepth = 10  # sequencing depth
    innerRedundancy = 114  # total redundancy for inner code

    codec = fftqspa.BCJRQSPA(parity_filename, max_iteration, mapping_filename)

    n_info = codec.info_bits_len()
    n_code = codec.code_bits_len()

    # 生成320条长度为n_0的随机信息比特
    rng = np.random.default_rng(seed=42)
    info_bits = rng.integers(0, 2, size=(n_info,k_2), dtype=np.uint8)

    # 分别对每一条信息比特进行编码
    code_bits = np.empty((n_code, k_2), dtype=np.uint8)
    for i in range(k_2):
        code_bits[:, i] = codec.encoder4bibo(info_bits[:, i])
    #========================Inner code + DNA channel========================
    """
    将inner code和DNA存储信道视为一个复合信道
    输入：k_2条LDPC码字的比特流
    输出：k_2条LDPC码字的比特流的LLR
    """

    voting_scores = DNAChannel(code_bits, Pe, sequencingDepth, innerRedundancy)  # shape: (n_code, k_2)

    #==================================================================
    
    # voting score to llr， 避免出现inf，将0和1分别映射为1e-3和1-1e-3
    eps = 1e-3
    voting_scores = np.clip(voting_scores, eps, 1 - eps)
    llr = np.log((1.0 - voting_scores) / voting_scores)
    rr_bits_prob = (1.0 / (1.0 + np.exp(llr))).astype(np.float64)

    # 分别对每一条码字进行译码
    decoded_bits = np.empty((n_code, k_2), dtype=np.uint8)
    iters = np.empty(k_2, dtype=int)
    for i in range(k_2):
        decoded_bits[:, i], iters[i] = codec.decode4bibo(rr_bits_prob[:, i])

    # 分别系统化编码：信息位在码字末尾（前面是校验位）
    sys_start = n_code - n_info
    bit_errors = np.count_nonzero(decoded_bits[sys_start:] != info_bits)

    print(f"n_info={n_info}, n_code={n_code}")
    print(f"iterations={iters}, bit_errors={bit_errors}")

if __name__ == "__main__":
    main()
