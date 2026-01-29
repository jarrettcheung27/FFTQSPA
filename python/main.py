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
    sequencingDepth = 5  # sequencing depth
    innerRedundancy = 114  # total redundancy for inner code

    codec = fftqspa.BCJRQSPA(parity_filename, max_iteration, mapping_filename)

    n_info = codec.info_bits_len()
    n_code = codec.code_bits_len()

    # 生成320条长度为n_0的随机信息比特
    rng = np.random.default_rng(seed=42)
    info_bits = rng.integers(0, 2, size=(n_info,k_2), dtype=np.uint8)

    # 分别对每一条信息比特进行编码
    print("Q-ary LDPC encoding...")
    code_bits = np.empty((n_code, k_2), dtype=np.uint8)
    for i in range(k_2):
        code_bits[:, i] = codec.encoder4bibo(info_bits[:, i])
    #========================Inner code + DNA channel========================
    """
    将inner code和DNA存储信道视为一个复合信道
    输入：k_2条LDPC码字的比特流
    输出：k_2条LDPC码字的比特流的LLR
    """

    # 简单的BPSK+AWGN信道，输出P(b=0)
    '''
    snr_db = 5.0
    rate = n_info / n_code
    snr_lin = 10 ** (snr_db / 10.0)
    sigma = np.sqrt(1.0 / (2.0 * snr_lin * rate))

    # BPSK: 0->-1, 1->+1
    rr = np.where(code_bits == 0, -1.0, 1.0).astype(np.float64)
    noise = rng.normal(0.0, 1.0, size=code_bits.shape)
    y = rr + sigma * noise

    # 计算P(b=0) (BPSK: 0->-1, 1->+1)
    # LLR = 2*y/sigma^2, p0 = 1/(1+exp(LLR))
    llr = 2.0 * y / (sigma ** 2)
    '''
    #==================================================================

    # 使用inner code + DNA存储信道复合信道
    # 输出voting scores 为多数投票得分, 也是就在一个簇中1的占比（约等于P(b=1)），loss sequence的得分为0.5。
    voting_scores = DNAChannel(code_bits, Pe, sequencingDepth, innerRedundancy)  # shape: (n_code, k_2)

    #==================================================================
    # 计算P(b=0)
    rr_bits_prob =  1 - voting_scores  # P(b=0)

    # 分别对每一条码字进行译码
    print("Q-ary LDPC decoding...")
    decoded_bits = np.empty((n_code, k_2), dtype=np.uint8)
    iters = np.empty(k_2, dtype=int)
    for i in range(k_2):
        decoded_bits[:, i], iters[i] = codec.decode4bibo(rr_bits_prob[:, i])

    # 分别系统化编码：信息位在码字末尾（前面是校验位）
    sys_start = n_code - n_info
    decoded_bits = decoded_bits[sys_start:,:]  # 信息位
    
    # 计算误比特率(BER)和帧错误率(FER), 以decoded_bits[0, :]为1帧,及info_bits[0, :]为原始信息
    total_bit_errors = 0
    total_frame_errors = 0
    for i in range(n_info):
        bit_errors = np.sum(decoded_bits[i, :] != info_bits[i, :])
        total_bit_errors += bit_errors
        if bit_errors > 0:
            total_frame_errors += 1
    BER = total_bit_errors / (n_info * k_2)
    FER = total_frame_errors / n_info
    print(f"BER: {BER}, FER: {FER}")

if __name__ == "__main__":
    main()
