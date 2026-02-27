import numpy as np
import fftqspa
from Inner_Code_DNA_Channel_Simulation import DNAChannel

def main():
    # code parameter file name
    parity_filename = "2080_320_16ary.dat"
    # parity_filename = "2080_320_16ary2.dat"
    mapping_filename = "SignalSet_BPSK-4.txt"
    max_iteration = 10
    k_2 = 320  # inner code length of data bits
    # DNA channel parameters
    PEs = np.linspace(0.05, 0.12, 8)  # different base error rates
    sequencingDepths = [15, 20]  # sequencing depth list
    innerRedundancy = 114  # total redundancy for inner code

    codec = fftqspa.BCJRQSPA(parity_filename, max_iteration, mapping_filename)

    n_info = codec.info_bits_len()
    n_code = codec.code_bits_len()

    for sequencingDepth in sequencingDepths:
        import os
        if not os.path.exists('results'):
            os.makedirs('results')
        results_filename = f'results/FFTQSPA_DNA_Channel_16-ary_SequencingDepth{sequencingDepth}_InnerRedundancy{innerRedundancy}_test.csv'
        with open(results_filename, 'w') as f:
            f.write('Pe,FER\n')

        for Pe in PEs:
            print(f"Simulating for base error rate Pe={Pe}...")
            # 生成320条长度为n_0的随机信息比特
            rng = np.random.default_rng(seed=42)
            info_bits = rng.integers(0, 2, size=(n_info,k_2), dtype=np.uint8)

            # 分别对每一条信息比特进行编码
            print("Q-ary LDPC encoding...")
            code_bits = np.empty((n_code, k_2), dtype=np.uint8)
            for i in range(k_2):
                code_bits[:, i] = codec.encoder4bibo(info_bits[:, i])
            # 使用inner code + DNA存储信道复合信道
            # 输出voting scores 为多数投票得分, 也是就在一个簇中1的占比（约等于P(b=1)），loss sequence的得分为0.5。
            voting_scores = DNAChannel(code_bits, Pe, sequencingDepth, innerRedundancy)  # shape: (n_code, k_2)

            # 计算P(b=0)
            rr_bits_prob =  1 - voting_scores  # P(b=0)

            # 为rr_bits_prob加一个微小值，防止出现0或1的概率，导致LLR无穷大
            epsilon = 1e-5
            rr_bits_prob = np.clip(rr_bits_prob, epsilon, 1 - epsilon)


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
            with open(results_filename, 'a') as f:
                f.write(f'{Pe},{FER}\n')
if __name__ == "__main__":
    main()
