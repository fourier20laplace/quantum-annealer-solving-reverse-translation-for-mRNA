import os.path

import python_codon_tables as pct
import numpy as np
import sys

sys.path.append(os.path.dirname(__file__))
from ga import GA
from qa import QA


class ReverseTranslation:
    def __init__(self, tar_seq, expression="e_coli_316407", constant_f=0.1, constant_gc=1, constant_r=0.1,
                 constant_epsilon=5, rou_t=0.5, big_num=50, verbose=False):
        self.aa_len = None
        self.codon_count_per_aa_list = None
        self.codon_list = None
        self.codon_len = None
        self.tar_seq = tar_seq
        self.expression = expression
        self.constant_f = constant_f
        self.constant_gc = constant_gc
        self.constant_r = constant_r
        self.constant_epsilon = constant_epsilon
        self.rou_t = rou_t
        self.big_num = big_num
        self.linear_array, self.quadratic_mat = self.get_linear_quadratic_array()
        # 指示是否打印详细信息
        self.verbose = verbose

    @staticmethod
    def max_consecutive_repeats(s):
        """
        计算字符串s中最大连续重复字符数
        :return: 最大连续重复字符数
        """
        max_count = 1  # 最大连续重复字符数，初始值设为1
        current_count = 1  # 当前字符的连续计数
        # 从第二个字符开始遍历
        for i in range(1, len(s)):
            if s[i] == s[i - 1]:
                current_count += 1
                max_count = max(max_count, current_count)
            else:
                current_count = 1

        return max_count

    def get_linear_quadratic_array(self):
        """
        获取线性项和二次项
        :return: 线性项和二次项
        """
        aa_len = len(self.tar_seq)
        self.aa_len = aa_len
        # 选择表达系统
        table = pct.get_codons_table(self.expression)

        # 获取标准codon列表 GC含量字典
        standard_codon_list = []
        gc_content_dict = {}
        for aa in table:
            for codon in table[aa]:
                standard_codon_list.append(codon)
                gc_content_dict[codon] = (codon.count('G') + codon.count('C')) / 3

        # codon_list 或者说量子的列表
        codon_list = []
        codon_count_per_aa_list = []
        # 密码子偏好列表 freq_list
        # 小量防止欠定义
        epsilon = 1e-6
        freq_list = []
        for aa in self.tar_seq:
            for codon in table[aa]:
                codon_list.append(codon)
                val = self.constant_f * np.log(1 / (table[aa][codon] + epsilon))
                freq_list.append(val)
            codon_count_per_aa_list.append(len(table[aa]))
        freq_arr = np.array(freq_list)
        codon_len = len(codon_list)
        self.codon_len = codon_len
        self.codon_list = codon_list
        self.codon_count_per_aa_list = codon_count_per_aa_list
        # GC含量
        s_list = []
        for codon in codon_list:
            s_list.append(gc_content_dict[codon])
        s_array = np.array(s_list) * (2 * self.rou_t * self.constant_gc / aa_len)
        s_array_pow = np.power(s_array, 2) * self.constant_gc / (aa_len ** 2)
        s_mat = np.outer(s_array, s_array)

        linear_array = freq_arr - s_array + s_array_pow - self.constant_epsilon

        # 损失项：避免连续重复核苷酸
        # 预先计算出矩阵

        pre_overlap_matrix = np.zeros((64, 64), dtype=int)
        for i in range(64):
            for j in range(64):
                val = self.max_consecutive_repeats(standard_codon_list[i] + standard_codon_list[j])
                pre_overlap_matrix[i][j] = val ** 2 - 1

        cur_overlap_matrix = np.zeros((codon_len, codon_len), dtype=int)
        start_i = 0
        start_j = codon_count_per_aa_list[0]
        for i in range(aa_len - 1):
            if i != 0:
                start_i += codon_count_per_aa_list[i - 1]
                start_j += codon_count_per_aa_list[i]
            for i_inter in range(start_i, start_i + codon_count_per_aa_list[i]):
                for j_inter in range(start_j, start_j + codon_count_per_aa_list[i + 1]):
                    codon1 = codon_list[i_inter]
                    codon2 = codon_list[j_inter]
                    idx1 = standard_codon_list.index(codon1)
                    idx2 = standard_codon_list.index(codon2)
                    cur_overlap_matrix[i_inter][j_inter] = pre_overlap_matrix[idx1][idx2]

        # 转化为对称阵
        cur_overlap_matrix = np.triu(cur_overlap_matrix) + np.tril(cur_overlap_matrix.T, -1)

        # 保证独热性 tao_mat
        tao_mat = np.zeros((codon_len, codon_len), dtype=int)
        start_i = 0
        for i in range(aa_len):
            if i != 0:
                start_i += codon_count_per_aa_list[i - 1]
            for i_inter in range(start_i, start_i + codon_count_per_aa_list[i]):
                for j_inter in range(start_i, start_i + codon_count_per_aa_list[i]):
                    if i_inter != j_inter:
                        tao_mat[i_inter][j_inter] = self.big_num

        quadratic_mat = 2 * self.constant_gc / (aa_len ** 2) * s_mat + self.constant_r * cur_overlap_matrix + tao_mat
        return linear_array, quadratic_mat

    def evaluate(self, individual):
        ind_arr = np.array(individual)
        val = np.dot(ind_arr, self.linear_array)
        upper_triangle = np.triu(self.quadratic_mat, k=1)
        val += ind_arr @ upper_triangle @ ind_arr.reshape(-1, 1)
        return (val,)

    def gen_res(self, qubits):
        vector_result = qubits
        start = 0
        my_dict = {}
        origin_dict = {}
        res_str = ""
        for i, num in enumerate(self.codon_count_per_aa_list):
            tmpList = []
            originList = []
            strTmp = ""
            for index in range(start, start + num):
                if vector_result[index] == 1:
                    tmpList.append(self.codon_list[index])
                    if not strTmp:
                        strTmp = self.codon_list[index]
                    else:
                        strTmp = strTmp + "/" + self.codon_list[index]
                originList.append(self.codon_list[index])
            res_str = res_str + str(i) + "_" + strTmp + "_"
            origin_dict[i] = (self.tar_seq[i], originList)
            my_dict[i] = (self.tar_seq[i], tmpList)
            start += num
        res_dict = {"qubits": qubits, "energy": self.evaluate(qubits)[0], "codon_seq": res_str,
                    "aa_seq": '_'.join([f"{i}_{char}" for i, char in enumerate(self.tar_seq)]),
                    "origin_dict": origin_dict, "my_dict": my_dict}
        return res_dict

    @staticmethod
    def print_verbose(dict):
        print("二进制优化结果")
        print(dict["qubits"])
        print("每个位置的备选密码子")
        print('\n'.join(f"{key}: {value}" for key, value in dict["origin_dict"].items()))
        print("每个位置的选中密码子")
        print('\n'.join(f"{key}: {value}" for key, value in dict["my_dict"].items()))
        print("能量")
        print(dict["energy"])
        print("目标氨基酸序列")
        print(dict["aa_seq"])
        print("最终的密码子序列")
        print(dict["codon_seq"])
        print("end")

    def compare(self, ga_n_times=2):
        print("start")
        print(f"Target sequence: {self.tar_seq}")
        ga = GA(codon_len=self.codon_len, aa_len=self.aa_len, codon_count_per_aa_list=self.codon_count_per_aa_list,
                evaluate=self.evaluate, gen_res=self.gen_res, print_verbose=self.print_verbose, verbose=self.verbose)
        qa = QA(self.linear_array, self.quadratic_mat, self.evaluate, self.gen_res, self.print_verbose,
                verbose=self.verbose)
        ga_dict = ga.ga_n_times(ga_n_times)
        qa_dict = qa.qa()
        qa_energy = qa_dict["energy"]
        ga_energy = ga_dict["energy"]
        if qa_energy < ga_energy:
            print(f"energy_contrast::qa:{qa_energy} < ga:{ga_energy}")
        else:
            print(f"energy_contrast::qa:{qa_energy} >= ga:{ga_energy}")

        qa_time = qa_dict["time"]
        ga_time = ga_dict["time"]
        if qa_time < ga_time:
            print(f"time_contrast::qa:{qa_time}s < ga:{ga_time}s")
        else:
            print(f"time_contrast::qa:{qa_time}s >= ga:{ga_time}s")

        print(f"qa:{qa_dict['codon_seq']}")
        print(f"ga:{ga_dict['codon_seq']}")

        parts = qa_dict['codon_seq'].strip('_').split('_')
        parts1 = [(int(parts[i]), parts[i + 1]) for i in range(0, len(parts), 2)]
        parts = ga_dict['codon_seq'].strip('_').split('_')
        parts2 = [(int(parts[i]), parts[i + 1]) for i in range(0, len(parts), 2)]

        differences = []
        for (index1, value1), (index2, value2) in zip(parts1, parts2):
            if value1 != value2:
                differences.append({
                    'index': index1,
                    'aa': self.tar_seq[index1],
                    'value1': value1,
                    'value2': value2
                })

        if differences:
            print("Differences found:")
            for diff in differences:
                print(f"Index {diff['index']}: {diff['aa']} - qa_codon: {diff['value1']}, ga_codon: {diff['value2']}")
        else:
            print("identical.")

        print("end")
        print()
        return ga_energy, qa_energy
