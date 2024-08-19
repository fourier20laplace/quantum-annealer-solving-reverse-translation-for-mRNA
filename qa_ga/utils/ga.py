from deap import base, creator, tools, algorithms
import random
import numpy as np
import time


class GA:
    def __init__(self, codon_len, aa_len, codon_count_per_aa_list, evaluate, gen_res, print_verbose, verbose=False):
        self.codon_len = codon_len
        self.aa_len = aa_len
        self.codon_count_per_aa_list = codon_count_per_aa_list
        self.evaluate = evaluate
        self.gen_res = gen_res
        self.print_verbose = print_verbose
        self.verbose = verbose

    # 初始化时生成自定义的个体
    def my_custom_individual(self):
        individual = np.zeros(self.codon_len)
        for i in range(self.aa_len):
            start = sum(self.codon_count_per_aa_list[:i])
            end = sum(self.codon_count_per_aa_list[:i + 1])
            rand_idx = random.randint(start, end - 1)
            individual[rand_idx] = 1
        return creator.Individual(individual.tolist())

    def custom_mutation(self, individual, indpb):
        if random.random() < indpb:
            idx_aa = random.randint(0, self.aa_len - 1)
            start = sum(self.codon_count_per_aa_list[:idx_aa])
            end = start + self.codon_count_per_aa_list[idx_aa]
            idx_mut = random.randint(start, end - 1)
            for i in range(start, end):
                individual[i] = 0
            individual[idx_mut] = 1

        return individual,

    def custom_crossover(self, ind1, ind2):
        left = random.randint(0, self.aa_len - 1)
        right = random.randint(0, self.aa_len - 1)

        if left > right:
            left, right = right, left

        start = sum(self.codon_count_per_aa_list[:left])
        end = sum(self.codon_count_per_aa_list[:right + 1])
        ind1[start:end], ind2[start:end] = ind2[start:end], ind1[start:end]
        return ind1, ind2

    def ga(self):

        # Check if the classes already exist and delete them if necessary
        if hasattr(creator, 'FitnessMax'):
            del creator.FitnessMax
        if hasattr(creator, 'Individual'):
            del creator.Individual
        # 创建一个最大化适应度的目标
        creator.create("FitnessMax", base.Fitness, weights=(-1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)

        # 定义个体和种群的生成方式
        toolbox = base.Toolbox()
        # toolbox.register("attr_bool", random.randint, 0, 1)  # 二进制编码
        # toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=self.codon_len)
        toolbox.register("individual", self.my_custom_individual)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        # 注册遗传算法的操作
        toolbox.register("evaluate", self.evaluate)
        # toolbox.register("mate", tools.cxTwoPoint)  # 交叉操作
        # toolbox.register("mutate", tools.mutFlipBit, indpb=0.2)  # 变异操作
        # 内置的交叉和变异不好 用自定义的
        toolbox.register("mate", self.custom_crossover)  # Custom crossover
        toolbox.register("mutate", self.custom_mutation, indpb=0.2)  # Custom mutation
        toolbox.register("select", tools.selTournament, tournsize=3)  # 选择操作

        population = toolbox.population(n=100)  # 种群大小
        NGEN = 50  # 代数
        CXPB = 0.5  # 交叉概率
        MUTPB = 0.2  # 变异概率

        # 运行遗传算法
        results = algorithms.eaSimple(population, toolbox, cxpb=CXPB, mutpb=MUTPB, ngen=NGEN, verbose=False)
        best_individual = tools.selBest(population, 1)[0]

        res_dict = self.gen_res(best_individual)

        if self.verbose:
            self.print_verbose(res_dict)
        return res_dict

    def ga_n_times(self, n):
        start_time = time.time()
        ga_res = 100
        res_dict = {}
        for i in range(n):
            ga_cur = self.ga()
            if ga_cur["energy"] < ga_res:
                ga_res = ga_cur["energy"]
                res_dict = ga_cur
        end_time = time.time()
        res_dict["time"] = end_time - start_time
        return res_dict
