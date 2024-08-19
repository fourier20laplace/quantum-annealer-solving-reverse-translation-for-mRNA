import numpy as np
import dimod
from dwave.system import LeapHybridSampler
import os
import time


class QA:
    def __init__(self, linear_array, quadratic_mat, evaluate, gen_res, print_verbose, verbose=False):
        self.linear_array = linear_array
        self.quadratic_mat = quadratic_mat
        self.verbose = verbose
        self.evaluate = evaluate
        self.gen_res = gen_res
        self.print_verbose = print_verbose

    def qa(self):
        start_time = time.time()
        os.environ['DWAVE_API_TOKEN'] = 'DEV-b11514f3df392bddec88c07fe3f47b571ed3079a'
        bqm = dimod.BinaryQuadraticModel(self.linear_array, np.triu(self.quadratic_mat, k=1), 0, dimod.BINARY)
        # Create a sampler instance
        sampler = LeapHybridSampler()

        # Sample the BQM
        sampleRes = sampler.sample(bqm)
        sample = sampleRes.first.sample

        # Convert the sample to a list or array
        vector_result = [sample[v] for v in bqm.variables]
        energy = sampleRes.first.energy

        # (val,)这一项的val形如array([-95.77537761])需要再取一下
        if round(energy, 2) != round(self.evaluate(vector_result)[0][0], 2):
            raise ValueError("Energy mismatch error")
        end_time = time.time()
        res_dict = self.gen_res(vector_result)
        res_dict["time"] = end_time - start_time
        if self.verbose:
            self.print_verbose(res_dict)
        return res_dict
