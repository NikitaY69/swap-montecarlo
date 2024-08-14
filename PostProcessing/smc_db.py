import numpy as np
import pandas as pd
import pickle as pkl
import pprint

class RunsFactory():
    def __init__(self, root):
        '''
        This class helps managing runs generated from the SMC module.
        The database is constructed upon 7 tables: 
        rootdir, algorithm, N, T, steps, linPoints and logPoints
        '''
        self.root = root
        try:
            self.load(root)
        except FileNotFoundError:
            print(f'Root does not exist. Create new DB with self.insert')
            self.db = []

    def load(self, root):
        self.db = pkl.load(open(root, 'rb'))

    def insert(self, rootdir, index=-1):

        # reading params values
        params = self.read_params(rootdir)

        if index != -1:
            print(f'Modifying DB at index {index}')
            self.db[index] = params
        else:
            print(f'Creating new run at index = {len(self.db)}')
            self.db.append({})
            self.db[-1] = params
            
        # saving
        with open(self.root, 'wb') as db:
            pkl.dump(self.db, db)

    # def get_index(self, params):
    def show(self):
        pprint.pprint(self.db)
        
    def set_run(self, idx):
        self.run = self.db[idx]
        print(f'Database now pointing at run {idx} with params:')
        pprint.pprint([self.db[idx]])

    @staticmethod
    def read_params(rootdir):
        params = pd.read_csv(f'{rootdir}/params.txt', delim_whitespace=True)
        return params.iloc[0].to_dict()
        # add error if doesn't follow smc module out format