import numpy as np
import pickle as pkl
import sys

class RunsFactory():
    def __init__(self, root, **kwargs):
        '''
        This class helps managing with runs generated from the SMC module.
        The database is constructed upon 8 tables: 
        rootdir, algorithm, N, T, steps, linPoints and logPoints
        '''
        self.tables = ['rootdir', 'algorithm', 'N', 'T', 'steps', \
                       'linPoints', 'logPoints']
        self.root = root
        try:
            self.load(root)
        except FileNotFoundError:
            self.path = kwargs.get('path')
            self.db = []

    def load(self, root):
        self.db = pkl.load(open(root, 'rb'))

    def create(self):
        self.db.append({})
    
    def insert(self, index, params):
        if 
        try:
            for table in self.tables:
                self.db[index][table].append(params[table])
        except IndexError:
            print('Index is invalid...')
            print(f'Creating new run at index = {len(self.db)}')
            self.create()
    def get_index(self, params):
        
    def open(self, idx):
        self.current = self.db