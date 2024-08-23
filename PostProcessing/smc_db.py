import numpy as np
import pandas as pd
import pickle as pkl
import pprint
from colorama import Fore, Style

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
            print(f'Root does not exist. {Fore.GREEN}Create new DB with {Fore.CYAN}self.insert{Style.RESET_ALL}')
            self.db = []

    def load(self, root):
        self.db = pkl.load(open(root, 'rb'))

    def insert(self, rootdir, index=-1):
    
        # reading params values
        params = self.read_params(rootdir)
        if params in self.db:
            print(f'{Fore.RED}IDLING....\nDatabase already contains a run with params:{Style.RESET_ALL}')
            idx = self.db.index(params)
            self.pprint(self.db[idx])
            print(f'at index {Fore.RED}{idx}{Style.RESET_ALL}')

        else:
            if index != -1:
                print(f'{Fore.GREEN}MODYFING DB at index {Fore.RED}{index}{Style.RESET_ALL}')
                print(f'{Fore.YELLOW}Previous parameters were:{Style.RESET_ALL}')
                self.pprint(self.db[index])
                print(f'{Fore.YELLOW}New parameters are:{Style.RESET_ALL}')
                self.pprint(params)
                modify = input("Do you want to continue? (Y/N): ").strip().upper()
                if modify == 'Y':
                    self.db[index] = params
            else:
                print(f'{Fore.GREEN}CREATING new run at index {Fore.RED}{len(self.db)}{Style.RESET_ALL} with parameters:')
                self.pprint(params)
                self.db.append({})
                self.db[-1] = params
            
        # saving
        with open(self.root, 'wb') as db:
            pkl.dump(self.db, db)

    # def get_index(self, params):
    def show(self):
        print(f'{Fore.RED}····································································\
                    {Style.RESET_ALL}')
        for i, run in enumerate(self.db):
            print(f'{Fore.RED}Index {i}{Style.RESET_ALL}')
            self.pprint(run)
            print(f'{Fore.RED}····································································\
                    {Style.RESET_ALL}')
            
    # def set_run(self, idx):
    #     self.run = self.db[idx]
    #     print(f'Database now pointing at run {Fore.RED}{idx}{Style.RESET_ALL} with params:')
    #     pprint.pprint([self.db[idx]])

    @staticmethod
    def pprint(dic):
        key_order = ['rootdir', 'algorithm', 'N', 'T', 'steps', 'linPoints', 'logPoints']
        for key in key_order:
            print(f'{Fore.CYAN}{key}{Style.RESET_ALL}: {dic[key]}')

    @staticmethod
    def read_params(rootdir):
        params = pd.read_csv(f'{rootdir}/params.txt', delim_whitespace=True)
        return params.iloc[0].to_dict()
        # add error if doesn't follow smc module out format