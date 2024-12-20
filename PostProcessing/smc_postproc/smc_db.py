import numpy as np
import pandas as pd
import pickle as pkl
from colorama import Fore, Style
import os 

class RunsFactory():
    def __init__(self, root):
        '''
        This class helps managing runs generated from the SMC module.
        The database is constructed upon 9 tables: 
        rootdir, N, T, tau, tw, cycles, linPoints, logPoints and p_swap
        '''
        self.root = root
        try:
            self.load(root)
        except FileNotFoundError:
            print(f'Root does not exist. {Fore.GREEN}Create new DB with ' + \
                  f'{Fore.CYAN}self.insert{Style.RESET_ALL}')
            self.db = []

    def load(self, root):
        self.db = pkl.load(open(root, 'rb'))
        self.show()

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
                print(f'{Fore.GREEN}CREATING new run at index ' + 
                      f'{Fore.RED}{len(self.db)}{Fore.GREEN} with parameters:{Style.RESET_ALL}')
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

    # def get_observable(self):
    #     return None
    
    # def get_cfg(self):
    #     return None
    
    @staticmethod
    def pprint(dic):
        key_order = ['rootdir', 'N', 'T', 'tau', 'tw', 'cycles', 'logPoints', \
                     'linPoints', 'p_swap']
        for key in key_order:
            print(f'{Fore.CYAN}{key}{Style.RESET_ALL}: {dic[key]}')

    @staticmethod
    def read_params(rootdir):
        if os.path.exists(rootdir):
            try:
                params = pd.read_csv(f'{rootdir}/params.txt', delim_whitespace=True)
                return params.iloc[0].to_dict()
            except FileNotFoundError:
                raise FileNotFoundError(f"There is no {Fore.CYAN}params.txt{Style.RESET_ALL} " + 
                                        f"file in '{rootdir}'")
        else:
            raise FileNotFoundError(f"{Fore.YELLOW}{rootdir}{Style.RESET_ALL} does not exist.")