class BdD(object):
    """docstring for BdD"""
    def __init__(self):
        super(BdD, self).__init__()


    def get_dados(componente):
        import csv
        from collections import namedtuple
        # Name;Formula;MW;Tc;Pc;Vc;Rho_c;Zc;w
        # alocando espaço para as listas

        import os
        rel_path = "DadosCriticos.csv"
        folder = os.path.join(os.path.dirname(__file__), 'critical_data')
        file = os.path.join(folder, 'Yaws Collection.csv')

        Especie = namedtuple('Componente','CASRN Name Tc Pc Vc w')
        saida = None
        with open(file, 'r') as csvfile:
            csvreader = csv.reader(csvfile,delimiter='\t',quoting=csv.QUOTE_NONNUMERIC)
            headers = next(csvreader)
            a = []
            for row in csvreader:
                if componente in (row[0],row[1]):
                    saida = Especie(*row)
                    print
                    break  

            if not saida:
                raise ValueError('Especie não encontrada no banco de dados')

        return saida

if __name__ == '__main__':
    _,_,Tc,Pc,Vc,w  = BdD.get_dados('hexamethyldisilazane')
    print(Tc)