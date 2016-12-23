class BdD(object):
    """docstring for BdD"""
    def __init__(self):
        super(BdD, self).__init__()


    def get_dados(componente,ret=None):
        import csv
        # Name;Formula;MW;Tc;Pc;Vc;Rho_c;Zc;w
        # alocando espaço para as listas

        import os
        rel_path = "DadosCriticos.csv"
        folder = os.path.join(os.path.dirname(__file__), 'critical_data')
        file = os.path.join(folder, 'Yaws Collection.csv')

        with open(file, 'r') as csvfile:
            csvreader = csv.reader(csvfile,delimiter='\t',quoting=csv.QUOTE_NONNUMERIC)

            headers = next(csvreader)

            # Atribuindo os valores as listas
            CASRN = []
            Name = []
            Dados = []
            
            for row in csvreader:
                CASRN.append(row[0])
                Name.append(row[1])
                Dados.append(row[2:])

        try:
            I = Name.index(componente)
        except ValueError: 
            try:
                I = CASRN.index(componente)
            except ValueError:
                raise ValueError('Componente fora do banco de dados')      

        saida = []
        for j in ret:
            saida.append(Dados[I][j])

        if all(isinstance(n,float) for n in saida):
            return saida
        else:
            raise ValueError('Dados incompletos no banco de dados')

if __name__ == '__main__':
    a = BdD.get_dados('benzene',ret=[1,2,3])
    print(a[2])