
from __future__ import division, print_function, absolute_import

__all__ = ['vdw_base', 'rk_mod', 'pr_mod', 'bwrs']


from math import sqrt, cos, acos, exp, log, pi
from database import BdD
from thermosolver.config import config_bwrs, config_vdw, config_rk_mod, config_pr_mod
from cubicsolver import solverAnalitic


class vdw_mod(object):
    """docstring for vdw"""
    def __init__(self,T=float,comp=None,*args,**kwargs):
        R, Omega_a, Omega_b = config_vdw.conf()
        _,_,Tc,Pc,Vc,w = BdD.get_dados(comp)

        self.T = T
        self.Tc = Tc
        self.Pc = Pc
        self.Vc = Vc
        self.Tr = T/Tc
        self.R = R
        self.comp = comp
        self.Omega_a = Omega_a
        self.Omega_b = Omega_b

        for k in kwargs.keys():
            if k in ('P', 'T', 'V', 'comp'):
                self.__setattr__(k,kwargs[k])

        if hasattr(self, 'P'):
            self.Pr = self.P/Pc 

    def __repr__(self):
        return str(self.vdw())


    # T e P independentes
    def vdw(self):
        Omega_a = self.Omega_a 
        Omega_b = self.Omega_b 
        
        if hasattr(self, 'P'):
            # Constantes caracteristicas
            Tr = self.Tr 
            Pr  = self.Pr 
            # Calculo das constantes A*=A_ e B*=B_
            A_  = Omega_a*Pr/Tr**2
            B_  = Omega_b*Pr/Tr
            # Coeficientes da equacao cubica
            self.A = [-A_*B_, A_,-(B_ + 1)]
            # Solucao analitica
            return solverAnalitic(self)

        elif hasattr(self, 'V'):
            T = self.T 
            V = self.V
            Tc = self.Tc 
            Pc = self.Pc 
            R = self.R 

            Alpha = Alpha_r = Alpha_c = 1
            # cálculo das constantes características ac, b, a
            ac = Omega_a*R**2*Tc**2/(Alpha_c*Pc)
            b = Omega_b*R*Tc/Pc
            a = ac*Alpha
            # cálculo da pressão em bar
            P = R*T/(V - b) - a/V**2
            # cálculo do fator de compressibilidade
            Z = P*V/(R*T)
            return {'Z':Z,'P':P}

class rk_mod(object):
    """docstring for rk_mod"""
    def __init__(self,T=float,comp=None,*args,**kwargs):
        R, Omega_a, Omega_b = config_rk_mod.conf()
        _,_,Tc,Pc,Vc,w = BdD.get_dados(comp)

        self.T = T
        self.Tc = Tc
        self.Pc = Pc
        self.Vc = Vc
        self.w = w
        self.Tr = T/Tc
        self.R = R
        self.comp = comp
        self.Omega_a = Omega_a
        self.Omega_b = Omega_b

        for k in kwargs.keys():
            if k in ('P', 'T', 'V', 'comp'):
                self.__setattr__(k,kwargs[k])

        if hasattr(self, 'P'):
            self.Pr = self.P/Pc 
    
    def __repr__(self):
        return str(self.RK())

    def baseTP(self):
        
        Tr = self.Tr
        Pr = self.Pr
        Omega_a = self.Omega_a
        Omega_b = self.Omega_b
        Alpha_r = self.Alpha_r
        
        # cálculo dos parâmetros adimensinais
        A_ = Pr*Alpha_r*Omega_a/Tr**2 
        B_ = Pr*Omega_b/Tr

        # Coeficientes da equacao cubica
        self.A = [-A_*B_, A_-B_**2-B_, -1]
        
        # Solucao analitica
        return solverAnalitic(self)

    def baseTV(self,T,V,componente,Alpha,Alpha_c):
        T = self.T
        Tc = self.Tc
        Pc = self.Pc
        V = self.V
        R = self.R
        Alpha = self.Alpha
        Alpha_c = self.Alpha_c
        Omega_a = self.Omega_a
        Omega_b = self.Omega_b
        
        # constantes caracteristicas
        ac = Omega_a*(R*Tc)**2/(Alpha_c*Pc)
        b = Omega_b*R*Tc/Pc
        a = ac*Alpha

        # calculo da pressão em bar
        P = R*T/(V-b)-a/(V*(V+b))

        # fator de compressibilidade
        Z = P*V/(R*T)

        return {'P':P,'Z':Z}

    def RK(self):
        # T e P independentes
        if hasattr(self, 'P'):
            Tr = self.Tr
            # cálculo do parâmetro Alpha
            self.Alpha_r = 1/sqrt(Tr)
            return rk_mod.baseTP(self)
        
        # T e V independentes
        elif hasattr(self, 'V'):
            T = self.T
            Tc = self.Tc
            # cálculo dos parâmetros Alpha
            self.Alpha_c = 1/sqrt(Tc)
            self.Alpha = 1/sqrt(T)
            return rk_mod.baseTV(self)

    def wilson(self):
        w = self.w
        Tr = self.Tr
        # cálculo do parâmetro Alpha
        m = 1.57 + 1.62*w
        self.Alpha_r = Tr*(1+m*(1/Tr-1))
        if hasattr(self, 'P'):
            return rk_mod.baseTP(self)
        elif hasattr(self, 'V'):
            self.Alpha_c = 1
            return rk_mod.baseTV(self)


    def SRK(self):
        w = self.w
        Tr = self.Tr
        # cálculo do parâmetro Alpha
        m = 0.480 + 1.574*w - 0.176*w**2 
        self.Alpha_r = (1+m*(1-sqrt(Tr)))**2
        if hasattr(self, 'P'):
            return rk_mod.baseTP(self)
        elif hasattr(self, 'V'):
            self.Alpha_c = 1
            return rk_mod.baseTV(self)

class pr_mod(object):
    """docstring for pr_mod"""
    def __init__(self,T=float,comp=None,*args,**kwargs):
        R, Omega_a, Omega_b = config_pr_mod.conf()
        _,_,Tc,Pc,Vc,w = BdD.get_dados(comp)

        self.T = T
        self.Tc = Tc
        self.Pc = Pc
        self.Vc = Vc
        self.w = w
        self.Tr = T/Tc
        self.R = R
        self.comp = comp
        self.Omega_a = Omega_a
        self.Omega_b = Omega_b

        for k in kwargs.keys():
            if k in ('P', 'T', 'V', 'comp'):
                self.__setattr__(k,kwargs[k])

        if hasattr(self, 'P'):
            self.Pr = self.P/Pc 
    
    def __repr__(self):
        return str(self.PR())

    def baseTP(self):
        Tr = self.Tr
        Pr = self.Pr
        Omega_a = self.Omega_a
        Omega_b = self.Omega_b
        Alpha_r = self.Alpha_r

        # calculo dos parametros reduzidos
        A_ = Omega_a*Alpha_r*Pr/Tr**2
        B_ = Omega_b*Pr/Tr
        # calculo dos coeficentes da equacao cubica
        self.A = [B_**3+B_**2-A_*B_, -(3*B_**2+2*B_-A_), B_-1]   
        # Solucao analitica

        return solverAnalitic(self)

    def baseTV(self):
        T = self.T
        Tc = self.Tc
        Pc = self.Pc
        V = self.V
        R = self.R
        Alpha = self.Alpha
        Alpha_c = self.Alpha_c
        Omega_a = self.Omega_a
        Omega_b = self.Omega_b
        
        # constantes caracteristicas
        ac = Omega_a*(R*Tc)**2/(Alpha_c*Pc)
        b = Omega_b*R*Tc/Pc
        a = ac*Alpha

        # calculo da pressão em bar
        P = R*T/(V-b)-a/(V**2+2*b*V-b**2)

        # fator de compressibilidade
        Z = P*V/(R*T)

        return {'P':P,'Z':Z}

    def baseTV(T,P,comp):
        pass

    def PR(self):
        Tr = self.Tr
        w = self.w
        m = 0.37464 + 1.54226*w - 0.26992*w**2
        self.Alpha_r = (1 + m*(1-sqrt(Tr)))**2
        self.Alpha_c = 1
        self.Alpha = self.Alpha_r
        if hasattr(self, 'P'):
            return pr_mod.baseTP(self)
        elif hasattr(self, 'V'):
            return pr_mod.baseTV(self)

class bwrs(object):
    """docstring for bwrs"""
    def __init__(self,T=float,comp=None,*args,**kwargs):
        R, A, B = config_bwrs.conf()
        _,_,Tc,Pc,Vc,w = BdD.get_dados(comp)

        self.T = T
        self.Tc = Tc
        self.Pc = Pc
        self.Vc = Vc
        self.w = w
        self.comp = comp
        self.Tr = T/Tc
        self.R = R
        self.A = A
        self.B = B
        # cálculo dos parâmetros caracteristicos de cada substância
        self.A0 = (A[1]+B[1]*w)*R*Tc*Vc #(bar.cm6/mol2)
        self.B0 = (A[0]+B[0]*w)*Vc #(cm3/mol)
        self.C0 = (A[2]+B[2]*w)*R*Tc**3*Vc #(bar.cm6.K2/mol2)
        self.D0 = (A[8]+B[8]*w)*R*Tc**4*Vc #(bar.cm6.K3/mol2)
        self.E0 = (A[10]+B[10]*w*exp(-3.8*w))*R*Tc**5*Vc #(bar.cm6.K4/mol2)
        self.a = (A[5]+B[5]*w)*R*Tc*Vc**2 #(bar.cm9/mol3)
        self.b = (A[4]+B[4]*w)*Vc**2 #(bar.cm9/mol3)
        self.c = (A[7]+B[7]*w)*R*Tc**3*Vc**2 #(bar.cm9.K2/mol3)
        self.d = (A[9]+B[9]*w)*R*Tc**2*Vc**2 #(bar.cm9.K/mol3)
        self.Alpha = (A[6]+B[6]*w)*Vc**3 #(cm9/mol3)
        self.Gamma = (A[3]+B[3]*w)*Vc**2 #(cm6/mol2)
        # cálculo dos parâmetros K1, K2, K3 e K4
        self.K1 = self.B0*R*T-self.A0-self.C0/T**2+self.D0/T**3-self.E0/T**4
        self.K2 = self.b*R*T-self.a-self.d/T
        self.K3 = self.Alpha*(self.a+self.d/T)
        self.K4 = self.c/T**2

        for k in kwargs.keys():
            if k in ('P', 'T', 'V', 'comp'):
                self.__setattr__(k,kwargs[k])

        if hasattr(self, 'P'):
            self.Pr = self.P/Pc 

    def __repr__(self):
        return str(self._bwrs())

    def _bwrs(self):
        if hasattr(self, 'P'):
            return bwrs.TP(self)
        elif hasattr(self, 'V'):
            return bwrs.TV(self)
        
    # bwrs para T e V independentes
    def TV(self):
        T=self.T
        V=self.V
        R=self.R
        Gamma=self.Gamma
        K1=self.K1
        K2=self.K2
        K3=self.K3
        K4=self.K4
        
        rho = 1/V
        A1 = R*T*rho
        A2 = K1*rho**2 
        A3 = K2*rho**3
        A4 = K3*rho**6
        A5 = K4*rho**3*(1+Gamma*rho**2)*exp(-Gamma*rho**2)
        P = A1 + A2 + A3 + A4 + A5
        Z = P*V/(R*T)
        return {'P':P,'Z':Z} 

    # bwrs para T e P independentes
    def TP(self,fase='vapor',Itemax=100,Tol=1e-10):
        T = self.T
        P = self.P
        R = self.R
        Gamma=self.Gamma
        K1=self.K1
        K2=self.K2
        K3=self.K3
        K4=self.K4
        Tc=self.Tc

        # estimativa inicial
        if fase == 'vapor':
            Z = 1
        elif fase == 'liquid':
            Z = 0.003

        # inicialização das ariáveis de controle do processo iterativo: 
        Z0 = float('inf')
        Ite = 1
        L = [False, False, False]

        # método de newton-raphson
        while True:
            # densidade molar através da presão
            rho = P/(Z*R*T)
            # cálculo dos coeficientes Omega_
            Omega_1 = K1/(R*T)*rho
            Omega_2 = K2/(R*T)*rho**2
            Omega_3 = K3/(R*T)*rho**5
            Omega_4 = K4/(R*T)*rho**2*(1+Gamma*rho**2)*exp(-Gamma*rho**2)
            # cálculo da função objetivo
            F = 1+Omega_1+Omega_2+Omega_3+Omega_4-Z
            # convergência
            L[0] = Ite>Itemax
            L[1] = abs(F)<Tol and abs(Z-Z0)<Tol
            if L[0] or L[1]: break
            # cálcular as derivadas
            d_Omega_1 = -Omega_1/Z
            d_Omega_2 = -2*Omega_2/Z
            d_Omega_3 = -5*Omega_3/Z
            d_Omega_4 = K4/(R*T)*(-2*rho**2-4*Gamma*rho**4)/Z*exp(-Gamma*rho**2)\
            +Omega_4*(2*Gamma*rho**2)/Z
            dF = d_Omega_1+d_Omega_2+d_Omega_3+d_Omega_4-1
            L[2] = abs(dF)<1e-20
            if L[2]: break #verificando pontos de mínimo ou máximo
            # cálcular um novo valor para Z via método de Newton-Raphson
            Z0 = Z
            Z = Z-F/dF
            Ite = Ite+1

        # cálculo do volume molar em cm3/mol
        V = Z*R*T/P
        return {'Z':Z, 'V':V}   


if __name__ == '__main__':
    a = bwrs(T=473.15,P=40,comp='carbon dioxide')
    print(a)
