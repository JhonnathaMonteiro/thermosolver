
from __future__ import division, print_function, absolute_import

__all__ = ['vdw_base', 'rk_mod', 'pr_mod', 'bwrs']


from math import sqrt, cos, acos, exp, log, pi
from .database import BdD
from thermosolver.config import config_bwrs, config_vdw, config_rk_mod, config_pr_mod
from .cubicsolver import solverAnalitic


class vdw_base(object):
    """docstring for vdw"""
    def __init__(self,T,comp,*args,**kwargs):
        R, Omega_a, Omega_b = config_vdw.conf()
        Tc,Pc,Vc,w = BdD.get_dados(comp,ret=[0,1,2,3])

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
        

    def Tboyle(self):
        from config_vdw import conf
        Omega_a, Omega_b = conf()[1:]
        Tc = BdD.get_dados(comp,ret=[0])[0]
        Tboyle_r = Omega_a/Omega_b
        Tboyle = Tboyle_r*Tc
        return {'Tboyle':Tboyle}

class rk_mod(object):
    """docstring for rk_mod"""
    def __init__(self,T,comp,*args,**kwargs):
        R, Omega_a, Omega_b = config_rk_mod.conf()
        Tc,Pc,Vc,w = BdD.get_dados(comp,ret=[0,1,2,3])

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
        A = [-A_*B_, A_-B_**2-B_, -1]
        
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
    def __init__(self,T,comp,*args,**kwargs):
        R, Omega_a, Omega_b = config_pr_mod.conf()
        Tc,Pc,Vc,w = BdD.get_dados(comp,ret=[0,1,2,3])

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
<<<<<<< HEAD
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
=======
        return pr_mod.solverAnalitic(self,T,P,A,R)

    def baseTV(T,P,comp):
        pass
>>>>>>> master


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

    def __init__(self):
        super(bwrs, self).__init__()

    
    def parametros(T,componente):
        R,A,B = config_bwrs.conf()
        Tc,Pc,Vc,w = BdD.get_dados(componente,ret=[0,1,2,3])

        # cálculo dos parâmetros caracteristicos de cada substância
        A0 = (A[1]+B[1]*w)*R*Tc*Vc #(bar.cm6/mol2)
        B0 = (A[0]+B[0]*w)*Vc #(cm3/mol)
        C0 = (A[2]+B[2]*w)*R*Tc**3*Vc #(bar.cm6.K2/mol2)
        D0 = (A[8]+B[8]*w)*R*Tc**4*Vc #(bar.cm6.K3/mol2)
        E0 = (A[10]+B[10]*w*exp(-3.8*w))*R*Tc**5*Vc #(bar.cm6.K4/mol2)
        a = (A[5]+B[5]*w)*R*Tc*Vc**2 #(bar.cm9/mol3)
        b = (A[4]+B[4]*w)*Vc**2 #(bar.cm9/mol3)
        c = (A[7]+B[7]*w)*R*Tc**3*Vc**2 #(bar.cm9.K2/mol3)
        d = (A[9]+B[9]*w)*R*Tc**2*Vc**2 #(bar.cm9.K/mol3)
        Alpha = (A[6]+B[6]*w)*Vc**3 #(cm9/mol3)
        Gamma = (A[3]+B[3]*w)*Vc**2 #(cm6/mol2)

        # cálculo dos parâmetros K1, K2, K3 e K4
        K1 = B0*R*T-A0-C0/T**2+D0/T**3-E0/T**4
        K2 = b*R*T-a-d/T
        K3 = Alpha*(a+d/T)
        K4 = c/T**2

        return A0,B0,C0,D0,E0,a,b,c,d,Alpha,Gamma,K1,K2,K3,K4,Tc

    # bwrs para T e V independentes
    def TV(self,T,V,componente):
        from config_bwrs import conf
        R,A,B = conf()
        A0,B0,C0,D0,E0,a,b,c,d,Alpha,Gamma,K1,K2,K3,K4,Tc = bwrs.parametros(T,componente)
        
        # Os coeficientes calculados aqui são variáveis locais
        rho = 1/V
        A1 = R*T*rho
        A2 = K1*rho**2 
        A3 = K2*rho**3
        A4 = K3*rho**6
        A5 = K4*rho**3*(1+Gamma*rho**2)*exp(-Gamma*rho**2)
        P = A1 + A2 + A3 + A4 + A5
        Z = P*V/(R*T)
        return (P,Z) 

    # bwrs para T e P independentes
    def TP(self,T,P,componente,fase='vapor',Itemax=100,Tol=1e-10):

        R,A,B = config_bwrs.conf()
        A0,B0,C0,D0,E0,a,b,c,d,Alpha,Gamma,K1,K2,K3,K4,Tc = bwrs.parametros(T,componente)

        # estimativa inicial
        if fase == 'vapor':
            Z = 1
        elif fase == 'liquido':
            Z = 0.003
        else:
            print('As fases devem ser: liquido ou vapor')
            exit(1)

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
            Omega_4 = K4/(R*T)*rho**2*(1+Gamma*rho**2)\
            *exp(-Gamma*rho**2)
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
        
    # cálculo do segundo e terceiro coeficiente do virial
    def Virial_23(self,T,componente):
        
        from config_bwrs import conf
        R,A,B = conf()
        A0,B0,C0,D0,E0,a,b,c,d,Alpha,Gamma,K1,K2,K3,K4,Tc = bwrs.parametros(T,componente)
        
        B = K1/(R*T)
        C = (K2 + K4)/(R*T)
        return {'B':B, 'C':C}
    

    # cálculo da temperatura de boyle
    # em fase de desenvolvimento
    def Tboyle(self,T,componente,fase='vapor',Itemax=100,Tol=1e-10,Tr=3):
        from config_bwrs import conf
        R,A,B = conf()
        A0,B0,C0,D0,E0,a,b,c,d,Alpha,Gamma,K1,K2,K3,K4,Tc = bwrs.parametros(T,componente)
        # cálculo dos coeficientes da função objetivo
        Aux0 = -E0/(B0*R*Tc**5)
        Aux1 = -D0/(B0*R*Tc**4)
        Aux2 = -C0/(B0*R*Tc**3)
        Aux4 = -A0/(B0*R*Tc)
        # Inicialização das variáveis de controle do processo iterativo
        T0 = float('inf')
        Ite = 1
        L = [False, False]
        # método de newton raphson
        while True:
            # função objetivo
            F = Tr**5 + Aux4*Tr**4 + Aux2*Tr**2 \
            + Aux1*Tr + Aux0
            L[0] = Ite>Itemax
            L[1] = abs(F)<Tol or abs(Tr-T0)<Tol
            if L[0] or L[1]: break
            dF = 5*Tr**4 + 4*Aux4*Tr**3 + 2*Aux2*Tr\
            + Aux1
            if abs(dF)<Tol:break
            T0 = Tr
            Tr = Tr - F/dF
            Ite += 1
            
        T_boyle = Tc*Tr

        return {'Tboyle':T_boyle}
