
from __future__ import division, print_function, absolute_import

__all__ = ['vdw', 'rk_mod', 'pr_mod', 'bwrs']


from math import sqrt, cos, acos, exp, log, pi
from .database import BdD
from thermosolver.config import config_bwrs, config_vdw, config_rk_mod, config_pr_mod


class cubic(object):
    """docstring for cubic"""
    
    def solverAnalitic(self):
        T = self.T
        P = self.P
        R = self.R
        A = self.A
        # Solucao analitica
        Q_2 = (A[0]+2*A[2]**3/27-A[2]*A[1]/3)/2
        P_3 = (3*A[1]-A[2]**2)/9
        Delta = P_3**3 + Q_2**2

        if Delta >= 0:
            u = (-Q_2+sqrt(Delta))**(1/3)
            v = -P_3/u
            x1 = -A[2]/3 + u + v

            # Fator de compressibilidade
            Z = x1
            V = Z*R*T/P
            return {'Z':Z,'V':V}

        else:
            Phi = acos(-Q_2/sqrt(-P_3**3))
            x1 = -A[2]/3 + 2*sqrt(-P_3)*cos(Phi/3)
            x2 = -A[2]/3 + 2*sqrt(-P_3)*cos((Phi-pi)/3)
            x3 = -A[2]/3 + 2*sqrt(-P_3)*cos((Phi+pi)/3)
            
            # Fator de compressibilidade
            Zv, Zl = max(x1,x2,x3), min(x1,x2,x3)
            Vv = Zv*R*T/P
            Vl = Zl*R*T/P
            return {'Zv':Zv,'Zl':Zl,'Vv':Vv,'Vl':Vl}

class vdw(cubic):
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
    def TP(self):
        if not hasattr(self, 'P'):
            raise AttributeError('Invalid method, inform T= , P= , comp= " " ')
        
        # Constantes caracteristicas
        Tr = self.Tr 
        Pr  = self.Pr 
        Omega_a = self.Omega_a 
        Omega_b = self.Omega_b 
                        
        # Calculo das constantes A*=A_ e B*=B_
        A_  = Omega_a*Pr/Tr**2
        B_  = Omega_b*Pr/Tr
        
        # Coeficientes da equacao cubica
        self.A = [-A_*B_, A_,-(B_ + 1)]


        # Solucao analitica
        return vdw.solverAnalitic(self)

    # T e V independentes
    def TV(self):
        T = self.T 
        V = self.V
        Tc = self.Tc 
        Pc = self.Pc 
        R = self.R 
        Omega_a = self.Omega_a 
        Omega_b = self.Omega_b 
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

class rk_mod(cubic):
    """docstring for rk_mod"""
    def __init__(self):
        super(rk_mod, self).__init__()
    
    def baseTP(self,T,P,componente,Alpha_r):
        R,Omega_a,Omega_b = config_rk_mod.conf()        
        Tc,Pc,w = BdD.get_dados(componente,ret=[0,1,3])
        # propriedadede reduzida
        Tr, Pr = T/Tc, P/Pc
        
        # cálculo dos parâmetros adimensinais
        A_ = Pr*Alpha_r*Omega_a/Tr**2 
        B_ = Pr*Omega_b/Tr

        # Coeficientes da equacao cubica
        A = [-A_*B_, A_-B_**2-B_, -1]
        
        # Solucao analitica
        return rk_mod.solverAnalitic(self,T,P,A,R)

    def baseTV(self,T,V,componente,Alpha,Alpha_c):
        R, Omega_a,Omega_b = config_rk_mod.conf()
        Tc,Pc,w = BdD.get_dados(componente,ret=[0,1,3])
        
        # propriedadede reduzida
        Tr = T/Tc
        
        # constantes caracteristicas
        ac = Omega_a*(R*Tc)**2/(Alpha_c*Pc)
        b = Omega_b*R*Tc/Pc
        a = ac*Alpha

        # calculo da pressão em bar
        P = R*T/(V-b)-a/(V*(V+b))

        # fator de compressibilidade
        Z = P*V/(R*T)

        return {'P':P,'Z':Z}

    def RK_TP(self,T,P,componente):
        # propriededades físicas da substância
        Tc = BdD.get_dados(componente,ret=[0])[0]

        # constantes características e R:
        R,Omega_a,Omega_b = config_rk_mod.conf()

        # propriedadede reduzida
        Tr = T/Tc

        # cálculo do parâmetro Alpha
        Alpha_r = 1/sqrt(Tr)

        return rk_mod.baseTP(self,T,P,componente,Alpha_r)

    def RK_TV(self,T,V,componente):
        # propriededades físicas da substância
        Tc = BdD.get_dados(componente,ret=[0])[0]

        # constantes características e R:
        R,Omega_a,Omega_b = config_rk_mod.conf()

        # propriedadede reduzida
        Tr = T/Tc

        # cálculo dos parâmetros Alpha
        Alpha_c = 1/sqrt(Tc)
        Alpha = 1/sqrt(T)

        return rk_mod.baseTV(self, T,V,componente,Alpha,Alpha_c)

    def wilson_TP(self,T,P,componente):
        # propriededades físicas da substância
        Tc,w= BdD.get_dados(componente,ret=[0,3])

        # constantes características e R:
        R,Omega_a,Omega_b = config_rk_mod.conf()

        # propriedadede reduzida
        Tr = T/Tc

        # cálculo do parâmetro Alpha
        m = 1.57 + 1.62*w
        Alpha_r = Tr*(1+m*(1/Tr-1))

        return rk_mod.baseTP(T,P,componente,Alpha_r)

    def wilson_TV(self,T,V,componente):
        # propriededades físicas da substância
        Tc,w = BdD.get_dados(componente,ret=[0,3])

        # constantes características e R:
        R,Omega_a,Omega_b = config_rk_mod.conf()

        # propriedadede reduzida
        Tr = T/Tc

        # cálculo dos parâmetros Alpha
        m = 1.57 + 1.62*w
        Alpha = Tr*(1+m*(1/Tr-1))
        Alpha_c = 1
        return rk_mod.baseTV(T,V,componente,Alpha,Alpha_c)

    def SRK_TP(self,T,P,componente):
        # propriededades físicas da substância
        Tc,w= BdD.get_dados(componente,ret=[0,3])

        # constantes características e R:
        R, Omega_a,Omega_b = config_rk_mod.conf()

        # propriedadede reduzida
        Tr = T/Tc

        # cálculo do parâmetro Alpha
        m = 0.480 + 1.574*w - 0.176*w**2 
        Alpha_r = (1+m*(1-sqrt(Tr)))**2

        return rk_mod.baseTP(T,P,componente,Alpha_r)

    def SRK_TV(self,T,V,componente):
        Tc,w = BdD.get_dados(componente,ret=[0,3])

        # constantes características e R:
        R, Omega_a,Omega_b = config_rk_mod.conf()

        # propriedadede reduzida
        Tr = T/Tc

        # cálculo dos parâmetros Alpha
        m = 0.480 + 1.574*w - 0.176*w**2 
        Alpha = (1+m*(1-sqrt(Tr)))**2
        Alpha_c = 1
        return rk_mod.baseTV(T,V,componente,Alpha,Alpha_c)

class pr_mod(cubic):
    """docstring for pr_mod"""
    def __init__(self):
        super(pr_mod, self).__init__()
    
    def baseTP(sef,T,P,componente,Alpha_r):
        R, Omega_a, Omega_b = config_pr_mod.conf()
        Tc,Pc = BdD.get_dados(componente,ret=[0,1])
        # propriedades reduzidas
        Tr = T/Tc
        Pr = P/Pc
        # calculo dos parametros reduzidos
        A_ = Omega_a*Alpha_r*Pr/Tr**2
        B_ = Omega_b*Pr/Tr
        # calculo dos coeficentes da equacao cubica
        A = [B_**3+B_**2-A_*B_, -(3*B_**2+2*B_-A_), B_-1]        
        # Solucao analitica
        return pr_mod.solverAnalitic(self,T,P,A,R)

    def baseTV(T,P,comp):
        pass


    def PR_TP(self,T,P,componente):
        Tc,w = BdD.get_dados(componente,ret=[0,3])
        Tr = T/Tc
        m = 0.37464 + 1.54226*w - 0.26992*w**2
        Alpha = (1 + m*(1-sqrt(Tr)))**2
        return pr_mod.baseTP(T,P,componente,Alpha)

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
    
    # cálculo da temperatura de Boyle
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
