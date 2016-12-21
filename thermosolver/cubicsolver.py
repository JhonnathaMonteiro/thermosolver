from math import acos, sqrt, cos

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
