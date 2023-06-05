import numpy as np
import matplotlib.pyplot as plt
#import sympy as sym


class Sistema1:
    def __init__(
        self,
        mass: float,
        k: float,
        Fo: float,
        w: float,
        c: float,
        xo: float,
        vo: float,
        ti:float,
        tf:float,
        dt:float
    ) -> None:

        # Inicializando variáveis
        self.mass = mass
        self.k = k
        self.Fo = Fo
        self.w = w
        self.c = c
        self.xo = xo
        self.vo = vo
        self.ti = ti
        self.tf = tf
        self.dt = dt

        # Variáveis calculadas
        self.wn = np.sqrt(self.k / self.mass)
        self.r = self.w / self.wn
        self.deltast = self.Fo / self.k
        self.xi = c/(2*np.sqrt(mass*k))
        if self.xi < 1:
            self.wd = np.sqrt(1 - self.xi**2) * self.wn
        else:
            self.wd = 1e-5

    def function(self) -> tuple[np.array, np.array]:
        # Sistema amortecido com força harmônica externa
        self.t = np.arange(self.ti, self.tf, self.dt)

        big_X = self.deltast / (np.sqrt((1 - self.r**2) ** 2 + (2 * self.xi * self.r) ** 2))


        if self.r == 1:
            phi = np.pi / 2
        elif self.r > 1:
            phi = np.arctan(2 * self.xi * self.r / (1 - self.r**2)) + np.pi
        else:
            phi = np.arctan(2 * self.xi * self.r / (1 - self.r**2))

        self.xp = big_X*np.cos(self.w*self.t-phi)
#       xhom = np.exp(-xi*wn*t) * (c1 * np.cos(wd*t) + c2 * np.sin(wd*t))

        self.c1 = self.xo - self.xp[0]
        self.c2 = ( self.vo - big_X*self.w*np.sin(phi) + self.xi*self.wn*self.c1 ) / self.wd
        self.xhom = np.exp(-self.xi*self.wn*self.t) * (self.c1 * np.cos(self.wd*self.t) + self.c2 * np.sin(self.wd*self.t))
        self.x = self.xhom + self.xp


        return (self.t,self.x)



    def plot_vibration_behaviour(self) -> None:
        x_axis, y_axis = self.function()

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 5), tight_layout=True)

        ax.plot(x_axis, y_axis)

        ax.set(title="Comportamento do Sistema", xlabel="Tempo [s]", ylabel="x(t)")

        #fig.savefig("Trabalho/Resposta_total.png")

        plt.show()


    def plot_transient_continuos(self) -> None:
        self.function()

        transient = self.xp
        homogeneous = self.xhom
        total = self.x




        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(10, 10), tight_layout=True)

        ax[0].plot(self.t, total)
        ax[1].plot(self.t, transient)
        ax[2].plot(self.t, homogeneous)

        ax[0].set(
            title="Resposta do sistema",
            xlabel="t",
            ylabel="x",
        )
        ax[1].set(
            title="Resposta Transiente",
            xlabel="t",
            ylabel="x",
        )
        ax[2].set(
             title="Resposta Homogênea",
             xlabel="t",
             ylabel="x",
        )

        #fig.savefig("Trabalho/Resposta_dividida.png")

        plt.show()


# Trabalho
sist = Sistema1(mass=1, k=1000, Fo=-100, w=50, c=5, xo=0.7, vo=30,ti = 0, tf = 5, dt = 0.02)
sist.plot_vibration_behaviour()
#sist.plot_transient_continuos()

# TESTE
# sist.function()
# print('x0:',sist.xo)
# print('x(0):',sist.x[0])
# print('v0:',sist.vo)
# print("x'(0):",(np.diff(sist.x)/sist.dt)[0])
# fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 5), tight_layout=True)
# F = ( sist.Fo * np.cos(sist.w * sist.t) ) [:-2]
# Teste = sist.mass * np.diff( np.diff(sist.x) )/(sist.dt ** 2) + (sist.c * np.diff(sist.x)/sist.dt )[:-1] + ( sist.k*sist.x )[:-2] 
# ax.plot(sist.t[:-2],np.abs(Teste - F))
# plt.show()