import numpy as np
import matplotlib.pyplot as plt
#import sympy as sym

import pandas as pd
# from sklearn.metrics import mean_absolute_error, mean_absolute_percentage_error

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
        ti: float,
        tf: float,
        dt: float
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
        self.t = np.arange(self.ti, self.tf+self.dt, self.dt)

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
sist_1 = Sistema1(mass=1, k=1000, Fo=-100, w=50, c=5, xo=0.7, vo=30,ti = 0, tf = 5, dt = 0.02)
sist_2 = Sistema1(mass=1, k=1000, Fo=-100, w=50, c=5, xo=0.7, vo=30,ti = 0, tf = 5, dt = 0.002)
# sist.plot_vibration_behaviour()
# sist.plot_transient_continuos()

sist_1.function()
t_anal_1 = sist_1.t
pos_anal_1 = sist_1.x
# print(x_anal.shape)
sist_2.function()
t_anal_2 = sist_2.t
pos_anal_2 = sist_2.x



x_anal_1 = np.vstack((pos_anal_1,t_anal_1)).T
x_anal_2 = np.vstack((pos_anal_2,t_anal_2)).T
# print(x_anal_1.shape)
# print(x_anal_2.shape)

x_conv_1 = np.loadtxt("conv_dt_02.csv",delimiter=",")
x_conv_2 = np.loadtxt("conv_dt_002.csv",delimiter=",")

x_dif_1 = np.loadtxt("dif_dt_02.csv",delimiter=",")
x_dif_2 = np.loadtxt("dif_dt_002.csv",delimiter=",")

# print(x_conv_1.shape)
# print(x_conv_2.shape)

# Plotando individualmente
plt.figure(figsize=(10,6))
#analitico

# plt.plot(x_anal_1[:,1],x_anal_1[:,0], "b.-", markersize = 3.5, label = "Resposta do sistema - Método Analitico")
# plt.plot(x_anal_2[:,1],x_anal_2[:,0], "b.-", markersize = 3.5, label = "Resposta do sistema - Método Analitico")

# plt.plot(x_conv_1[:,1],x_conv_1[:,0], "b.-", markersize = 3.5, label = "Resposta do sistema - Integral de Convolução")
# plt.plot(x_conv_2[:,1],x_conv_2[:,0], "b.-", markersize = 3.5, label = "Resposta do sistema - Integral de Convolução")

# plt.plot(x_dif_1[:,1],x_dif_1[:,0], "b.-", markersize = 3.5, label = "Resposta do sistema - Diferenças Finitas")
# plt.plot(x_dif_2[:,1],x_dif_2[:,0], "b.-", markersize = 3.5, label = "Resposta do sistema - Diferenças Finitas")


#configurando plot
# plt.grid()
# plt.legend()
# plt.xlim(0,5)
# plt.ylim(-1,1.2)
# plt.ylabel("Posição [m]")
# plt.xlabel("Tempo [s]")

# plt.savefig("x_t_analitico_dt_02.png", dpi = 1000)
# plt.savefig("x_t_analitico_dt_002.png", dpi = 1000)

# plt.savefig("x_t_conv_dt_02.png", dpi = 1000)
# plt.savefig("x_t_conv_dt_002.png", dpi = 1000)

# plt.savefig("x_t_dif_dt_02.png", dpi = 1000)
# plt.savefig("x_t_dif_dt_002.png", dpi = 1000)

plt.close()

# plotando tudo junto
# plt.plot(x_anal_1[:,1],x_anal_1[:,0], "b", markersize = 3.5, label = "Resposta do sistema - Método Analitico")
# plt.plot(x_conv_1[:,1],x_conv_1[:,0], "r", markersize = 3.5, label = "Resposta do sistema - Integral de Convolução")
# plt.plot(x_dif_1[:,1],x_dif_1[:,0]  , "c", markersize = 3.5, label = "Resposta do sistema - Diferenças Finitas")

plt.plot(x_anal_2[:,1],x_anal_2[:,0], "b", markersize = 3.5, label = "Resposta do sistema - Método Analitico")
plt.plot(x_conv_2[:,1],x_conv_2[:,0], "r", markersize = 3.5, label = "Resposta do sistema - Integral de Convolução")
plt.plot(x_dif_2[:,1],x_dif_2[:,0]  , "c", markersize = 3.5, label = "Resposta do sistema - Diferenças Finitas")

plt.grid()
plt.legend()
plt.xlim(0,5)
plt.ylim(-1,1.2)
plt.ylabel("Posição [m]")
plt.xlabel("Tempo [s]")
# plt.savefig("x_t_compare_dt_02.png", dpi = 500)
plt.savefig("x_t_compare_dt_002.png", dpi = 500)

# MAE = mean_absolute_error(x1,x_ref)
# MAPE = 1e2*mean_absolute_percentage_error(x1,x_ref)
# print(f'Mean absolute error: {1e3*MAE:.3} mm')
# print(f'Mean absolute percentage error: {1e2*MAPE:.3}%')
# textcov = f'Metodo Integral de Convolução: \n MAPE: {MAPE:.2f}%'
# Dif1 = x_ref - x1
# Dif2 = x_ref - x3

# MAE = mean_absolute_error(x3,x_ref)
# MAPE = 1e2*mean_absolute_percentage_error(x3,x_ref)
# print(f'Mean absolute error: {1e3*MAE:.3} mm')
# print(f'Mean absolute percentage error: {1e2*MAPE:.3}%')
# textdif = f'Metodo Diferenças Finitas: \n MAPE: {MAPE:.2f}%'

# # Plotar o gráfico no Python
# import matplotlib.pyplot as plt

# fig, axis = plt.subplots(ncols=1, nrows=2, sharex=True, figsize = (6,4))
# fig.subplots_adjust(hspace=0)

# axis[0].plot(t_ref, x_ref,'ro-', label='Analítico', markersize = 2, markerfacecolor='None')
# axis[0].plot(t1, x1, 'bx--', label='Metodo Integral de Convolução', markersize = 3, dashes=(5, 10))
# axis[0].plot(t3, x3,'ko--', label='Diferenças Finitas', markersize = 2, markerfacecolor='None', dashes=(5, 10))
# axis[1].plot(t1, Dif1,'bx--', label='Metodo Integral de Convolução', markersize = 3, dashes=(5, 10))
# axis[1].plot(t_ref, Dif2,'ko--', label='Diferenças Finitas', markersize = 2, markerfacecolor='None', dashes=(5, 10))

# props = dict(facecolor='none', edgecolor='none')
# axis[0].text(0.50, 0.225, textcov, transform=axis[0].transAxes, verticalalignment='top', bbox=props)
# axis[0].text(0.50, 0.975, textdif, transform=axis[0].transAxes, verticalalignment='top', bbox=props)


# axis[1].set_xlabel('Tempo t [s]')
# axis[0].set_ylabel('Deslocamento x(t) [m]')
# axis[1].set_ylabel('Erro [m]')
# # axis[0].set_ylim(-1.15,1.15)
# axis[0].set_xlim(0,5)


# # plt.title('Sobreposição de gráficos')
# plt.legend()
# plt.show()

# fig.savefig('superprop.png', dpi=600, bbox_inches='tight', pad_inches=0)
