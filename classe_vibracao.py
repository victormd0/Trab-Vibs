import numpy as np
import matplotlib.pyplot as plt
import sympy as sym


class Sistema:
    def __init__(
        self,
        mass: float,
        k: float,
        Fo: float,
        w: float,
        psi: float,
        Y: float,
        xo: float,
        vo: float,
    ) -> None:
        # Botar uns assert
        assert mass > 0, "Massa de valor negativo não existe!"
        assert k > 0, "Rigidez de valor negativo não existe!"
        assert w >= 0, "Não existe frequência de valor negativo!"

        # Inicializando variáveis
        self.mass = mass
        self.k = k
        self.Fo = Fo
        self.w = w
        self.psi = psi
        self.Y = Y
        self.xo = xo
        self.vo = vo

        # Variáveis calculadas
        self.wn = np.sqrt(self.k / self.mass)
        self.r = self.w / self.wn
        self.deltast = self.Fo / self.k
        if self.psi < 1:
            self.wd = np.sqrt(1 - self.psi**2) * self.wn
        else:
            self.wd = 1e-5

    def function(self) -> tuple[np.array, np.array]:
        t = np.arange(0, self.wn * 5 / np.pi, 0.01)
        if self.psi == 0:  # Sistema não amortecido
            if self.Fo == 0:  # Sistema não amortecido sem força harmônica externa
                return (
                    t,
                    (
                        self.xo * np.cos(self.wn * t)
                        + (self.vo * np.sin(self.wn * t)) / (self.wn)
                    ),
                )
            else:  # Sistema não amortecido com força externa
                if (
                    self.w == self.wn
                ):  # Sistema não amortecido excitado por uma força harmônica com ressonância
                    return (
                        t,
                        (
                            self.xo * np.cos(self.wn * t)
                            + (self.vo * np.sin(self.wn * t)) / (self.wn)
                            + (self.deltast * self.wn * t * np.sin(self.wn * t)) / 2
                        ),
                    )
                else:  # Sistema não amortecido excitado por uma força harmônica sem ressonância
                    return (
                        t,
                        (
                            (self.xo - self.Fo / (self.k - self.mass * self.w**2))
                            * np.cos(self.wn * t)
                            + (self.vo * np.sin(self.wn * t)) / (self.wn)
                            + (
                                self.Fo
                                * np.cos(self.wn * t)
                                / (self.k - self.mass * self.w**2)
                            )
                        ),
                    )
        else:  # Sistema amortecido
            if self.Fo == 0:  # Sistema amortecido sem força harmônica externa
                if (
                    self.Y != 0
                ):  # Sistema amortecido com base vibrando harmonicamente (Somente resposta permanente)
                    big_X = self.Y * np.sqrt(
                        (1 + (2 * self.psi * self.r) ** 2)
                        / ((1 - self.r**2) ** 2 + (2 * self.psi * self.r) ** 2)
                    )
                    phi = np.arctan(
                        (2 * self.psi * self.r**3)
                        / (1 + (4 * self.psi**2 - 1) * self.r**2)
                    )
                    return (
                        t,
                        big_X * np.sin(self.w * t - phi),
                    )
                else:
                    if (
                        self.psi < 1
                    ):  # Sistema sub amortecido sem força harmônica externa
                        return (
                            t,
                            (
                                self.xo * np.cos(self.wd * t)
                                + (self.vo + self.psi * self.wn * self.xo)
                                / (self.wd)
                                * np.sin(self.wd * t)
                            )
                            * np.exp(-self.psi * self.wn * t),
                        )
                    elif (
                        self.psi == 1
                    ):  # Sistema criticamente amortecido sem força harmônica externa
                        return (
                            t,
                            (self.xo + (self.vo + self.wn * self.xo) * t)
                            * np.exp(-self.wn * t),
                        )
                    else:  # Sistema superamortecido sem força harmônica externa
                        c_1 = (
                            self.xo * self.wn * (self.psi + np.sqrt(self.psi**2 - 1))
                            + self.vo
                        ) / (2 * self.wn * np.sqrt(self.psi**2 - 1))
                        c_2 = (
                            -1
                            * self.xo
                            * self.wn
                            * (self.psi - np.sqrt(self.psi**2 - 1))
                            - self.vo
                        ) / (2 * self.wn * np.sqrt(self.psi**2 - 1))
                        elevated_pos = np.exp(
                            (-1 * self.psi + np.sqrt(self.psi**2 - 1)) * self.wn * t
                        )
                        elevated_neg = np.exp(
                            (-1 * self.psi - np.sqrt(self.psi**2 - 1)) * self.wn * t
                        )
                        return (t, (c_1 * elevated_pos + c_2 * elevated_neg))
            else:  # Sistema amortecido com força harmônica externa
                Xo = np.sqrt(
                    self.xo**2 + (self.psi * self.wn * self.xo / self.wd) ** 2
                )
                if self.xo == 0:
                    phi_o = np.pi / 2
                else:
                    phi_o = np.arctan(
                        -1
                        * (self.vo + self.psi * self.wn * self.xo)
                        / (self.wd * self.xo)
                    )
                envelope = np.exp(-1 * self.psi * self.wn * t)
                big_X = self.deltast / (
                    np.sqrt((1 - self.r**2) ** 2 + (2 * self.psi * self.r) ** 2)
                )
                if self.r == 1:
                    phi = np.pi / 2
                else:
                    phi = np.arctan(2 * self.psi * self.r / (1 - self.r**2))
                return (
                    t,
                    (
                        Xo * envelope * np.cos(self.wd * t - phi_o)
                        + big_X * np.cos(self.w * t - phi)
                    ),
                )

    def plot_vibration_behaviour(self) -> None:
        # TODO: Adicionar informações tipo Amplitude máxima, angulo phi em forma de texto em algum lugar no plot
        x_axis, y_axis = self.function()

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 5), tight_layout=True)

        ax.plot(x_axis, y_axis)

        ax.set(title="Comportamento do Sistema", xlabel="Tempo [s]", ylabel="x(t)")

        fig.savefig("Output/Vibration_Behaviour.png")

        plt.clf()

        return

    def plot_transmissibility_of_movement(self) -> None:
        assert self.Y > 0, "Precisa ter um valor de Y"
        r = np.arange(0, 10, 0.01)
        Td = np.sqrt(
            (1 + (2 * self.psi * r) ** 2)
            / ((1 - r**2) ** 2 + (2 * self.psi * r) ** 2)
        )
        phi = np.arctan(
            (2 * self.psi * r**3) / (1 + (4 * self.psi**2 - 1) * r**2)
        )
        phi_degree = phi * 180 / np.pi

        rmax = round(np.sqrt(np.sqrt(1 + 8 * self.psi**2) - 1) / (2 * self.psi), 3)
        Tdmax = round(
            np.sqrt(
                (1 + (2 * self.psi * rmax) ** 2)
                / ((1 - rmax**2) ** 2 + (2 * self.psi * rmax) ** 2)
            ),
            3,
        )

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(20, 20), tight_layout=True)

        ax[0].plot(r, Td)
        ax[1].plot(r, phi)
        ax2 = ax[1].twinx()
        ax2.plot(r, phi_degree)

        ax[0].text(
            3,
            Tdmax,
            "O valor máximo ocorre em r={}, com Td={}".format(rmax, Tdmax),
            {"color": "C0", "fontsize": 13},
        )

        ax[0].set(
            title="Transmissibilidade do deslocamento da base",
            xlabel="r=w/wn",
            ylabel="Td",
        )
        ax[1].set(
            title="Ângulo de fase entre os movimentos da base e da massa",
            xlabel="r=w/wn",
            ylabel=r"$\phi$",
        )

        fig.savefig("Output/Td and phi x r Plot.png")

        plt.clf()

        return

    def plot_transmissibility_of_forces(self) -> None:
        assert self.Y > 0, "Precisa ter um valor de Y"
        r = np.arange(0, 10, 0.01)
        Force = r**2 * np.sqrt(
            (1 + (2 * self.psi * r) ** 2)
            / ((1 - r**2) ** 2 + (2 * self.psi * r) ** 2)
        )

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 5), tight_layout=True)

        ax.plot(r, Force)

        ax.set(
            title="Transmissibilidade de Força para a base do sistema",
            xlabel="r=w/wn",
            ylabel="Tf",
        )

        fig.savefig("Output/Tf x r Plot.png")

        plt.clf()

        return

    def plot_transient_continuos(self) -> None:
        assert self.Fo > 0, "Precisa de uma força sendo aplicada"
        t = np.arange(0, self.wn * 5 / np.pi, 0.01)

        Xo = np.sqrt(self.xo**2 + (self.psi * self.wn * self.xo / self.wd) ** 2)
        phi_o = np.arctan(
            -1 * (self.vo + self.psi * self.wn * self.xo) / (self.wd * self.xo)
        )
        envelope = np.exp(-1 * self.psi * self.wn * t)
        big_X = self.deltast / (
            np.sqrt((1 - self.r**2) ** 2 + (2 * self.psi * self.r) ** 2)
        )
        if self.r == 1:
            phi = np.pi / 2
        else:
            phi = np.arctan(2 * self.psi * self.r / (1 - self.r**2))

        transient = Xo * envelope * np.cos(self.wd * t - phi_o)
        continuos = big_X * np.cos(self.w * t - phi)
        combined = transient + continuos

        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(20, 20), tight_layout=True)

        ax[0].plot(t, combined)
        ax[1].plot(t, transient)
        ax[2].plot(t, continuos)

        ax[0].set(
            title="Resposta do sistema",
            xlabel="t",
            ylabel="x(t)",
        )
        ax[1].set(
            title="Resposta das condições iniciais",
            xlabel="t",
            ylabel="x(t)",
        )
        ax[2].set(
            title="Resposta da força",
            xlabel="t",
            ylabel="x(t)",
        )

        fig.savefig("Output/Transient x Continuos.png")

        plt.clf()

        return

    def periodic_not_harmonic_force(
        self,
        piecewise: sym.Piecewise,
        inferior_limit: float,
        superior_limit: float,
        n_terms: int,
    ) -> None:
        t = sym.symbols("t")
        fourier_series = sym.fourier_series(
            piecewise,
            (t, inferior_limit, superior_limit),
        )

        if self.w == 0:
            period = superior_limit - inferior_limit
            self.w = 2 * np.pi / period
            self.r = self.w / self.wn

        time = np.arange(0, self.w * 5 / np.pi, 0.01)

        an, bn = [], []

        for n in range(1, n_terms):
            denominator = np.sqrt(
                (1 - n**2 * self.r**2) ** 2 + (2 * self.psi * n * self.r) ** 2
            )
            phi = np.arctan(2 * self.psi * n * self.r / (1 - n**2 * self.r**2))
            an.append(
                fourier_series.an.coeff(n).subs(t, n)
                * np.cos(n * self.w * time - phi)
                / (denominator * self.k)
            )
            bn.append(
                fourier_series.bn.coeff(n).subs(t, n)
                * np.sin(n * self.w * time - phi)
                / (denominator * self.k)
            )
        res_an = np.sum(np.array(an), axis=0)
        res_bn = np.sum(np.array(bn), axis=0)
        xp = res_an + res_bn + fourier_series.a0 / (2 * self.k)

        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 5), tight_layout=True)

        ax.plot(time, xp)

        ax.set(title="Comportamento do Sistema", xlabel="Tempo [s]", ylabel="f(t)")

        fig.savefig("Output/Periodic_not_harmonic_plot.png")

        plt.clf()

        return
