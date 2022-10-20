from math import exp, log, sqrt
from pytherm import constants
from pytherm.activity.activity_model import ActivityModel
import numpy as np
from .eos import EOS
R = constants.R


class PSRK(EOS):
    """Class for solving the Predictive Soave-Redlich-Kwong equation"""

    def __init__(self,
                 system: dict[str, float],
                 cr_params,
                 ms_params,
                 activity_model: ActivityModel,
                 omegas=None) -> None:
        """PSRK class constructor

        Parameters
        ----------
        system : dict[str, float]
            Dictionary with component concentrations
        cr_params : dict
            Critical parameters
        ms_params : dict
            Dictionary with Mathias-Copeman parameters
        activity_model : Activity_model
            Model for activity coefficient calculation,
            UNIFAC in original model
        omegas : dict
            Acentric factors
        """
        self.ms_params = ms_params
        self.system = system
        self.activity_model = activity_model
        self.cr_params = cr_params
        self.omegas = omegas

        self.bi = self.get_bi()

    def get_p(self, system: dict[str, float], T: float, V: float) -> float:
        """Calculate p using PSRK equation

        Parameters
        ----------
        system : dict[str, float]
            Dictionary with component concentrations
        T : float
            Temperature, [K]
        V : float
            Molar volume, [m^3/mol]

        Returns
        -------
        float
            Pressure, [Pa]
        """
        alphas = self.get_alphas(T)
        ai = self.get_ai(alphas)
        ge = self.get_ge(system, T)

        b = self.get_b(system)
        a = self.get_a(system, T, b, ai, ge)
        return (R * T) / (V - b) - a / (V * (V + b))

    def get_alphas(self, T: float) -> dict[str, float]:
        """Calculate alpha coefficients using Mathias-Copeman equation

        Parameters
        ----------
        T : float
            Temperature, [K]

        Returns
        -------
        dict[str, float]
            alpha coefficients
        """
        alphas = {}
        for i in self.system:
            Tr = T / self.cr_params[i]['Tc']
            sqrt_T = sqrt(Tr)
            if Tr < 1:
                c1 = self.ms_params[i]['c1']
                c2 = self.ms_params[i]['c2']
                c3 = self.ms_params[i]['c3']
                alphas[i] = (1
                             + c1 * (1 - sqrt_T)
                             + c2 * (1 - sqrt_T) ** 2
                             + c3 * (1 - sqrt_T) ** 3) ** 2
            else:
                c1 = self.ms_params[i]['c1']
                alphas[i] = (1
                             + c1 * (1 - sqrt_T)) ** 2
        return alphas

    def get_ai(self, alphas: dict[str, float]) -> dict[str, float]:
        """Calculate a_i coefficient for each component in system

        Parameters
        ----------
        alphas : dict[str, float]
            Mathias-Copeman alpha coefficients

        Returns
        -------
        dict[str, float]
            a_i coefficient for each component in system
        """
        ai = {}
        for i in self.system:
            Tc = self.cr_params[i]['Tc']
            Pc = self.cr_params[i]['Pc']
            ai[i] = 0.42748 * R ** 2 * Tc ** 2 * alphas[i] / Pc
        return ai

    def get_bi(self) -> dict[str, float]:
        """Calculate b_i coefficient for each component in system.
        Calculated only once during the initialization of the object.

        Returns
        -------
        dict[str, float]
            b_i coefficient for each component in system
        """
        bi = {}
        for i in self.system:
            Tc = self.cr_params[i]['Tc']
            Pc = self.cr_params[i]['Pc']
            bi[i] = 0.08664 * R * Tc / Pc
        return bi

    def get_a(self, system: dict[str, float],
              T: float,
              b: float,
              ai: dict[str, float],
              ge: float,
              A=-0.64663) -> float:
        """Calculate a parameter for SRK using modified Huron-Vidal mixing rule

        Parameters
        ----------
        system : dict[str, float]
            Dictionary with component concentrations
        T : float
            Temperature, [K]
        b : float
            b parameter for SRK
        ai : dict[str, float]
            a_i parameters dictionary for SRK
        ge : float
            Excess Gibbs free energy
        A : float, optional
            MHV constant, by default -0.64663

        Returns
        -------
        float
            a parameter for SRK
        """
        bi = self.bi
        s1 = 0
        s2 = 0
        for i in self.system:
            s1 += system[i] * ai[i] / bi[i]
            s2 += system[i] * log(b / bi[i])

        return b * (ge / A + s1 + R * T / A * s2)

    def get_b(self, system: dict[str, float]) -> float:
        """Calculate b parameter for SQR equation using linear mixing rule

        Parameters
        ----------
        system : dict[str, float]
            Dictionary with component concentrations

        Returns
        -------
        float
            b parameter for SQR equation
        """
        bi = self.bi
        b = 0
        for i in self.system:
            b += system[i] * bi[i]
        return b

    def get_ge(self, system: dict[str, float], T: float) -> float:
        """_summary_

        Parameters
        ----------
        system : dict[str, float]
            Dictionary with component concentrations
        T : float
            Temperature, [K]

        Returns
        -------
        float
            Excess Gibbs free energy
        """
        return self.activity_model.get_ge(system, T)

    def get_cubic_coef(self, system: dict[str, float], T: float, P: float) -> tuple:
        """Return coefficient in cubic form of SQR V^3 + a1 * V^2 + a2*V + a3 =0

        Parameters
        ----------
        system : dict[str, float]
            Dictionary with component concentrations
        T : float
            Temperature, [K]
        P : float
            Pressure, [Pa]

        Returns
        -------
        tuple
            (a1, a2, a3)
        """
        alphas = self.get_alphas(T)
        ai = self.get_ai(alphas)
        ge = self.get_ge(system, T)

        b = self.get_b(system)
        a = self.get_a(system, T, b, ai, ge)

        return (1,
                - R * T / P,
                a / P - b * R * T / P - b ** 2,
                - a * b / P)

    def get_roots(self, system: dict[str, float], P: float, T: float) -> tuple:
        """Solving PSRK equation for P and T

        Parameters
        ----------
        system : dict[str, float]
            Dictionary with component concentrations
        P : float
            Pressure, [Pa]
        T : float
            Temperature, [K]

        Returns
        -------
        tuple
            One or two roots
        """
        cs = self.get_cubic_coef(system=system, T=T, P=P)
        r = np.roots(cs)
        out = []
        for i in r:
            if i.imag == 0:
                out.append(i)
        if len(out) == 1:
            return [float(out[0])]
        else:
            return [float(min(r)), float(max(r))]

    def get_f(self,
              system: dict[str, float],
              P: float,
              V: float,
              T: float) -> dict[str, float]:
        """Calculate fugacity coefficients for given P, V, T

        Parameters
        ----------
        system : dict[str, float]
            Dictionary with component concentrations
        P : float
            Pressure, [Pa]
        V : float
            Molar volume, [m^3/mol]
        T : float
            Temperature, [K]

        Returns
        -------
        dict[str, float]
            Dictionary with fugacity coefficients for all components
        """
        alphas = self.get_alphas(T)
        ai = self.get_ai(alphas)
        bi = self.get_bi()

        alp = {}
        for i in system:
            alp[i] = ai[i] / (bi[i] * R * T)

        b = self.get_b(system)

        yi = self.activity_model.get_y(system, T)
        ge_RT = self.activity_model.get_ge_RT(system, T)
        al = self.__get_al(system, ge_RT, b, bi, alp)

        der_ai = self.__get_der_ai(yi, b, bi, alp)

        def f(v):
            lnf = {}
            for i in system:
                lnf[i] = (log(R * T / (P * (v - b)))
                          + (1 / (v - b) - al/(v + b)) * bi[i]
                          - der_ai[i] * log((v + b) / v))
            fi = {}
            for i in lnf:
                fi[i] = exp(lnf[i])
            return fi
        return f(V)

    def __get_der_ai(self, y, b, bi, alp, A=-0.64663):
        der_ai = {}
        for i in y:
            der_ai[i] = 1 / A * (log(y[i]) + log(b / bi[i]) +
                                 bi[i] / b - 1) + alp[i]
        return der_ai

    def __get_al(self, system, ge_RT, b, bi, alp, A=-0.64663):
        s1 = 0
        s2 = 0
        for i in bi:
            s1 += system[i] * log(b / bi[i])
            s2 += system[i] * alp[i]
        return 1 / A * (ge_RT + s1) + s2

    def get_system(self):
        return self.system
