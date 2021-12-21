class LinearElasticity:
    def __init__(self):
        self.__E = 0.0        # Young's modulus
        self.__nu = 0.0       # Poisson's ratio
        self.__la = 0.0
        self.__G = 0.0        # shear modulus
        self.__K = 0.0        # bulk modulus

    def set_lambda_G(self, lambd, G):
        self.__la = lambd
        self.__G = G
        self.__K = lambd + 0.66666667 * G
        self.__E = G * ((3.0 * lambd + 2.0 * G) / (lambd + G))
        self.__nu = 0.5 * lambd / (lambd + G)

    def set_E_G(self, E, G):
        self.__E = E
        self.__G = G
        self.__K = 0.33333333 * E * G / (3.0 * G - E)
        self.__la = G * (E - 2.0 * G) / (3.0 * G - E)
        self.__nu = 0.5 * E / G - 1.0

    def set_K_lambda(self, K, lambd):
        self.__K = K
        self.__la = lambd
        self.__E = 9.0 * K * (K - lambd) / (3.0 * K - lambd)
        self.__G = 1.5 * (K - lambd)
        self.__nu = lambd / (3.0 * K - lambd)

    def set_K_G(self, K, G):
        self.__K = K
        self.__G = G
        self.__E = 9.0 * K * G / (3.0 * K + G)
        self.__la = K - 0.66666667 * G
        self.__nu = (1.5 * K - G) / (3.0 * K + G)

    def set_lambda_nu(self, lambd, nu):
        self.__la = lambd
        self.__nu = nu
        self.__K = 0.33333333 * lambd * (1.0 + nu) / nu
        self.__E = lambd * (1.0 + nu) * (1.0 - 2.0 * nu) / nu
        self.__G = 0.5 * lambd * (1.0 - 2.0 * nu) / nu

    def set_G_nu(self, G, nu):
        self.__G = G
        self.__nu = nu
        self.__K = 0.66666667 * G * (1.0 + nu) / (1.0 - 2.0 * nu)
        self.__E = 2.0 * G * (1.0 + nu)
        self.__la = 2.0 * G * nu / (1 - 2 * nu)

    def set_E_nu(self, E, nu):
        self.__E = E
        self.__nu = nu
        self.__K = 0.33333333 * E / (1.0 - 2.0 * nu)
        self.__la = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
        self.__G = 0.5 * E / (1.0 + nu)

    def set_K_nu(self, K, nu):
        self.__K = K
        self.__nu = nu
        self.__E = 3.0 * K * (1.0 - 2.0 * nu)
        self.__la = 3.0 * K * nu / (1.0 + nu)
        self.__G = 1.5 * K * (1.0 - 2.0 * nu) / (1.0 + nu)

    def set_K_E(self, K, E):
        self.__K = K
        self.__E = E
        self.__la = 3.0 * K * (3.0 * K - E) / (9.0 * K - E)
        self.__G = 3.0 * K * E / (9.0 * K - E)
        self.__nu = 0.16666667 * (3.0 * K - E) / K

    def get_E(self):
        return self.__E

    def get_nu(self):
        return self.__nu

    def get_lambda(self):
        return self.__la

    def get_G(self):
        return self.__G

    def get_K(self):
        return self.__K
