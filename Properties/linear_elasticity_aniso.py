import sys
import numpy


class LinearElasticityAniso:
    def __init__(self):
        self.__Cijkl = numpy.zeros((6, 6), dtype=float)    # stiffness tensor
        self.__E = numpy.zeros(3, dtype=float)             # Young's moduli
        self.__nu = numpy.zeros((3, 2), dtype=float)       # Poisson's ratios
        self.__G = numpy.zeros(3, dtype=float)             # shear moduli

    def compute_Cijkl_from_ortho(self):
        c1 = 1.0 / self.__E[0]
        c2 = 1.0 / self.__E[1]
        c3 = 1.0 / self.__E[2]

        c12 = -self.__nu[0, 0] / self.__E[0]
        c13 = -self.__nu[0, 1] / self.__E[0]
        c23 = -self.__nu[1, 1] / self.__E[1]

        det = c1 * c2 * c3 - c1 * c23 * c23 - c12 * c12 * c3 + 2.0 * c12 * c13 * c23 - c13 * c13 * c2

        self.__Cijkl[0, 0] = (c2 * c3 - c23 * c23) / det      # C_1111
        self.__Cijkl[1, 1] = (c1 * c3 - c13 * c13) / det      # C_2222
        self.__Cijkl[2, 2] = (c1 * c2 - c12 * c13) / det      # C_3333

        self.__Cijkl[0, 1] = (-c12 * c3 + c13 * c23) / det    # C_1122
        self.__Cijkl[1, 0] = self.__Cijkl[0, 1]

        self.__Cijkl[0, 2] = (c12 * c23 - c13 * c2) / det     # C_1133
        self.__Cijkl[2, 0] = self.__Cijkl[0, 2]

        self.__Cijkl[1, 2] = (-c1 * c23 + c12 * c13) / det    # C_2233
        self.__Cijkl[2, 1] = self.__Cijkl[1, 2]

        self.__Cijkl[3, 3] = self.__G[0]                             # C_1212
        self.__Cijkl[4, 4] = self.__G[1]                             # C_1313
        self.__Cijkl[5, 5] = self.__G[2]                             # C_2323

    def compute_ortho_from_Cijkl(self):
        Dij = numpy.zeros((6, 6), dtype=float)
        for i in range(0, 6):
            Dij[i, i] = self.__Cijkl[i, i]
        Dij[0, 1] = 0.5 * (self.__Cijkl[0, 1] + self.__Cijkl[1, 0])
        Dij[0, 2] = 0.5 * (self.__Cijkl[0, 2] + self.__Cijkl[2, 0])
        Dij[1, 2] = 0.5 * (self.__Cijkl[1, 2] + self.__Cijkl[2, 1])
        Dij[1, 0] = Dij[0, 1]
        Dij[2, 0] = Dij[0, 2]
        Dij[2, 1] = Dij[1, 2]
        # for 2D case, when G13 = G23 = 0
        if Dij[4, 4] == 0.0:
            Dij[4, 4] = sys.float_info.max
        if Dij[5, 5] == 0.0:
            Dij[5, 5] = sys.float_info.max

        Sijkl = numpy.linalg.inv(Dij)    # compliance tensor

        for i in range(0, 3):
            self.__E[i] = 1.0 / Sijkl[i, i]
            self.__G[i] = 1.0 / Sijkl[i + 3, i + 3]
            for j in range(0, 2):
                self.__nu[i, j] = -self.__E[i] * Dij[i, j + (j >= i)]

    def print_Cijkl(self):
        for i in range(0, 6):
            for j in range(0, 6):
                print('C_', i + 1, j + 1, ' = ', self.__Cijkl[i, j], sep='')

    def print_ortho(self):
        for i in range(0, 3):
            print('E_', i + 1, ' = ', self.__E[i], sep='')
        for i in range(0, 3):
            for j in range(0, 2):
                print('nu_', i + 1, j + 1 + (j >= i), ' = ', self.__nu[i, j], sep='')
        print('G_12 = ', self.__G[0], sep='')
        print('G_13 = ', self.__G[1], sep='')
        print('G_23 = ', self.__G[2], sep='')
