from Material import isotropic_material, anisotropic_material


class Composite:
    def __init__(self, gamma):
        self.__gamma = gamma
        self.matr = isotropic_material.IsotropicMaterial()
        self.incl = isotropic_material.IsotropicMaterial()
        self.effective = anisotropic_material.AnisotropicMaterial()

    def calc_eff_elast(self):
        raise NotImplementedError('Error! Function calc_eff_elast is not implemented in base class Composite!\n')
