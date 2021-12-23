from Material import isotropic_material, anisotropic_material


class Composite:
    def __init__(self, gamma):
        if gamma < 0.0 or gamma > 1.0:
            raise ValueError('Error! Concentration of inclusions must be greater than zero and less than one!\n')
        self._gamma = gamma
        self.matr = isotropic_material.IsotropicMaterial()
        self.incl = isotropic_material.IsotropicMaterial()
        self.eff = anisotropic_material.AnisotropicMaterial()

    def compute_eff_elast(self):
        raise NotImplementedError('Error! Function calc_eff_elast is not implemented in base class Composite!\n')
