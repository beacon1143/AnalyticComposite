from Material import isotropic_material
from Composite import composite


class DispersedComposite(composite.Composite):
    def __init__(self, gamma):
        super().__init__(gamma)
        self.__eff2 = isotropic_material.IsotropicMaterial()

    def compute_eff_elast(self):
        nu_m = self.matr.el_props.get_nu()

        # auxiliary variable
        frac = self.incl.el_props.get_G() / self.matr.el_props.get_G()

        G_eff = self.matr.el_props.get_G() * (
            1.0 - 15.0 * (1.0 - nu_m) * (1.0 - frac) * self._gamma / (
                7.0 - 5.0 * nu_m + 2.0 * (4.0 - 5.0 * nu_m) * frac
            )
        )

        # auxiliary variable
        diff = self.incl.el_props.get_K() - self.matr.el_props.get_K()

        K_eff = self.matr.el_props.get_K() + (
            self._gamma * diff / (
                1.0 + diff / (self.matr.el_props.get_K() + 4.0 * self.matr.el_props.get_G() / 3.0)
            )
        )
        print('K_eff =', K_eff)

        self.__eff2.el_props.set_K_G(K_eff, G_eff)

        for i in range(0, 3):
            self.eff.el_props.E[i] = self.__eff2.el_props.get_E()
            self.eff.el_props.G[i] = self.__eff2.el_props.get_G()
            for j in range(0, 2):
                self.eff.el_props.nu[i, j] = self.__eff2.el_props.get_nu()

        self.eff.el_props.compute_Cijkl_from_ortho()
