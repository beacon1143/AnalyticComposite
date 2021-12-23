from Composite import composite


class FiberComposite(composite.Composite):
    def compute_eff_elast(self):
        gam_fib = self._gamma
        gam_matr = 1.0 - gam_fib

        # auxiliary variable
        denom = gam_matr * self.matr.el_props.get_G() / (
            self.incl.el_props.get_K() + 0.33333333 * self.incl.el_props.get_G()
        ) + gam_fib * self.matr.el_props.get_G() / (
            self.matr.el_props.get_K() + 0.33333333 * self.matr.el_props.get_G()
        ) + 1.0

        # Young's modulus along fibers, E_L
        self.eff.el_props.E[0] = self.incl.el_props.get_E() * gam_fib + self.matr.el_props.get_E() * gam_matr + \
            4.0 * gam_fib * gam_matr * (self.incl.el_props.get_nu() - self.matr.el_props.get_nu()) ** 2 * \
            self.matr.el_props.get_G() / denom

        # Poisson's ratio, nu_LT
        self.eff.el_props.nu[0, 1] = gam_fib * self.incl.el_props.get_nu() + gam_matr * self.matr.el_props.get_nu() + \
            4.0 * gam_fib * gam_matr * (self.incl.el_props.get_nu() - self.matr.el_props.get_nu()) * (
                self.matr.el_props.get_G() / (
                    self.matr.el_props.get_K() + 0.33333333 * self.matr.el_props.get_G()
                ) - self.matr.el_props.get_G() / (
                    self.incl.el_props.get_K() + 0.33333333 * self.incl.el_props.get_G()
                )
            ) / denom

        self.eff.el_props.nu[0, 0] = self.eff.el_props.nu[0, 1]

        # bulk modulus
        k23 = self.matr.el_props.get_K() + 0.33333333 * self.matr.el_props.get_G() + gam_fib / (
            1.0 / (self.incl.el_props.get_K() - self.matr.el_props.get_K() + 0.33333333 * (
                self.incl.el_props.get_G() - self.matr.el_props.get_G()
            )) + gam_matr / (
                self.matr.el_props.get_K() + 1.33333333 * self.matr.el_props.get_G()
            )
        )

        # shear module, G_LT
        self.eff.el_props.G[0] = self.matr.el_props.get_G() * (
            (1.0 + gam_fib) * self.incl.el_props.get_G() + (1.0 - gam_matr) * self.matr.el_props.get_G()
        ) / (
            (1.0 - gam_fib) * self.incl.el_props.get_G() + (1.0 + gam_matr) * self.matr.el_props.get_G()
        )

        self.eff.el_props.G[1] = self.eff.el_props.G[0]

        # shear module, G_TT
        self.eff.el_props.G[2] = self.matr.el_props.get_G() * (1.0 + gam_fib / (
            self.matr.el_props.get_G() / (self.incl.el_props.get_G() - self.matr.el_props.get_G()) +
            (self.matr.el_props.get_K() + 2.33333333 * self.matr.el_props.get_G()) /
            (2.0 * self.matr.el_props.get_K() + 2.66666667 * self.matr.el_props.get_G())
        ))

        # auxiliary variables
        coef = 4.0 * self.eff.el_props.nu[0, 0] ** 2 * self.eff.el_props.G[2] * k23
        frac = coef / self.eff.el_props.E[0]

        # Young's modulus across fibers, E_T
        self.eff.el_props.E[1] = 4.0 * self.eff.el_props.G[2] * k23 / (
            k23 + self.eff.el_props.G[2] + frac
        )

        self.eff.el_props.E[2] = self.eff.el_props.E[1]

        # Poisson's ratio, nu_TT
        self.eff.el_props.nu[1, 1] = (k23 - self.eff.el_props.G[2] - frac) / (
            k23 + self.eff.el_props.G[2] + frac
        )

        self.eff.el_props.nu[2, 1] = self.eff.el_props.nu[1, 1]

        # Poisson's ratio, nu_TL
        self.eff.el_props.nu[1, 0] = 4.0 * coef / (
            self.eff.el_props.E[0] * (k23 + self.eff.el_props.G[2]) + coef
        )

        self.eff.el_props.nu[2, 0] = self.eff.el_props.nu[1, 0]

        self.eff.el_props.compute_Cijkl_from_ortho()
