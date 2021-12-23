from Properties import linear_elasticity


class IsotropicMaterial:
    def __init__(self):
        self.el_props = linear_elasticity.LinearElasticity()
