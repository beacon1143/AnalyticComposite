from Properties import linear_elasticity


class IsotropicMaterial:
    def __init__(self):
        self.elastic_props = linear_elasticity.LinearElasticity()
