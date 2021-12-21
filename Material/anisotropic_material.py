from Properties import linear_elasticity_aniso


class AnisotropicMaterial:
    def __init__(self):
        self.elastic_props = linear_elasticity_aniso.LinearElasticityAniso()
