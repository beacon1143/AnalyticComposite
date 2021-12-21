# from Material import isotropic_material, anisotropic_material
from Composite import composite

# im = isotropic_material.IsotropicMaterial()
# im.elastic_props.set_E_nu(1.0, 0.25)
# print(im.elastic_props.get_K())
#
# am = anisotropic_material.AnisotropicMaterial()
# am.elastic_props.print_Cijkl()

try:
    c = composite.Composite(0.1)
    c.calc_eff_elast()
except BaseException as e:
    print(e)
