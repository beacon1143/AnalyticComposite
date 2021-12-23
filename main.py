# from Material import isotropic_material, anisotropic_material
from Composite import fiber_composite

# im = isotropic_material.IsotropicMaterial()
# im.elastic_props.set_E_nu(1.0, 0.25)
# print(im.elastic_props.get_K())
#
# am = anisotropic_material.AnisotropicMaterial()
# am.elastic_props.print_Cijkl()

try:
    fc = fiber_composite.FiberComposite(0.1)
    fc.incl.el_props.set_E_nu(2000.0, 0.2)
    fc.matr.el_props.set_E_nu(2.0, 0.3)
    fc.compute_eff_elast()
    fc.eff.el_props.print_Cijkl()
except BaseException as e:
    print(e)
