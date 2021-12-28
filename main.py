# from Material import isotropic_material, anisotropic_material
from Composite import fiber_composite, dispersed_composite

try:
    # fc = fiber_composite.FiberComposite(0.1)
    # fc.incl.el_props.set_E_nu(2000.0, 0.2)
    # fc.matr.el_props.set_E_nu(2.0, 0.3)
    # fc.compute_eff_elast()
    # fc.eff.el_props.print_Cijkl()

    dc = dispersed_composite.DispersedComposite(0.05)
    dc.incl.el_props.set_E_nu(1.0e-6, 0.25)
    dc.matr.el_props.set_E_nu(1.0, 0.4)
    print('K_matr =', dc.matr.el_props.get_K())
    dc.compute_eff_elast()
except BaseException as e:
    print(e)
