import isotropic_material
import anisotropic_material

im = isotropic_material.isotropic_material()
im.elastic_props.set_E_nu(1.0, 0.25)
print(im.elastic_props.get_K())

am = anisotropic_material.anisotropic_material()
am.elastic_props.print_Cijkl()
