import linear_elasticity
import linear_elasticity_aniso

le = linear_elasticity.linear_elasticity()
le.set_E_nu(1.0, 0.25)
print(le.get_K())

lea = linear_elasticity_aniso.linear_elasticity_aniso()
lea.print_Cijkl()
