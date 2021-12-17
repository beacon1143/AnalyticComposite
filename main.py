import linear_elasticity

le = linear_elasticity.linear_elasticity()
le.set_E_nu(1.0, 0.25)
print(le.get_K())
