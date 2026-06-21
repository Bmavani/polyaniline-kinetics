
def initial_solverfinal(t, ANI2, k1, k_prime_2, APS0, ANI0):
    """
    ODE right-hand side for ANI polymerization.

    Parameters
    t         : current time (min)
    ANI2      : current ANI concentration (mol/L)
    k1        : propagation rate constant
    k_prime_2 : reverse/termination rate constant
    APS0      : initial APS (initiator) concentration (mol/L)
    ANI0      : initial ANI concentration (mol/L)

    Returns
    dANI2_dt : rate of change of ANI concentration
    """
    beta_1 = 0.8 * k1 - k_prime_2
    beta_2 = (k1 * APS0) - (0.8 * k1 * ANI0) + (k_prime_2 * ANI0)

    dANI2_dt = -ANI2 * ((ANI2 * beta_1) + beta_2)

    return dANI2_dt
