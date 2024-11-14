import numpy as np

def earth_pressures(
    friction: float=30, cohesion: float=10,
    adhesion: float=0,
    delta: float=0, beta: float=0, theta: float=0,
    active: bool=False,
    surcharge: float=0, inclination: float=0,
    log: bool=False,
    ):
    """Computes earth pressure coefficients according to Eurocode 7 - Part 1 - Annex D.

    Args:
        friction (float, optional): Friction angle, phi (in degrees). Defaults to 30.
        cohesion (float, optional): Cohesion, c. Defaults to 10.
        adhesion (float, optional): Adhesion, a (only for undrained conditions). Defaults to 0.
        delta (float, optional): Wall-soil friction angle (in degrees). Defaults to 0.
        beta (float, optional): Slope of soil (in degrees), positive when rising away from the wall. Defaults to 0.
        theta (float, optional): Inclination of wall (in degrees), positive when soil overhangs the wall. Defaults to 0.
        surcharge (float, optional): Surcharge, q. Defaults to 0.
        inclination (float, optional): Inclination of surcharge (in degrees), positive when towards the wall. Defaults to 0.
        active (bool, optional): Compute active pressure. Defaults to False.
        log (bool, optional): Print log. Defaults to False.

    Returns:
        float: Bearing resistance, R.
    """
    # Check inputs
    if friction < 0 or friction > 90: raise ValueError("Friction angle must be between 0 and 90 degrees.")
    if cohesion < 0: raise ValueError("Cohesion must be non-negative.")
    if delta < 0 or delta > 90: raise ValueError("Wall-soil friction angle must be between 0 and 90 degrees.")
    if beta < -90 or beta > 90: raise ValueError("Slope must be between -90 and 90 degrees.")
    if theta < -90 or theta > 90: raise ValueError("Inclination must be between -90 and 90 degrees.")
    if surcharge < 0: raise ValueError("Surcharge must be non-negative.")
    if inclination < -90 or inclination > 90: raise ValueError("Inclination must be between -90 and 90 degrees.")
    # Soil
    _f = np.deg2rad(friction)
    _c = cohesion
    # Angles
    _d = np.deg2rad(delta)
    _b = np.deg2rad(beta)
    _t = np.deg2rad(theta)
    # Adhesion
    _a = adhesion if _f==0 else np.tan(_d)/np.tan(_f)
    # Active
    if active:
        _f = -_f
        _c = -_c
        _d = -_d
    # Surcharge
    _q = surcharge
    _i = np.deg2rad(inclination)
    _p = _q * np.cos(_i)
    if (_c==0 and (_q==0 or _i==0)) or active: _bo = _b
    elif _f==0: _bo = 0
    else: _bo = np.arctan( _q*np.sin(_b+_i) / ( _q*np.cos(_b+_i) + _c/np.tan(_f) ) )
    # Boundary conditions
    if _f==0:
        _mt = np.arccos(- _p/_c * np.sin(_b) * np.cos(_b)) / 2
        _mw = np.arccos(_a/_c) / 2
    else:
        _mt = ( np.arccos(- np.sin(_bo) / np.sin(_f)) - _f - _bo ) / 2
        _mw = ( np.arccos(np.sin(_d) / np.sin(_f)) - _f - _d ) / 2
    # Tangent rotation
    _v = _mt + _b - _mw - _t
    # Coefficients
    _Kn = (1 + np.sin(_f)*np.sin(2*_mw+_f)) / (1 - np.sin(_f)*np.sin(2*_mt+_f)) * np.exp(2*_v*np.tan(_f))
    _Kq = (np.cos(_b))**2 if _f==0 else _Kn * (np.cos(_b))**2
    _Kc = 2*_v + np.sin(2*_mt) + np.sin(2*_mw) if _f==0 else (_Kn - 1) / np.tan(_f)
    _Kg = np.cos(_t) + np.sin(_b)*np.cos(_mw)/np.sin(_mt) if _f==0 else _Kn * np.cos(_b) * np.cos(_b-_t)
    # Log
    if log:
        print("Input parameters:")
        print("phi =", friction, 'deg')
        print("c =", _c, 'kPa')
        print("a =", _a, 'kPa')
        print("delta =", delta, 'deg')
        print("beta =", beta, 'deg')
        print("theta =", theta, 'deg')
        print("q =", _q, 'kPa')
        print("p =", _p, 'kPa')
        print("\nBoundary conditions:")
        print("mw =", np.rad2deg(_mw), 'deg')
        print("v =", np.rad2deg(_v), 'deg')
        print("mt =", np.rad2deg(_mt), 'deg')
        print("\nCoefficients:")
        print("Kn =", _Kn)
        print("Kq =", _Kq) 
        print("Kc =", _Kc)
        print("Kg =", _Kg)
    return (_Kn,_Kq,_Kc,_Kg)