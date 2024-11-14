import numpy as np

def bearing_resistance(
    width: float, length: float=0,
    embedment: float=0, inclination: float=0,
    unit_weight: float=20, friction: float=30, cohesion: float=10,
    vertical: float=None,
    horizontal: float=0, longitudinal_force: float=0,
    moment: float=0, longitudinal_moment: float=0,
    surcharge: float=0,
    log: bool=False,
    ):
    """Computes bearing resistance according to Eurocode 7 - Part 1 - Annex D.

    Args:
        width (float): Foundation width, B.
        length (float, optional): Foundation length, L. Use zero for strip. Defaults to 0.
        embedment (float, optional): Embedment depth, D. Defaults to 0.
        inclination (float, optional): Inclination, alpha (in degrees). Defaults to 0.
        unit_weight (float, optional): Effective weight density, gamma. Defaults to 20.
        friction (float, optional): Friction angle, phi (in degrees). Defaults to 30.
        cohesion (float, optional): Cohesion. Defaults to 10.
        vertical (float, optional): Vertical load, V. Use None when horizontal loads and moments are not specified. Defaults to None.
        horizontal (float, optional): Horizontal load, HB, in width direction. Defaults to 0.
        longitudinal_force (float, optional): Horizontal load, HL, in length direction. Defaults to 0.
        moment (float, optional): Moment, MB, about length axis. Defaults to 0.
        longitudinal_moment (float, optional): Moment, ML, about width axis. Defaults to 0.
        surcharge (float, optional): Surcharge pressure, q, in addition to gamma*D. Defaults to 0.
        log (bool, optional): Print log. Defaults to False.

    Returns:
        float: Bearing resistance, R.
    """
    # Check inputs
    if width <= 0: raise ValueError("Width must be positive.")
    if length < 0: raise ValueError("Length must be non-negative.")
    if embedment < 0: raise ValueError("Embedment must be non-negative.")
    if inclination < 0: raise ValueError("Inclination must be non-negative.")
    if inclination >= 90: raise ValueError("Inclination must be less than 90 degrees.")
    if unit_weight <= 0: raise ValueError("Unit weight must be positive.")
    if friction < 0: raise ValueError  ("Friction angle must be non-negative.")
    if cohesion < 0: raise ValueError("Cohesion must be non-negative.")
    if vertical is not None and vertical < 0: raise ValueError("Vertical load must be positive.")
    if horizontal < 0: raise ValueError("Horizontal load must be non-negative.")
    if longitudinal_force < 0: raise ValueError("Longitudinal force must be non-negative.")
    if moment < 0: raise ValueError("Moment must be non-negative.")
    if longitudinal_moment < 0: raise ValueError("Longitudinal moment must be non-negative.")
    if surcharge < 0: raise ValueError("Surcharge must be non-negative.")
    # Geometry
    _B = width
    _L = length
    _D = embedment
    _a = np.deg2rad(inclination)
    # Soil
    _gamma = unit_weight
    _phi = np.deg2rad(friction)
    _c = cohesion
    # Load
    _V = vertical
    _HB = horizontal
    _HL = longitudinal_force
    _H = 0 if _V==None else np.sqrt(_HB**2 + _HL**2)
    _MB = moment
    _ML = longitudinal_moment
    # Surcharge
    _q = surcharge + _gamma*_D
    # Effective Dimensions
    if _V is None:
        _eB = 0
        _eL = 0
    else:
        _eB = _MB / _V
        _eL = _ML / _V
    _Be = _B - 2*_eB
    _Le = _L - 2*_eL
    _Ae = _Be if _L==0 else _Be * _Le
    # Bearing Resistance
    _Nq = np.exp(np.pi * np.tan(_phi)) * np.tan(np.pi/4+_phi/2)**2
    _Nc = np.pi+2 if _phi==0 else (_Nq-1) / np.tan(_phi)
    _Ng = 2 * (_Nq-1) * np.tan(_phi)
    # Inclination of the Foundation Base
    _bq = (1 - _a * np.tan(_phi))
    _bg = _bq
    _bc = 1 - 2*_a / (np.pi+2) if _phi==0 else _bq - (1 - _bq) / (_Nc*np.tan(_phi))
    # Shape of the Foundation
    if _Le == 0:
        _sq = 1
        _sg = 1
        _sc = 1
    else:
        _sq = 1 + (_Be/_Le) * np.sin(_phi)
        _sg = 1 - 0.3*(_Be/_Le)
        _sc = 1 + 0.2*(_Be/_Le) if _phi==0 else (_sq*_Nq-1) / (_Nq-1)
    # Inclination of the Load
    if _H==0:
        _iq = 1
        _ig = 1
        _ic = 1
    else:
        _mB = (2 + _Be/_Le)/(1 + _Be/_Le)
        _mL = (2 + _Le/_Be)/(1 + _Le/_Be)
        _theta = 0 if np.max(_H)==0 else np.arctan2(_HB,_HL)
        _m = _mL * np.cos(_theta)**2 + _mB * np.sin(_theta)**2
        _iq = 1 if _phi==0 else (1 - _H/(_V+_Ae*_c/np.tan(_phi)))**_m
        _ig = 1 if _phi==0 else (1 - _H/(_V+_Ae*_c/np.tan(_phi)))**(_m+1)
        _ic = 0.5 * (1 + np.sqrt(1-_H/(_Ae*_c))) if _phi==0 else _iq - (1-_iq)/(_Nc*np.tan(_phi))
    # Result
    _qu = _c*_Nc*_bc*_sc*_ic + _q*_Nq*_bq*_sq*_iq + 0.5*_gamma*_Be*_Ng*_bg*_sg*_ig
    _R = _Ae*_qu
    # Log
    if log:
        print('Geometry:')
        print('B =', _B, 'm')
        print('L =', _L, 'm')
        print('D =', _D, 'm')
        print('a =', inclination, 'deg')
        print('\nSoil Properties:')
        print('gamma =', _gamma, 'kN/m3')
        print('phi =', friction, 'deg')
        print('c =', _c, 'kPa')
        print('\nLoading:')
        print('V =', _V, 'kN')
        print('HB =', _HB, 'kN')
        print('HL =', _HL, 'kN')
        print('H =', _H, 'kN')
        print('MB =', _MB, 'kNm')
        print('ML =', _ML, 'kNm')
        print('q =', _q, 'kPa')
        print('\nEffective Dimensions:')
        print('eB =', _eB, 'm')
        print('eL =', _eL, 'm')
        print('Be =', _Be, 'm')
        print('Le =', _Le, 'm')
        print('Ae =', _Ae, 'm2')
        print('\nBearing Capacity Factors:')
        print('Nq =', _Nq)
        print('Nc =', _Nc)
        print('Ng =', _Ng)
        print('\nInclination Factors:')
        print('bq =', _bq)
        print('bc =', _bc)
        print('bg =', _bg)
        print('\nShape Factors:')
        print('sq =', _sq)
        print('sc =', _sc)
        print('sg =', _sg)
        print('\nLoad Inclination Factors:')
        print('iq =', _iq)
        print('ic =', _ic)
        print('ig =', _ig)
        print('\nResults:')
        print('qu =', _qu, 'kPa')
        print('R =', _R, 'kN')
    return _R