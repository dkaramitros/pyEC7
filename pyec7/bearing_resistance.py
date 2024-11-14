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
    _HB = np.abs(horizontal)
    _HL = np.abs(longitudinal_force)
    _H = 0 if _V==None else np.sqrt(_HB**2 + _HL**2)
    _MB = np.abs(moment)
    _ML = np.abs(longitudinal_moment)
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
        _mB = 2 if _Le==0 else (2 + _Be/_Le)/(1 + _Be/_Le)
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
        print('B =', f'{_B:.4g}', 'm')
        print('L =', f'{_L:.4g}', 'm')
        print('D =', f'{_D:.4g}', 'm')
        print('a =', f'{inclination:.4g}', 'deg')
        print('\nSoil Properties:')
        print('gamma =', f'{_gamma:.4g}', 'kN/m3')
        print('phi =', f'{friction:.4g}', 'deg')
        print('c =', f'{_c:.4g}', 'kPa')
        print('\nLoading:')
        print('V =', f'{_V:.4g}', 'kN')
        print('HB =', f'{_HB:.4g}', 'kN')
        print('HL =', f'{_HL:.4g}', 'kN')
        print('H =', f'{_H:.4g}', 'kN')
        print('MB =', f'{_MB:.4g}', 'kNm')
        print('ML =', f'{_ML:.4g}', 'kNm')
        print('q =', f'{_q:.4g}', 'kPa')
        print('\nEffective Dimensions:')
        print('eB =', f'{_eB:.4g}', 'm')
        print('eL =', f'{_eL:.4g}', 'm')
        print('Be =', f'{_Be:.4g}', 'm')
        print('Le =', f'{_Le:.4g}', 'm')
        print('Ae =', f'{_Ae:.4g}', 'm2')
        print('\nBearing Capacity Factors:')
        print('Nq =', f'{_Nq:.4g}')
        print('Nc =', f'{_Nc:.4g}')
        print('Ng =', f'{_Ng:.4g}')
        print('\nInclination Factors:')
        print('bq =', f'{_bq:.4g}')
        print('bc =', f'{_bc:.4g}')
        print('bg =', f'{_bg:.4g}')
        print('\nShape Factors:')
        print('sq =', f'{_sq:.4g}')
        print('sc =', f'{_sc:.4g}')
        print('sg =', f'{_sg:.4g}')
        print('\nLoad Inclination Factors:')
        print('iq =', f'{_iq:.4g}')
        print('ic =', f'{_ic:.4g}')
        print('ig =', f'{_ig:.4g}')
        print('\nResults:')
        print('qu =', f'{_qu:.4g}', 'kPa')
        print('R =', f'{_R:.4g}', 'kN')
    return _R