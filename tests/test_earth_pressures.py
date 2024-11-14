import pytest
import numpy as np
from .context import pyec7
from pyec7 import earth_pressures

def test_basic_earth_pressures():
    result = earth_pressures()
    assert isinstance(result, tuple)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

def test_active_pressures():
    result = earth_pressures(active=True)
    assert isinstance(result, tuple)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

def test_zero_friction():
    result = earth_pressures(friction=0)
    assert isinstance(result, tuple)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

def test_zero_cohesion():
    result = earth_pressures(cohesion=0)
    assert isinstance(result, tuple)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

def test_with_surcharge():
    result = earth_pressures(surcharge=10)
    assert isinstance(result, tuple)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

def test_with_inclination():
    result = earth_pressures(inclination=15)
    assert isinstance(result, tuple)
    assert len(result) == 4
    assert all(isinstance(x, float) for x in result)

@pytest.mark.parametrize("friction", [-1, 91])
def test_invalid_friction(friction):
    with pytest.raises(ValueError):
        earth_pressures(friction=friction)

@pytest.mark.parametrize("cohesion", [-1, -10])
def test_invalid_cohesion(cohesion):
    with pytest.raises(ValueError):
        earth_pressures(cohesion=cohesion)

@pytest.mark.parametrize("delta", [-1, 91])
def test_invalid_delta(delta):
    with pytest.raises(ValueError):
        earth_pressures(delta=delta)

@pytest.mark.parametrize("beta", [-91, 91])
def test_invalid_beta(beta):
    with pytest.raises(ValueError):
        earth_pressures(beta=beta)

@pytest.mark.parametrize("theta", [-91, 91])
def test_invalid_theta(theta):
    with pytest.raises(ValueError):
        earth_pressures(theta=theta)

@pytest.mark.parametrize("surcharge", [-1, -10])
def test_invalid_surcharge(surcharge):
    with pytest.raises(ValueError):
        earth_pressures(surcharge=surcharge)

@pytest.mark.parametrize("inclination", [-91, 91])
def test_invalid_inclination(inclination):
    with pytest.raises(ValueError):
        earth_pressures(inclination=inclination)

def test_logging_output(capsys):
    earth_pressures(log=True)
    captured = capsys.readouterr()
    assert "Input parameters:" in captured.out
    assert "Boundary conditions:" in captured.out
    assert "Coefficients:" in captured.out
