import pytest
from .context import pyec7
from pyec7 import bearing_resistance

def test_basic_bearing_resistance():
    result = bearing_resistance(width=2.0, length=4.0)
    assert isinstance(result, float)
    assert result > 0

def test_strip_foundation():
    result = bearing_resistance(width=2.0, length=0)
    assert isinstance(result, float)
    assert result > 0

def test_zero_friction():
    result = bearing_resistance(width=2.0, length=4.0, friction=0)
    assert isinstance(result, float)
    assert result > 0

def test_zero_cohesion():
    result = bearing_resistance(width=2.0, length=4.0, cohesion=0)
    assert isinstance(result, float)
    assert result > 0

def test_inclined_foundation():
    result = bearing_resistance(width=2.0, length=4.0, inclination=15)
    assert isinstance(result, float)
    assert result > 0

def test_with_embedment():
    result = bearing_resistance(width=2.0, length=4.0, embedment=1.5)
    assert result > bearing_resistance(width=2.0, length=4.0, embedment=0)

def test_with_surcharge():
    result = bearing_resistance(width=2.0, length=4.0, surcharge=10)
    assert result > bearing_resistance(width=2.0, length=4.0, surcharge=0)

def test_with_horizontal_loads():
    result = bearing_resistance(
        width=2.0, length=4.0,
        vertical=100,
        horizontal=10,
        longitudinal_force=5
    )
    assert isinstance(result, float)
    assert result > 0

def test_with_moments():
    result = bearing_resistance(
        width=2.0, length=4.0,
        vertical=100,
        moment=20,
        longitudinal_moment=10
    )
    assert isinstance(result, float)
    assert result > 0

@pytest.mark.parametrize("width,length", [(0, 4.0),(-1, 4.0),(2.0, -1),])
def test_invalid_dimensions(width, length):
    with pytest.raises(Exception):
        bearing_resistance(width=width, length=length)

def test_logging_output(capsys):
    bearing_resistance(width=2.0, length=4.0, log=True)
    captured = capsys.readouterr()
    assert "Geometry" in captured.out
    assert "Soil Properties" in captured.out
    assert "Results" in captured.out
