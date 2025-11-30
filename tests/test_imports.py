import pytest

def test_import_spherets():
    """Check that importing works."""
    from sphereTS import sphere_ts  # noqa: F401

def test_calculation_resultsdata(tmp_path):
    """Check a few TS results."""

    from sphereTS import sphere_ts
    
    ts = sphere_ts.sphere_ts(38e3, 0.0381, 1490.3, 6853.0, 4171.0, 1027.1, 14900.0)
    assert ts == pytest.approx(-33.896379, abs=1e-4)
    