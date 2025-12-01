# Sphere target strength calculator

[![Static badge](https://img.shields.io/pypi/v/spherets.svg)](https://pypi.org/project/spherets/)
[![python](https://img.shields.io/pypi/pyversions/spherets.svg?logo=python&logoColor=white)](https://pypi.org/project/spherets/)
[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/gavinmacaulay/spherets/publish-to-pypi.yml)](https://github.com/gavinmacaulay/spherets/actions/workflows/publish-to-pypi.yml)

The sphereTS package calculates the acoustic target strength (TS) of elastic spheres immersed in a fluid. It is intended for calculating the TS of the spheres used for calibrating echosounders used for quantitative estimates of marine organism backscatter.

Install from pypi:

```bash
pip install spherets
```

and run the GUI from a terminal (type `spherets`) to get the main window:

![main](docs/main%20screen.png "Main screen")

TS values are shown in a separate window:

![results](docs/plot%20screen.png "Results window")

The underlying target strength functions are also available, for example:

```py
import matplotlib.pyplot as plt
from sphereTS import sphere_ts as sts

wc = sts.material_properties()['Tungsten carbide']

c = 1470  # [m/s] sound speed in water
rho = 1027  # [kg/m^3] density of water
a = 0.0381/2  # [m] radius of sphere

ts = sts.sphere_ts(38e3, a, c, wc['c1'], wc['c2'], rho, wc['rho1'])

print(ts)

f, tsf = sts.freq_response(12e3, 200e3, a, c, 
                           wc['c1'], wc['c2'], rho, wc['rho1'])

plt.plot(f/1e3, tsf)
```
