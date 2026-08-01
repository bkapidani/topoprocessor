from __future__ import annotations
from netgen.libngpy._meshing import _Redraw
from pathlib import Path
from pyngcore.pyngcore import Timer
import sys
from . import config
from . import libngpy
from . import version
__all__: list[str] = ['Path', 'Redraw', 'TimeFunction', 'Timer', 'config', 'libngpy', 'load_occ_libs', 'v', 'version']
def Redraw(*args, **kwargs):
    ...
def TimeFunction(func, name = None):
    ...
def _check_python_version():
    ...
def _get_diagnostics():
    ...
def load_occ_libs():
    ...
__diagnostics_template: str = '\nNetgen diagnostics:\n    sys.platform:          {sys.platform}\n    sys.executable:        {sys.executable}\n    sys.version:           {sys.version}\n    Netgen python version: {config.PYTHON_VERSION}\n    Netgen path            {__file__}\n    Netgen config          {config.__file__}\n    Netgen version         {config.NETGEN_VERSION}\n    sys.path: {sys.path}\n'
__package_name__: str = 'netgen-mesher'
__version__: str = '6.2.2606'
_netgen_bin_dir: str = 'C:\\Windows\\Temp\\tmpi7p95r_f\\wheel\\platlib\\netgen'
_netgen_lib_dir: str = 'C:\\Windows\\Temp\\tmpi7p95r_f\\wheel\\platlib\\netgen'
v: sys.version_info  # value = sys.version_info(major=3, minor=12, micro=9, releaselevel='final', serial=0)
