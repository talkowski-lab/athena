from pkg_resources import get_distribution
from athena.cli import master

__version__ = get_distribution('athena').version
