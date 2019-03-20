from pkg_resources import get_distribution
import athena.cli.master

__version__ = get_distribution('athena').version
