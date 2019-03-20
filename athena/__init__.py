from pkg_resources import get_distribution
import athena.cli

__all__ = []

__version__ = get_distribution('athena').version
