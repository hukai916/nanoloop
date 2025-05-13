from .cli import main

__version__ = None

# Use version from pyproject.toml
with open('pyproject.toml', 'r') as f:
  for line in f:
    if 'version' in line:
      __version__ = line.split('=')[1].strip().strip('"')
      break
