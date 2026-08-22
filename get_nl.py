#!/usr/bin/env python3
"""Print the number of fine-mode L points (used by pipeline.sh)."""
from config_tsmcN40 import get_config

print(len(get_config('tt', coarse=False)['LENGTH']))
