"""Free-Energy Landscape Tool"""
from __future__ import print_function, absolute_import, division
import sys
from ..cmdline import App
from ..commands import *
from ..version import version
# the commands register themselves when they're imported

# Load external commands which register themselves
# with entry point felt.commands
from pkg_resources import iter_entry_points

for ep in iter_entry_points("felt.commands"):
    external_command = ep.load()
    # Some groups start with numbers for ordering
    # Some start with descriptions e.g. "FELT"
    # Let's set the group to start with ZZZ to put plugins last.
    external_command._group = "ZZZ-External_" + external_command._group
    print(ep)


class FELTApp(App):
    pass


def main():
    try:
        app = FELTApp(name='FELT', description=__doc__)
        app.start()
    except RuntimeError as e:
        sys.exit("Error: %s" % e)
    except Exception as e:
        message = """\
An unexpected error has occurred with FELT (version %s)
"""
        print(message % version, file=sys.stderr)
        raise

if __name__ == '__main__':
    main()




