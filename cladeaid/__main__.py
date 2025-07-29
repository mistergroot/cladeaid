# MIT License
# 
# Copyright (c) 2025 Nicolas Locatelli
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import argparse
from . import core_lca
from . import refinery

def main():
    parser = argparse.ArgumentParser(
        description="LCA assignment and downstream refinement tools"
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # CladeAid subcommand
    lca_parser = subparsers.add_parser("lca", help="Run LCA assignment")
    core_lca.add_arguments(lca_parser)

    # Refinery subcommand
    refine_parser = subparsers.add_parser("refinery", 
                                          help="Run postprocessing/refinement")
    refinery.add_arguments(refine_parser)

    args = parser.parse_args()

    if args.command == "lca":
        core_lca.run(args)
    elif args.command == "refinery":
        refinery.run(args)

if __name__ == "__main__":
    main()
