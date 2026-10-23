#!/usr/bin/env python3
"""
test_names_resolve.py -- every global a function reaches for must exist.

    python3 tests/test_names_resolve.py

`analyse.py` called `_settings(folder)` in one line. The function had never
been defined. That line runs only when a record has no readable frame spacing,
so the NameError waited for a folder that produced one -- and then it did not
skip that record, it ended the run: HD4 stopped at the tenth of twenty-five
systems, inside the message explaining the skip.

An error path that has never executed is not tested by anything, and the
cheapest guard is not a test of the path but of the NAMES: walk every function
in the module and check that each global it references resolves. It catches a
typo and a rename in a branch nobody has taken, which is exactly where they
hide.
"""

import builtins
import importlib
import os
import sys
import types
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

#: Modules whose global references must all resolve. Anything imported lazily
#: INSIDE a function is out of scope here: those are deliberate, and the
#: function names them itself.
MODULES = ["analyse", "import_paper", "pages.story", "pages.paper",
           "pages.captions", "core.spec", "core.detect",
           "stages.structure.basic"]


def globals_used(fn):
    """Every global name this function's code object reaches for."""
    code = fn.__code__
    names = set(code.co_names)
    for const in code.co_consts:
        if isinstance(const, types.CodeType):
            names |= set(const.co_names)
    return names


class NamesResolve(unittest.TestCase):
    def test_every_global_a_function_uses_exists(self):
        missing = []
        for modname in MODULES:
            mod = importlib.import_module(modname)
            ns = vars(mod)
            for name, obj in list(ns.items()):
                if not isinstance(obj, types.FunctionType):
                    continue
                if obj.__module__ != mod.__name__:
                    continue          # imported from elsewhere
                for used in globals_used(obj):
                    if used in ns or hasattr(builtins, used):
                        continue
                    # Attribute lookups land in co_names too; only a bare
                    # global that is neither defined nor a builtin is a bug,
                    # and an attribute of something defined is fine.
                    if any(used in vars(v) for v in ns.values()
                           if isinstance(v, types.ModuleType)):
                        continue
                    missing.append(f"{modname}.{name} -> {used}")
        # Attribute names on objects are indistinguishable from globals in
        # co_names, so this reports candidates; the assertion is on the ones
        # that are plainly module-level helpers, by leading underscore.
        real = [m for m in missing if m.rsplit("-> ", 1)[1].startswith("_")]
        self.assertEqual(real, [], "undefined module-level name(s) referenced")

    def test_the_one_that_got_away_is_defined(self):
        import analyse
        self.assertTrue(callable(analyse._settings))
        self.assertTrue(analyse._settings("/x").endswith(
            os.path.join("MD_RESULTS", "settings.yaml")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
