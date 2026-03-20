"""
Tests for ospgrillage.ospgui.

The GUI module requires PyQt6, which is an optional dependency.  All tests
here are designed to run in headless / CI environments where PyQt6 is absent.
They verify that:

  1. The module imports cleanly even when PyQt6 is not installed.
  2. main() reports the missing dependency and exits with code 1.
"""

import importlib
import sys
import os

import pytest

sys.path.insert(0, os.path.abspath("."))
sys.path.insert(0, os.path.abspath("../"))


def test_ospgui_importable_without_pyqt6():
    """ospgui must import cleanly even when PyQt6 is not installed."""
    mod = importlib.import_module("ospgrillage.ospgui")
    assert hasattr(mod, "_PYQT6_AVAILABLE")
    # When PyQt6 happens to be installed the flag is True – nothing to
    # assert about absence, but the module still imported cleanly.
    if mod._PYQT6_AVAILABLE:
        pytest.skip("PyQt6 is installed; cannot verify the absent-flag path")


def test_ospgui_main_exits_when_pyqt6_absent(capsys):
    """main() must call sys.exit(1) and write a helpful message to stderr."""
    from ospgrillage.ospgui import main, _PYQT6_AVAILABLE

    if _PYQT6_AVAILABLE:
        pytest.skip("PyQt6 is installed; this test only applies when it is absent")

    with pytest.raises(SystemExit) as exc_info:
        main()

    assert exc_info.value.code == 1
    captured = capsys.readouterr()
    assert "pip install" in captured.err
