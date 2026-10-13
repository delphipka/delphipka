#!/usr/bin/env python3
# SPDX-License-Identifier: LGPL-3.0-or-later
# Modified on branch nucleic-acids by noshou: new script, clang-format wrapper (Linux, macOS, Windows).

"""Format the C/C++ sources with clang-format, using the style in .clang-format (4 spaces, no tabs, Allman braces).

    python scripts/format.py              # format everything under src/ in place
    python scripts/format.py --check      # list the files that are not formatted and exit 1 if there are any; changes nothing
    python scripts/format.py FILE ...     # only these files

clang-format comes with LLVM: `sudo dnf install clang-tools-extra` (Fedora), `sudo apt install clang-format` (Debian, Ubuntu),
`brew install clang-format` (macOS), `pacman -S mingw-w64-x86_64-clang-tools-extra` (MSYS2) or the LLVM installer on Windows.
Only whitespace and the like changes; the tokens of the program stay the same.
"""
import argparse
import os
import shutil
import subprocess
import sys

EXT = (".cpp", ".h", ".hpp", ".c")


def find_clang_format():
    names = ["clang-format"] + ["clang-format-%d" % v for v in range(40, 12, -1)]
    for name in names:
        path = shutil.which(name)
        if path:
            return path
    return None


def format_text(tool, path, text):
    """clang-format settles in one pass on nearly every file; a few constructs (trailing comments on #endif lines) need a second
    or third, so run until the output stops changing. That is what makes --check reliable."""
    cur = text
    for _ in range(4):
        out = subprocess.run([tool, "-style=file", "--assume-filename=" + path], input=cur, stdout=subprocess.PIPE, check=True).stdout
        if out == cur:
            break
        cur = out
    return cur


def sources(root):
    for base, dirs, files in os.walk(os.path.join(root, "src")):
        dirs.sort()
        for f in sorted(files):
            if f.endswith(EXT):
                yield os.path.join(base, f)


def main():
    ap = argparse.ArgumentParser(description="Format the C/C++ sources with clang-format.")
    ap.add_argument("--check", action="store_true", help="only report files that would change")
    ap.add_argument("files", nargs="*", help="files to format (default: everything under src/)")
    a = ap.parse_args()

    tool = find_clang_format()
    if not tool:
        sys.exit("clang-format not found on the PATH (see the top of this script for how to install it)")
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if not os.path.exists(os.path.join(root, ".clang-format")):
        sys.exit("no .clang-format in %s" % root)
    files = a.files or list(sources(root))

    changed = []
    for path in files:
        old = open(path, "rb").read()
        new = format_text(tool, path, old)
        if new != old:
            changed.append(path)
            if not a.check:
                with open(path, "wb") as fh:
                    fh.write(new)
    rel = lambda p: os.path.relpath(p, os.getcwd())
    if a.check:
        for p in changed:
            print("not formatted:", rel(p))
        print("%d of %d files need formatting" % (len(changed), len(files)))
        return 1 if changed else 0
    print("formatted %d of %d files" % (len(changed), len(files)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
