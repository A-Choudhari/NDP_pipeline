#!/usr/bin/env bash
# trust_check.sh -- refuse to run shared code that anyone but its owner could
# have changed.
#
#   source trust_check.sh
#   trust_check DIR [DIR ...] || exit 2
#
# WHY. Everyone in the lab runs one shared copy of this toolkit, and its owner
# runs it unattended every night. Whoever can change that code decides what
# runs under those accounts. Three ways of changing it were found here, each
# closed by hand, none of them visible to anything that ran the code:
#
#   2026-09-10  the toolkit directories were group-writable;
#   2026-09-16  the folder ABOVE them was group-writable, so any group member
#               could rename the pipeline folder and put a directory of their own
#               in its place, which the 03:00 job would then have run;
#   2026-09-16  eleven compiled caches (lib/__pycache__/*.cpython-39.pyc)
#               belonged to another account. A file's owner can rewrite it in
#               any directory, and Python 3.9 runs a cache in place of its
#               source whenever the cache's header matches that source.
#
# So the check runs every time, before any code is executed.
#
# WHAT IT CHECKS, for each DIR, after resolving symbolic links:
#
#   1. Every file and directory beneath DIR belongs to DIR's owner and is not
#      writable by its group or by others. `.git` is skipped. Symbolic links
#      are skipped too: their own mode is always 777 and means nothing, and
#      whatever they point at is checked where it lives when that directory is
#      passed in as well.
#   2. Every directory above DIR, up to /, belongs to root or to DIR's owner,
#      and does not let its group or others rename or replace DIR -- unless its
#      sticky bit is set, which forbids exactly that (as on /tmp).
#
# Prints one line per problem and returns 1; returns 0 when there are none.
# It reads permission bits. Filesystem ACL entries that grant write access
# without showing in those bits are not seen.

trust_check() {
    local d real owner p pown line found=0
    for d in "$@"; do
        if ! real="$(cd "$d" 2>/dev/null && pwd -P)"; then
            echo "  cannot resolve $d"
            found=1
            continue
        fi
        owner="$(stat -c %U "$real")"

        # 1. everything beneath it
        while IFS= read -r line; do
            echo "  $line"
            found=1
        done < <(find "$real" -path '*/.git' -prune -o ! -type l \
                    \( ! -user "$owner" -o -perm -g=w -o -perm -o=w \) \
                    -printf '%M %u  %p\n' 2>/dev/null | head -n 20)

        # 2. everything above it
        p="$real"
        while [[ "$p" != "/" ]]; do
            p="$(dirname "$p")"
            pown="$(stat -c %U "$p")"
            if [[ "$pown" != "root" && "$pown" != "$owner" ]]; then
                echo "  $(stat -c '%A %U' "$p")  $p  (owned by another account, which could replace $real)"
                found=1
            fi
            if [[ -n "$(find "$p" -maxdepth 0 \( -perm -g=w -o -perm -o=w \) \
                        ! -perm -1000 -print 2>/dev/null)" ]]; then
                echo "  $(stat -c '%A %U' "$p")  $p  (its group or others could rename or replace $real)"
                found=1
            fi
        done
    done
    return "$found"
}

# Run directly, it checks the directories given, or its own.
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    set -u
    if [[ $# -eq 0 ]]; then
        set -- "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    fi
    if trust_check "$@"; then
        echo "  OK: only the owner can change $*"
    else
        exit 1
    fi
fi
