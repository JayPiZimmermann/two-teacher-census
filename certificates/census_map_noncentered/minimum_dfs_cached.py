"""Bounded exact-cache entry point for :mod:`minimum_dfs`.

This wrapper changes no interval formula, search order, split, verdict, spec,
or certificate format.  It only remembers the value of ``noncentered.atoms``
for an exactly identical interval at the current interval precision.  The
returned ``Atoms`` records are immutable after construction in the DFS code,
so reusing one is the same computation with duplicate work removed.

The cache is deliberately bounded.  A hard singleton can visit millions of
nodes, and an unbounded cache would trade the DFS node budget for an implicit
memory budget.  Final artifacts remain ordinary ``minimum_dfs`` forests and
must be replayed by the unchanged, uncached ``minimum_dfs_check.py``.

Usage is identical to ``minimum_dfs.py``::

  python3 minimum_dfs_cached.py b0 b1 y0 y1 delta dmax tag \
      n_parts lo_idx hi_idx [budget] [minw] [absolute|relative] [scope]
"""
import minimum_dfs as D


CACHE_LIMIT = 4096
_ORIGINAL_ATOMS = D.J.atoms
_ATOM_CACHE = {}
_HITS = 0
_MISSES = 0
_EVICTIONS = 0


def cached_atoms(interval):
    """Return exactly ``noncentered.atoms(interval)``, with FIFO memoization."""
    global _HITS, _MISSES, _EVICTIONS
    # The endpoint tuples are mpmath's exact internal directed endpoints.
    # Precision is part of the key so a caller changing ``iv.prec`` cannot
    # reuse a value rounded under a different interval context.
    key = (D.mpmath.iv.prec, interval._mpi_)
    try:
        result = _ATOM_CACHE[key]
    except KeyError:
        _MISSES += 1
        result = _ORIGINAL_ATOMS(interval)
        if len(_ATOM_CACHE) >= CACHE_LIMIT:
            _ATOM_CACHE.pop(next(iter(_ATOM_CACHE)))
            _EVICTIONS += 1
        _ATOM_CACHE[key] = result
        return result
    _HITS += 1
    return result


def install_atom_cache():
    """Install the exact cache once in this generator process."""
    if D.J.atoms is not cached_atoms:
        if D.J.atoms is not _ORIGINAL_ATOMS:
            raise RuntimeError("noncentered.atoms was already replaced")
        D.J.atoms = cached_atoms


def cache_info():
    """Small diagnostic record used by the regression test."""
    return {"limit": CACHE_LIMIT, "size": len(_ATOM_CACHE),
            "hits": _HITS, "misses": _MISSES,
            "evictions": _EVICTIONS}


def main():
    install_atom_cache()
    return D.main()


if __name__ == "__main__":
    raise SystemExit(main())
