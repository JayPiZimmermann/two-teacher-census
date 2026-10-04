"""Fresh-process selected-minimum replay with bounded exact atom memoization.

This entry point installs the cache from :mod:`minimum_dfs_cached` and then
calls :func:`minimum_dfs_check.main` unchanged.  It changes no reader, root,
split, proof method, target-verdict recomputation, coverage check, or failure
condition.  The cache key contains the current interval precision and
mpmath's exact directed endpoint tuple; final artifacts remain ordinary
``minimum_dfs`` forests and can always be replayed by the uncached entry point.

Usage: python3 minimum_dfs_check_cached.py minimum_dfs_<tag>_pNNN.json
"""
import minimum_dfs_cached as C
import minimum_dfs_check as K


def main():
    C.install_atom_cache()
    return K.main()


if __name__ == "__main__":
    raise SystemExit(main())
