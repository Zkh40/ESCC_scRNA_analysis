# Compatibility source attribution

R/compat_monocle2.R adapts the Monocle2 DDRTree projection and ordering implementation. Monocle is developed by Cole Trapnell and contributors; upstream source: https://github.com/cole-trapnell-lab/monocle-release. The installed version used by this analysis is 2.34.0. These adaptations retain the upstream license applicable to Monocle; they are not represented as wholly original implementation. Changes replace deprecated igraph traversal access, use deterministic dense Prim MST construction for projected distances, and resolve helper functions in a local closure rather than mutating the package namespace.

The project analysis wrapper and dense Prim implementation were recovered from the frozen stage1 script. The release must include the upstream license text and attribution for adapted portions.
