# Documentation style

## Lie-group notation

In prose that is rendered by Documenter, write Lie-group names as inline
LaTeX using double backticks and upright group letters: ``\mathrm{SU}(N)``,
``\mathrm{SO}(N)``, ``\mathrm{Sp}(N)``, ``\mathrm{U}(1)``, and
``\mathrm{Z}_N``. Use the same form for concrete examples, such as
``\mathrm{SU}(2)``.

In Julia docstrings, double each LaTeX backslash so the Julia string parser
preserves it; for example, write ``\\mathrm{SU}(N)`` in the source.

Keep Julia identifiers, type signatures, code examples, filenames, and cache
keys in code formatting, e.g. `SU{N}`, `SO{N}`, `Sp{N}`, and `U1`. Those are
program syntax rather than mathematical prose.
