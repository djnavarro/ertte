---
name: write-news-entries
description: Guidance for writing entries in NEWS.md for the ertte R package. Use whenever adding, updating, or reviewing a NEWS.md entry, e.g. after implementing a new feature, fixing a bug, or making an API change.
---

# Writing NEWS.md Entries

`NEWS.md` is a user-facing changelog, not a commit log or a design document.
Unreleased changes accumulate under the top heading of the file, which
carries the current in-development version number matching `DESCRIPTION`
(e.g. `# ertte 0.0.0.9000`), grouped into `## New features`,
`## Improvements`, `## Bug fixes`, `## Testing`, and `## Documentation`
subsections. When a release is cut, that heading's version number drops the
`.9000` development suffix (e.g. `# ertte 0.1.0`) and a fresh heading with
the next in-development version is added above it for the next cycle.

ertte does not have a `NEWS.md` file yet, and is GitHub-only/pre-CRAN (see
`AGENTS.md`) -- no version of it has ever been released. That means every
change made so far is, by construction, unreleased, and the "same
development cycle" carve-out in rule 3 below doesn't currently have anything
to bite on. Start the file (at the package root, alongside `DESCRIPTION`)
the first time a change genuinely warrants a changelog entry, with a single
`# ertte 0.0.0.9000` heading matching `DESCRIPTION`'s current `Version`
field. Sibling package [erglm](https://github.com/djnavarro/erglm)'s
`NEWS.md` shows the pattern this eventually converges on once CRAN releases
begin: an "Initial CRAN release." heading (e.g. `# erglm 0.1.0`) listing
everything shipped in that first release, and subsequent version headings
(`# erglm 0.1.1`) for what changed since.

Two failure modes recur when writing entries here: restating information
that already lives elsewhere (function documentation, a linked GitHub
issue), and describing changes to code that was never in a CRAN release in
the first place. The rules below exist to avoid both.

## Rules

1. **Keep it short.** One sentence is usually enough; two only if the change
   genuinely needs it. A bare function name like `` ertte_scm_forward() `` is
   auto-linked by pkgdown straight to its help page, so do not restate
   parameter names, defaults, algorithmic detail, or examples that are already
   covered in that function's `@details`/`@examples`. The NEWS entry's job is
   to say *that* something changed and give just enough context to decide
   whether to click through — not to duplicate the documentation.

2. **NEWS.md is strictly user-facing.** Never reference or link to
   agent-facing or contributor-facing material: skills, `AGENTS.md`,
   `.agents/HISTORY.md`/`.agents/PLAN.md`, internal dot-prefixed helper
   functions (e.g. `.ertte_refit()`, `.ertte_simulate_draws()`), CI
   configuration, or "how we implemented this" narrative. If a change has no
   visible effect on the public API, documented behavior, or output (e.g. an
   internal helper was refactored, `.ertte_check_dist()` gained a case
   nothing exported surfaces differently), it does not belong in NEWS.md at
   all — skip it rather than finding a way to phrase it.

3. **Only describe what changed since the last CRAN release.** Check the
   heading at the top of the file: everything under the current
   `.9000`-suffixed development heading (e.g. `# ertte 0.0.0.9000`) is
   unreleased. If a bug being fixed was introduced by a feature added
   earlier in the *same* development cycle (i.e. that feature has never
   shipped to CRAN), it is not a user-visible "bug fix" — it's the feature
   working correctly. Don't add a `## Bug fixes` entry for it; either revise
   the original feature's own bullet if it needs correcting, or just fix the
   code silently. Only regressions in code that already shipped in a
   previous numbered release deserve their own `## Bug fixes` bullet. Since
   ertte has never had a CRAN (or any tagged) release, this distinction
   currently never triggers — every existing change is same-cycle by
   definition — but keep applying the rule once a first release ships.

4. **Point to the issue/PR instead of re-explaining it.** Append `(#N)` to the
   end of the bullet when a GitHub issue or PR number exists. Don't restate
   the issue's background, discussion, or design rationale — that history is
   already written down in the issue itself, one click away.

5. **Match the existing structure and voice.** Group each bullet under the
   closest matching heading (`## New features`, `## Improvements`,
   `## Bug fixes`, `## Testing`, `## Documentation`), only adding a new
   heading if none fits. Start each bullet with a past-tense verb ("Added",
   "Fixed", "Changed", "Removed", "Deprecated").

## Example

Too long — restates documentation detail and mentions an internal helper:

```md
- Added a `criterion` argument to `ertte_scm_forward()` and
  `ertte_scm_backward()` that selects the significance test used at each
  step: `"p-value"` (the existing LRT-based test using `threshold`),
  `"aic"`, or `"bic"`. Internally, `.ertte_refit()` dispatches the refit to
  the matching engine constructor rather than using `stats::update()`, since
  `survreg`/`coxph` objects carry a `$call` that doesn't resolve in a
  caller's frame. `ertte_scm_history()`'s returned tibble now also records
  which criterion drove each step, in a new `criterion` column (#12).

- `.ertte_refit()` now accepts an `ertte_coxph` object as well as
  `ertte_aft`, dispatching on class (#12).
```

Right level of detail — says what changed, links to the function for the
rest, skips the internal-only change entirely:

```md
- Adds a `criterion` argument to `ertte_scm_forward()`/`ertte_scm_backward()`
  that selects between `"p-value"` (the default, LRT-based), `"aic"`, and
  `"bic"` significance testing; `ertte_scm_history()` records which
  criterion drove each step (#12).
```

## Checklist before adding an entry

- [ ] Is this visible to a user of the package, not just to future
      contributors? If not, skip NEWS.md entirely.
- [ ] If this is a bug fix, did the bug exist in a previously released
      version, rather than being a same-cycle regression in unreleased code?
      (Currently always "same-cycle", since ertte has no release yet.)
- [ ] Is the bullet one or two sentences, with parameter/implementation
      detail left to the function's own documentation?
- [ ] Does it rely on pkgdown's auto-linking of bare function names rather
      than re-describing what that link leads to?
- [ ] Does it append `(#N)` instead of re-explaining the linked issue/PR?
- [ ] Does it avoid mentioning skills, `AGENTS.md`, `.agents/*.md`, internal
      dot-prefixed helpers, or other agent/contributor-facing material?
