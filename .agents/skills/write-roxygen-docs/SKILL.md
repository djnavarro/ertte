---
name: write-roxygen-docs
description: Guidance for writing and reviewing roxygen2 documentation comments in the ertte R package. Use whenever adding a new exported function, editing an existing @param/@returns/@details/@examples block, or reviewing a roxygen comment before running devtools::document().
---

# Writing roxygen2 Documentation

Roxygen comments become the content of `?function` and the pkgdown reference
site — they are read by a human user deciding whether and how to call a
function, not by a future contributor or an agent. Getting the mechanics
right (tags, `@export`, blank lines) is the easy part; agents reliably get
that right already. The failure modes worth guarding against are about
*content*: putting the wrong thing in a section, writing at the wrong level
of detail, or leaking information that shouldn't be there at all.

This skill assumes familiarity with roxygen2 basics. For deeper background on
any topic below, see [R Packages (2e), ch. 16](https://r-pkgs.org/man.html).

## What goes where

Each part of the introduction has a distinct job. Don't let content drift
into the wrong one:

- **Title** (first sentence, sentence case, no full stop): what the function
  does, distinguishing it from sibling functions. ertte's two engines,
  `ertte_aft()` and `ertte_coxph()`, are distinct constructors (not one
  constructor with an `engine =` switch), so their titles should each name
  the model family *and* the underlying `survival::` function they wrap
  (e.g. "Fit an exposure-response time-to-event AFT model based on
  `survreg()`" vs the Cox PH equivalent based on `coxph()`), not a generic
  verb that could describe either. The same applies to shared-name
  generic/method pairs like `ertte_predict()`/`ertte_predict.ertte_aft()`/
  `ertte_predict.ertte_coxph()`: the generic's title describes the shared
  contract, while each method's title (if documented separately) should
  name what's distinctive about that engine's computation.
- **Description** (next paragraph): one paragraph on *this* function's
  purpose, in different words than the title — not a restatement of it, and
  not a template copied from a sibling function. It's easy to start the
  Cox PH sibling of an existing AFT function (or vice versa) by copying its
  docs and forget to re-target the description; the tell is a description
  that reads correctly for *either* engine but never actually says what
  makes this one distinct (e.g. closed-form CDF evaluation vs
  `survival::survfit()`/baseline-hazard machinery, or an analytic
  delta-method SE vs an exact step-function sum) — that's a sign the real
  description got left in `@details` instead. It's just as easy to skip the
  description paragraph entirely — roxygen2 doesn't warn you, it silently
  reuses the title as the description, which is *always* a restatement by
  construction.
- **`@details`**: everything else — default behaviour, edge cases, how an
  argument being `NULL` is treated (e.g. `censor_time = NULL` in
  `simulate.ertte_model()` falling back to capping censored rows at their
  observed exit time while leaving event rows uncensored; `ref = NULL` in
  `ertte_power()` defaulting to `median(x, na.rm = TRUE)`), interactions with
  other parts of the fitting/SCM pipeline (e.g. how `criterion` changes
  which significance test `ertte_scm_forward()`/`ertte_scm_backward()` use).
  It's fine for this to be a few sentences to a short paragraph. Details
  render *after* arguments and return value on the help page, so don't put
  anything here that a reader needs before they can parse `@param`.
- **`@param`**: a succinct summary of what the argument controls and, if it
  has a fixed set of values (like `dist` or `criterion`), what they are.
  State the default inline (e.g. "Defaults to `\"weibull\"`") since the
  usage block and the argument description are far apart on the rendered
  page. When the default is a sentinel like `NULL` whose *effect* isn't
  obvious from the value itself, say what it does rather than just naming it
  (e.g. "If `NULL` (the default), one is chosen automatically and reported
  via a message" reads better than "the default is `seed = NULL`", which
  tells the reader nothing until they go read `@details` too).
- **`@returns`**: the shape of the return value — for the fitting functions,
  "A survreg object with extra `ertte_aft`/`ertte_model` classes" or the
  `coxph`/`ertte_coxph` equivalent; for `ertte_predict()`/`ertte_rmst()`/
  `ertte_landmark()`, the tibble's row/column structure; for
  `ertte_scm_history()`, the class of the returned history object. Every
  exported function must have this tag.
- **`@examples`**: runnable code showing typical usage. Not a place to
  re-explain arguments already covered in `@param`. Fitting on
  `ertte_aft()`/`ertte_coxph()` against the bundled `ertte_data` is fast, so
  (unlike packages with an iterative-optimiser fitting step) there's no
  options object or time limit to pass just to keep examples from hanging —
  a plain `ertte_aft(Surv(time, event) ~ aucss, ertte_data)` is a complete,
  fast example.

## Calibrating detail

Match documentation density to how novel the content actually is:

- Across a shared-name generic/method pair — e.g.
  `ertte_predict.ertte_aft()`/`ertte_predict.ertte_coxph()`, or
  `ertte_fun.ertte_aft()`/`ertte_fun.ertte_coxph()` — the argument shape is
  identical (`object`, `newdata`, `time`/`param`, `conf_level`). That
  similarity belongs in the generic's own title/description, not restated
  at length in every method's docs. Spend the words on what's actually
  distinctive: closed-form base-distribution CDF vs `survfit()`'s baseline
  hazard + linear predictor and log-transform CI, or `extend = TRUE`'s
  effect on predictions beyond the last observed follow-up time.
- A dense wall of text is harder to scan than the same information broken
  into a sentence or two per idea, or a short bullet list (roxygen2 markdown
  supports `* item` lists in any prose section). Prefer that when an
  argument has more than two or three possible values or behaviours (e.g.
  the four `dist` choices for `ertte_aft()`, or the three `criterion`
  choices — `"p-value"`, `"aic"`, `"bic"` — for the SCM functions).
- Don't pad a short, genuinely simple function's documentation just to make
  it look thorough. `ertte_power()` or `ertte_landmark()` are compact by
  design (the latter reduces to a single `ertte_predict()` call with the CI
  bounds swapped); if the description already says everything, an empty or
  one-line `@details` — or omitting the tag — is correct.
- Once `@details` covers more than three or four distinct sub-topics (e.g.
  `simulate.ertte_model()` documenting coefficient sampling, the per-engine
  event-time draw method, and the `censor_time` fallback logic together),
  break it into markdown headings or `@section` blocks, one per sub-topic,
  instead of one long unbroken block of paragraphs. A reader looking for one
  specific fact shouldn't have to read the whole section serially to find
  it.
- `@section` titles must be capitalized (R Core's own
  [Rd file guidelines](https://developer.r-project.org/Rds.html) state this
  explicitly for both `\title` and `\section` titles). Don't just reuse a
  lowercase identifier verbatim as a heading — prefer a short, readable
  capitalized phrase, and refer to the actual identifier in the body text
  instead, in backticks.

## Keep it user-facing

Roxygen documentation ships to end users via `?function` and pkgdown. It is
governed by the same boundary as `NEWS.md` (see the `write-news-entries`
skill):

- **Never reference agent- or contributor-facing material.** No mentions of
  skills, `AGENTS.md`, `.agents/HISTORY.md`/`.agents/PLAN.md`, or CI
  configuration.
- **Don't name internal dot-prefixed helper functions or otherwise describe
  implementation details that could change.** Explain behaviour in terms of
  what the function does and what the user observes (inputs accepted,
  outputs produced, what the fitted object contains), not the private
  helper functions or code paths used to get there. This applies even when
  the temptation is to point at *how* a documented result is computed — e.g.
  document that `simulate.ertte_coxph()` draws event times by inverting the
  fitted baseline cumulative hazard, but not that this happens via
  `.ertte_simulate_draws.ertte_coxph()`; document that `ertte_scm_forward()`
  refits candidate models without `stats::update()`, but not that this goes
  through `.ertte_refit()`. Naming a dot-prefixed internal function in
  documentation also invites users to reach for it with `:::`, which is
  best avoided. If you find yourself writing "internally, this calls..." or
  "`.ertte_refit()` is used to...", cut the function name and keep only the
  behavioural consequence.

  Note: several existing roxygen blocks in this package (e.g.
  `simulate.ertte_model()`, `ertte_rmst()`) currently *do* name internal
  helpers like `.ertte_simulate_draws()` or `.ertte_rmst_pfun_delta()`
  in `@details`. That predates this rule and is tracked as cleanup, not a
  pattern to extend — don't use existing docs as a precedent when writing
  new ones.
- **Write for a reader who has never seen the source.** Avoid phrasing that
  only makes sense with the R script open (e.g. "as shown above", "the
  parameter vector described earlier" referring to code, not prose already
  in the same doc).
- **Cross-reference with square brackets, not just backticks.** Writing
  `` `ertte_predict()` `` renders as code but produces no link;
  `[ertte_predict()]` (or `[erplots::er_plot_add_model()]` for another
  package) is what roxygen2/pkgdown turn into an actual hyperlink. A
  backtick-only mention anywhere a function is referenced — in `@details`,
  `@seealso`, or prose — is a broken cross-reference, not a stylistic
  choice. When the natural link text is a function call but the useful
  target is a different topic, use `` [`function()`][topic] `` for custom
  link text rather than linking to a less useful default target.
- **Check the compiled `.Rd` for stray aliases, not just the prose.** A
  `@name`/`@rdname` block (shared docs for several functions/methods)
  attaches to whichever R object immediately follows it in the file. If an
  internal, unrelated object sits between the block and its intended first
  function, roxygen2 silently attaches the block to that object instead and
  gives it a public `\alias`, which can surface as nonsense in the rendered
  `\usage{}`. Reading the roxygen comments won't reveal this; open the
  generated `man/*.Rd` and check that every `\alias{}` is something the
  topic is actually meant to document.

## Checklist before finishing a roxygen block

- [ ] Title names what's distinctive about this function (engine, or
      generic vs specific method); description says the same thing in
      different words, not a restatement of the title.
- [ ] Description actually describes *this* function/method, not a template
      inherited from the AFT/Cox sibling that describes the pair in general
      instead.
- [ ] An explicit `@description` paragraph was actually written, distinct
      from the title — not left to roxygen2's default of silently reusing
      the title verbatim.
- [ ] Everything in `@details` is genuinely additional to the description,
      not filler to make the section non-empty.
- [ ] If `@details` covers more than three or four sub-topics, it's broken
      into headings/`@section`s rather than left as one long block.
- [ ] Every `@section` title is capitalized, and reads as a short phrase
      rather than a bare lowercase identifier copied from the code.
- [ ] Every `@param` states the default where one exists, and enumerates
      fixed value sets (e.g. `dist`, `criterion`).
- [ ] `@returns` is present and names the concrete class/shape returned.
- [ ] Every reference to another function (in `@details`, `@seealso`, or
      prose) uses square brackets so it actually renders as a link, not
      backticks alone.
- [ ] No mention of skills, `AGENTS.md`, `.agents/*.md`, CI, or other
      agent/contributor-facing material.
- [ ] No description of implementation details a user would need `:::` to
      verify, and no named `.ertte_*()` internal helper — only observable
      inputs/outputs/behavior.
- [ ] For `@name`/`@rdname` topics, every `\alias{}` in the compiled
      `man/*.Rd` is something the topic is actually meant to document — no
      internal object picked up by accident.
- [ ] Ran `devtools::document()` and skimmed the rendered `man/*.Rd` (or
      `?function` output) rather than just the roxygen comment source. If
      the fix touched code (not just comments) — e.g. reordering
      definitions to fix an alias leak — also ran `devtools::test()`.
