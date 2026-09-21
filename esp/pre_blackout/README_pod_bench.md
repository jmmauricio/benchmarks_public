# A marginally-damped Iberian benchmark for POD controller development

`esp_preblackout_pod_bench.hjson`

## What this case is, and what it is not

**It is** a test system for developing and evaluating power oscillation damping
(POD) controllers on converter-interfaced plant. It is built on a representation of
the Spanish transmission system, and it is deliberately parameterised so that one
inter-area mode is *stable but lightly damped* — a controller has something to act
on, and the system does not fall over without one.

**It is not** a representation of any actual operating point of the Spanish power
system, at any time, and no result obtained with it should be read as a statement
about the real system's stability margins.

**It is not** a statement about any system operator's practice. In particular, the
light damping in this case is **not** produced by disabling or de-tuning power
system stabilisers. Every synchronous unit committed in this case carries a
stabiliser tuned by the residue method (Section "Stabiliser tuning" below), and the
sensitivity table shows what happens if they are switched off — which is that the
system becomes violently unstable, not marginally damped. The margin in this
benchmark comes from the **dispatch**.

## Where the light damping comes from

Three ingredients, none of which is a statement about anyone's engineering.

**1. A forward-looking converter-penetration scenario.** Synchronous plant is
committed at a fraction of the base case, with photovoltaic generation taking up the
released energy. The synchronous machines in this model are plant aggregates, so
reducing a machine's rated power together with its scheduled output represents
*fewer units of that plant committed*; the per-unit machine data, and therefore the
physics of each running unit, are untouched.

**2. An N-1 outage of one ES–FR interconnection circuit.** The interconnection is
two circuits; the case runs with one out. This is a credible contingency that every
operator plans for, and it acts directly on the mode of interest — it raises the tie
reactance and lowers both the frequency and the damping of the Iberia-versus-external
mode. The scheduled exchange is reduced accordingly, as it would be in operation:
with one circuit out the remaining one has a transfer limit of about 2,650 MW, and
the pre-contingency exchange would not be secure against it.

**3. Stabilisers tuned to the applicable grid-code criterion, not beyond it.** The
NTS-SEPE criterion for small-signal performance is that no existing mode falls below
**5%** damping. The stabilisers here are tuned to exactly that, which is what a
compliant system looks like. Tuning them harder is possible — against a 10% target
the residue method damps the inter-area mode to about 14% in this system — but a
benchmark tuned to 10% is tuned beyond the requirement, and the margin it shows says
more about the target chosen than about the system.

This is the mechanism the benchmark rests on, and it is a property of
converter-dominated operation rather than of anyone's engineering:

- As synchronous plant is displaced, the stabilisers that damp inter-area
  oscillations are displaced with it. There are simply fewer of them, on less
  capacity.
- Synchronous machines have weaker authority over inter-area modes than over their
  own local modes — in this system, residues roughly five times smaller — so as
  committed capacity falls, inter-area damping is the first thing stabilisers lose
  the ability to hold.

## What this case does and does not demonstrate

Worth stating precisely, because the obvious stronger claim is not true here.

It is **not** the case that synchronous machines cannot damp this mode. Re-tuning
the stabilisers by the residue method against a 10% target damps it to about 14%,
at a range of committed capacities, with and without the interconnector outage.
Given a target, the method reaches it. A benchmark that looked marginal because its
stabilisers had been left tuned for a *different* dispatch would be an artefact, and
this case is deliberately not that.

What the case does show is that with stabilisers tuned to the **grid-code criterion
rather than beyond it**, the system settles at that criterion — around 5% — and the
inter-area mode is the one sitting on it. The damping in this benchmark is set by the
tuning target, honestly and by construction, not by a claimed limitation of
synchronous plant.

That is the right framing for a POD study. The controller's job is not to rescue a
system that stabilisers have failed; it is to provide margin **above the regulatory
minimum** under a stressed, converter-dominated, N-1 condition — where the
alternative would be tuning conventional stabilisers past what the code requires, on
plant that is being displaced anyway.

## Stabiliser tuning

All stabilisers in this case are tuned with
`pydae.bps.utils.pss_tuner.tune_psss`, which implements the residue method: for the
worst-damped electromechanical mode it selects the machine with the largest residue
on the loop the stabiliser closes (`v_ref` → `omega`), fits the lead-lag network
across the modes that machine has authority over, and chooses the gain that
maximises the *worst* damping in the system rather than that of the targeted mode.

Machines whose stabilisers are not committed have their `pss` block commented out
rather than left at zero gain, so the model carries no inert stabiliser states. The
tuner's report records, for every machine it declined, the reason — most often that
the required phase compensation exceeded what a two-stage lead-lag can supply at
that frequency, or that the machine had no meaningful authority over the mode.

## The case as shipped

| | |
|---|---|
| Committed synchronous capacity | 52% of the base case, by stored energy |
| Iberian stored kinetic energy | 55,512 MW·s (46% of the 119,474 MW·s reported at 12:30 on 28 April 2025) |
| Synchronous share of generation | 15.9% — 4,411 MW against 19,970 MW PV and 3,398 MW wind |
| ES–FR interconnection | one circuit of two in service (N-1) |
| ES–FR exchange | 988 MW export |
| Stabilisers | tuned by the residue method to the 5% criterion; 2 committed, the rest commented out |

Resulting electromechanical modes:

| mode | frequency | damping |
|---|---|---|
| **Inter-area, Iberia vs external** | **0.262 Hz** | **5.14%** |
| Worst mode in the system | 1.113 Hz | 4.24% |

Nothing is unstable, and the inter-area mode is the target for a damping controller.

### Sensitivity to the stabilisers

| stabiliser state | worst mode | inter-area mode |
|---|---|---|
| As shipped, tuned to the 5% criterion | 4.24% | 0.262 Hz, **5.14%** |
| Removed | **−16.22%** | 0.296 Hz, −16.22% |

This table is the point. The benchmark's modest margin is **not** obtained by
removing or weakening stabilisers — without them the same case is unstable by 16%.
It is obtained by a stressed operating condition with stabilisers tuned to, and
satisfying, the applicable criterion.

### Controllability of the target mode

Largest residues on the 0.262 Hz mode from synchronous plant, on the loop a
stabiliser closes: `hydro_ES432_hydro` 0.0227, `hydro_ES113_hydro` 0.0143. These are
an order of magnitude below the residues the same machines have on their local modes,
which is why holding this mode above the criterion costs a gain of ~15 on one unit
and why converter-interfaced plant, distributed and electrically closer to the load,
is worth investigating for the duty.

## Provenance of the data

| Quantity | Source |
|---|---|
| Network, machine, AVR and governor data | `esp_preblackout_dev.hjson` in this directory |
| External system | Single area equivalent behind the ES–FR interconnection |
| Dispatch | Constructed for this benchmark; see above |

Note on the external equivalent: it is an *area aggregate*, and its rotor flux time
constants are set long deliberately so that its fluxes are effectively constant over
the period of a sub-Hz mode — the classical machine form conventionally used for
area equivalents. Leaving an aggregate with a single turbo-generator's round-rotor
data is not neutral here: a `T1q0` of 0.8 s places the q-axis corner at 1.25 rad/s,
directly on top of the modes of interest, and that one parameter alone contributed
about twelve points of damping to a 0.2 Hz mode in testing.

## Suggested wording for a paper

> The test system is a dynamic model of the Spanish transmission network comprising
> N buses and M synchronous units, with photovoltaic and wind generation represented
> by aggregated converter models. To obtain a test case in which inter-area damping
> is low enough to evaluate a damping controller, the system is dispatched in a
> forward-looking high-converter-penetration scenario in which synchronous plant is
> committed at a reduced level and photovoltaic generation supplies the balance.
> This dispatch is constructed for the purposes of this study and is not intended to
> represent any observed or planned operating point.
>
> The scenario additionally considers the outage of one of the two ES–FR
> interconnection circuits, with the scheduled exchange reduced to approximately
> 1,000 MW as the remaining circuit's transfer capability requires.
>
> All synchronous units committed in this scenario are equipped with power system
> stabilisers tuned by the residue method against the small-signal criterion
> applicable in this system, namely that no mode falls below 5% damping. Under this
> tuning every electromechanical mode satisfies the criterion, with the inter-area
> mode at 0.26 Hz retaining a damping ratio of 5.1% — that is, the system operates
> at the minimum required margin rather than above it. With the stabilisers removed
> the same case is unstable at −16%, confirming that the low margin is a property of
> the operating condition and the tuning target, not of any absence of conventional
> damping control. The case therefore represents a stressed condition in which
> conventional stabilisers, correctly tuned to the applicable criterion, leave no
> margin above it, motivating the use of converter-interfaced plant to provide
> additional damping.

The essential points to preserve if this is reworded: the scenario is *constructed*,
forward-looking and contingency-stressed; the stabilisers are *tuned to the
applicable criterion*, not disabled or de-tuned; and the modest margin follows from
operating at that criterion rather than from any claimed shortcoming of synchronous
plant. Avoid any phrasing suggesting stabilisers are missing, off, or badly set.

## Reproducing

```
python main.py bench
```

Reports the target mode, the stabiliser sensitivity table, and a step response.
