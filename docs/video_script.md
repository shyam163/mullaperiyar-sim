# "I Flooded Kerala on My PC" — 6-minute video script

Narration ≈ 850 words ≈ 5:40 at a relaxed 150 wpm (headroom to 6:00). All visuals reference
real files in the repo. VO = voiceover. Text on screen in [brackets].

---

## 0:00 – 0:38 · COLD OPEN — the scary version

**VISUAL:** Black. One second of silence. Then straight into
`outputs/cascade_142/animation_3d.mp4` (the sea-view 3D film), full screen,
low ominous drone under it. As the flood spreads across the coastal plain,
timestamp overlay ticking.

**VO:**
This is a simulation of the Mullaperiyar dam failing. The gorge fills in
under an hour. By hour nine there's a forty-metre wave at Neriamangalam. By
hour twenty-two, the water is walking into the suburbs of Kochi.

**VISUAL:** Freeze frame at the widest flood extent. Stamp across it:
[THIS IS THE VERSION EVERYONE SHOWS YOU]

**VO:**
It's real physics. I computed it myself. And it is — quietly, crucially —
built on an assumption almost nobody mentions. Stay till the end, because
the honest version of this story is stranger than the scary one.

---

## 0:38 – 0:58 · THE SETUP — why I did this

**VISUAL:** Old survey-map aesthetic (the blog's gauge graphic works —
screen-record the page). Quick cuts: the dam (archival photo), a gauge with
136 / 142 / 152 ft marked, a newspaper headline montage.

**VO:**
Quick background. The Mullaperiyar dam is a hundred and thirty years old,
made of lime-surkhi masonry, sits in Kerala, is operated by Tamil Nadu, on
a lease signed in 1886 — for nine hundred and ninety-nine years. Kerala
says it's a time bomb upstream of three and a half million people. Tamil
Nadu says it's fine. The Supreme Court referees a fight over sixteen
vertical feet of water level. Everything the public knows traces back to
one famous study — so I decided to run the physics myself.

---

## 0:58 – 1:38 · THE BUILD — sixty seconds of method

**VISUAL:** Fast montage: Copernicus tiles downloading (terminal), the
hillshaded DEM (`outputs/*/max_depth.png` background layer), code scrolling
(solver.py), the smoke-test output "ALL SMOKE TESTS PASSED", the mass
ledger line "+0.000%".

**VO:**
The recipe: real terrain, measured from orbit at thirty metres. A
shallow-water solver I wrote from scratch — four hundred lines, every cubic
metre of water accounted for, to zero point zero zero zero percent. And a
breach model from the standard engineering regression — Froehlich
2008, fitted to every historical dam failure we have data for.

**VISUAL:** The canopy problem: zoom into the gorge on the DEM, overlay
[the satellite measured the TREES].

**VO:**
One catch: the satellite measured the treetops. The Periyar's canyon is
covered in forest, so my first flood hit a wall of oaks-as-geology and
refused to move. A day of terrain surgery later, the valley agreed to be a
valley. Remember that — it matters at the end.

---

## 1:38 – 2:38 · THE REALISTIC RUN — and the plot twist

**VISUAL:** `outputs/baseline_142/animation.mp4` (2D film) full screen.
Pause at Vandiperiyar as the gorge fills; cut to the town gauge readout.

**VO:**
Here's the realistic scenario: today's permitted water level, the
standard breach. First — and I want to say this without any irony — the
gorge is not okay. Vandiperiyar, a real town, takes twenty to thirty
metres of water within an hour or two, in every scenario, under every
assumption. For the people in that valley this dispute isn't abstract.
It's how much warning their families would get.

**VISUAL:** The film continues: the surge slides north into the Idukki
arms… and stops. Hold. Cut to `outputs/baseline_142/animation_3d.mp4`,
the "nothing reaches the coast" film. On screen: [+5.2 m … and it holds]

**VO:**
Then, the plot twist. The entire flood — three hundred and eighty million
cubic metres — arrives at the Idukki reservoir. And fits. The great lake
rises about five metres and holds all of it. In the realistic scenario,
nothing north of Idukki gets wet at all. The apocalypse you saw at the
start requires a sequel: a second, much larger, modern dam has to fail
too. That's the assumption nobody mentions.

---

## 2:38 – 3:13 · FACT-CHECK — the district maps

**VISUAL:** The cascade max-depth map (`outputs/cascade_142/max_depth.png`)
with district outlines sketched over it. Highlight the thin flooded
corridor inside the huge Ernakulam outline. Then pan far south to an
empty Pathanamthitta outline; a ridge symbol between it and the Periyar.
[worst case, both dams gone]

**VO:**
Which brings me to the maps from the panic years — the ones with entire
districts erased. All of Ernakulam. Even Pathanamthitta. My worst case —
both dams gone — puts one to three metres in Ernakulam's river corridor
and backwaters. Genuinely dangerous. Not a district erased. And
Pathanamthitta? It's in a different river basin. There is a mountain
range in the way. The water would have to flow uphill. I can't simulate
it flooding, because physics declined.

---

## 3:13 – 3:53 · "MY MODEL IS TOO POLITE" — matching the famous study

**VISUAL:** Split screen: my arrival times vs published. Numbers in mono
font. [them: Idukki in ~2 h · me: ~7 h]

**VO:**
Now the self-audit — because in the other direction, my model looks too
slow. The famous IIT Roorkee study gets the flood to Idukki in about two
hours. Mine takes seven. Their headline number rests on one input: the
dam is assumed to disintegrate, completely, in twelve minutes. The
historical regression says two and a half hours.

**VISUAL:** The blog's animated hydrograph (screen-record it drawing):
teal Froehlich curve, then the red 12-minute spike. [same water, angrier
assumption]

**VO:**
So I gave my model their assumption. Peak outflow jumped fifty percent —
to within thirteen percent of their number. Arrival times fell by more
than half. Most of the famous urgency isn't a finding. It's an input.

---

## 3:53 – 4:23 · DIGGING A RIVER — the experiment that backfired

**VISUAL:** Terminal montage of the carve runs; the arrival-profile table;
a diagram of the trench undercutting the lake (simple animated sketch).
[experiment: carve a 100 m channel]

**VO:**
The rest of the gap is my terrain — remember the trees? I tried carving a
synthetic hundred-metre channel along the river. My first attempt built a
moat. My second made the flood slower. The best version helped by four
percent. You cannot origami your way out of bad terrain data — someone
should lidar that gorge. Until then, the truth sits between my seven
hours and their two.

---

## 4:23 – 4:48 · THE INDEPENDENT JURY

**VISUAL:** `outputs/comparison_maxdepth.png` — the three-panel
mine-vs-LISFLOOD figure. Then a quick stat card:
[two engines · same river · arrival within 5%]

**VO:**
Am I marking my own homework? I re-ran everything through LISFLOOD-FP —
a completely independent academic flood model. It agreed with me within
five percent. The published studies of this dam disagree with each other
by a factor of six. That spread is the dispute, in one number.

---

## 4:48 – 5:33 · CONCLUSION — the honest ledger

**VISUAL:** Three title cards over slow pans of the max-depth maps.

**VO:**
So, after flooding Kerala a few dozen times, here's my ledger.

**VISUAL:** [REAL] over the gorge map.

**VO:**
Real: the gorge. Every version of this failure devastates the valley below
the dam, fast. Their safety should be the fixed point of this argument —
not a bargaining chip inside it.

**VISUAL:** [CONDITIONAL] over the cascade map.

**VO:**
Conditional: drowned Kochi. It requires the cascade, and the cascade
requires a modern arch dam to fail because its reservoir rose five metres.
Possible? Worth planning for. But it's a compound hypothetical — and the
most safety-critical object in this story might actually be Idukki's
operating margin, not the old dam everyone photographs.

**VISUAL:** [INFLATED] over the two hydrograph curves.

**VO:**
Inflated: the certainty — on both sides. One side's two-hour wall of water
is a twelve-minute assumption. The other side's serenity leans on the same
missing knowledge, pointed the other way. Nobody knows how a
hundred-thirty-year-old masonry dam actually fails. That's the whole
dispute.

---

## 5:33 – 5:43 · OUTRO

**VISUAL:** The blog page scrolling slowly; the gauge fills. End card:
[blog.quantumautomata.in · github.com/shyam163/mullaperiyar-sim ·
not a hazard map. not engineering advice.]

**VO:**
Everything's reproducible — solver, data, every number — links below. And
if you live in that valley: your district's disaster plan outranks every
pixel of mine. The water, at least, has no politics. It just goes
downhill.

---

### Production notes
- Total VO ≈ 900 words → ~6:00 at 150 wpm; trim the build section first
  if running long.
- Register: dry, first-person, one sincere drop at Vandiperiyar and again
  at [REAL] — no music under those two beats.
- All sim footage exists: `outputs/{cascade_142,baseline_142}/animation.mp4`
  and `animation_3d.mp4`; maps in the same folders; LISFLOOD figure at
  `outputs/comparison_maxdepth.png`; the animated hydrograph and gauge can
  be screen-recorded from https://blog.quantumautomata.in/the-dam/.
- Fact-check beat data (verified from sudden_152 max_depth.tif, worst
  case): Ernakulam lowland corridor median 1.8 m, p95 3.1 m, max 4.6 m;
  zero water south of ~9.95 N in the lowlands; Pathanamthitta is in the
  Pamba basin, ~70 km south of the southernmost flooding, across the
  watershed divide.
- Caption every simulation shot with "simulation — order of magnitude
  only" in small text; keep the disclaimer card ≥3 s.
