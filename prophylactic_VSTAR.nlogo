;;==========================================================================
;; FMD SEIRC-V AGENT-BASED MODEL  --  PROPHYLACTIC VACCINATION AT REQUIRED COVERAGE (Cv = V* = HIT / VP)
;; File: prophylactic_VSTAR.nlogo   |   engine v3 (review revision)
;; The engine below is byte-for-byte identical in all strategy files; only the
;; strategy-specific block (vaccination rules, zone geometry) differs.
;;==========================================================================

;; ===== BEGIN COMMON ENGINE v3 =====

globals [
  ;; scenario (set by set-scenario / set-disease-parameters)
  current-serotype Ro VP
  exposure-period subclinical-period clinical-period carrier-duration
  carrier-probability n-inf-days

  ;; transmission
  beta-daily           ;; per-contact daily infection probability (calibrated, see calibrate-beta)
  calibrated-R0        ;; expected secondary infections per case at the calibration target
  lifetime-R0          ;; expected secondary infections per case, acute + carrier phase
  acute-R0             ;; expected secondary infections per case, acute phase only
  sample-klists        ;; kernel values of a random sample of animals (calibration)
  mean-contacts        ;; mean number of potential contacts per animal

  ;; domain
  r3 area-km2 density population simulation-days
  used-seed            ;; random seed actually used in this run (recorded for reproducibility)
  sim-day              ;; days completed (equals ticks after setup)
  focus-xs focus-ys    ;; coordinates (km) of the index animal(s)

  ;; vaccination bookkeeping
  HIT                  ;; 1 - 1/Ro                                           Eq. (1)
  vaccination-coverage ;; proportion of eligible animals vaccinated (prophylactic: V* = HIT/VP, Eq. 2)
  protected-count      ;; animals made protected by the primary dose
  booster-protected    ;; animals made (re-)protected by a booster dose after a lapse
  doses-given          ;; doses administered (primary + booster, responders and non-responders)
  first-vacc-day       ;; day of the primary vaccination (-1 = no vaccination)

  ;; compartment counts
  S-count E-count I-sub-count I-clin-count C-count R-count V-count

  ;; outputs
  cumulative-infections    ;; S/V -> E events after the index case(s)
  infections-from-carriers ;; infections caused by animals in the carrier state
  infections-in-vaccinated ;; infections of animals in state V (partial protection only)
  index-offspring          ;; infections caused directly by the index animal(s)
  cum-inf-d182 cum-inf-d364 ;; cumulative infections after 182 and 364 days (-1 if not reached)
  peak-clinical peak-day epidemic-duration
  attack-rate              ;; (cumulative-infections + n-foci) / population
  major-outbreak           ;; 1 if attack-rate > 0.05
  balance-error            ;; population - sum of compartments (must be 0)
]

turtles-own [
  state                ;; S, E, I_sub, I_clin, C, R, V
  days-in-state
  initial-source?
  vaccinated?          ;; has received a primary dose (responder or not)
  scheduled?           ;; selected for the reactive campaign (ring/barrier)
  protect-at           ;; day on which a dose becomes protective (-1 = none pending)
  dose-round           ;; 0 = primary dose, k = k-th booster
  nbr-list             ;; potential contacts
  nbr-k                ;; kernel K(d) = exp(-alpha d) for each contact
]

;;--------------------------------------------------------------------------
;; SETUP
;;--------------------------------------------------------------------------
to setup
  clear-all
  ifelse run-seed >= 0 [ set used-seed run-seed ] [ set used-seed new-seed ]
  random-seed used-seed
  set-scenario
  if global-sensitivity? [
    randomize-common-parameters
    strategy-randomize
  ]
  set-disease-parameters
  set simulation-days sim-days
  set sim-day 0
  set cumulative-infections 0
  set infections-from-carriers 0
  set infections-in-vaccinated 0
  set index-offspring 0
  set cum-inf-d182 -1
  set cum-inf-d364 -1
  set peak-clinical 0
  set peak-day 0
  set epidemic-duration 0
  set protected-count 0
  set booster-protected 0
  set doses-given 0
  set first-vacc-day -1
  set vaccination-coverage 0

  set-domain
  create-turtles population [ initialise-animal ]
  place-animals
  build-contacts
  calibrate-beta                               ;; herd is fully susceptible here
  strategy-setup                               ;; prophylactic vaccination happens here (t = 0)
  introduce-index
  paint-domain                                 ;; display only (needs the focus coordinates)

  update-compartment-counts
  reset-ticks
  update-outputs
  log-setup
end

to initialise-animal
  set size 0.6
  set shape "cow"
  set state "S"
  set color brown
  set days-in-state 0
  set initial-source? false
  set vaccinated? false
  set scheduled? false
  set protect-at -1
  set dose-round 0
  set nbr-list []
  set nbr-k []
end

to set-scenario
  ;; serotype-choice (chooser: "O" "A" "Asia1") and R0-choice (chooser: 2 4)
  set current-serotype serotype-choice
  set Ro R0-choice
end

;; durations (days) and vaccine protectiveness per serotype
to set-disease-parameters
  if current-serotype = "O" [
    set exposure-period 1.2  set subclinical-period 1.7  set clinical-period 6.6  set VP 0.823
  ]
  if current-serotype = "A" [
    set exposure-period 1.0  set subclinical-period 2.8  set clinical-period 7.5  set VP 0.767
  ]
  if current-serotype = "Asia1" [
    set exposure-period 1.5  set subclinical-period 2.0  set clinical-period 7.0  set VP 0.787
  ]
  set carrier-duration 320
  set carrier-probability 0.5
  ;; integer number of days an animal transmits in the acute phase: progression is tested as
  ;; days-in-state >= duration on an integer counter, so 1.7 d -> 2 d, 6.6 d -> 7 d, etc.
  set n-inf-days (ceiling subclinical-period) + (ceiling clinical-period)
end

;; global (multi-parameter) sensitivity analysis: every run draws its own parameter set from the
;; ranges below, using the run seed; the drawn values are exported as metrics
to randomize-common-parameters
  set alpha (0.5 + random-float 2.5)                          ;; 0.5 - 3.0 per km
  set n-contacts (10 + random 51)                             ;; 10 - 60
  set immunity-duration (90 + random 276)                     ;; 90 - 365 days
  set protection-delay (random 15)                            ;; 0 - 14 days
  set vaccine-leak (random-float 0.5)                         ;; 0 - 0.5
  set carrier-infectiousness (10 ^ (-3 + (random-float 3)))   ;; 0.001 - 1 (log-uniform)
end

;;--------------------------------------------------------------------------
;; DOMAIN: circular, r3 = 10 km (identical for all strategies, including prophylactic)
;;--------------------------------------------------------------------------
to set-domain
  set r3 10
  set area-km2 pi * (r3 ^ 2)
  if density-choice = "reference" [ set population 18536 ]
  if density-choice = "moderate"  [ set population 9268 ]
  if density-choice = "sparse"    [ set population 2001 ]
  set density population / area-km2
end

to place-animals
  ifelse herd-size > 1 [ place-clustered ] [ place-uniform ]
end

;; homogeneous (uniform random) placement in the disc
to place-uniform
  ask turtles [
    let rr r3 * (sqrt (random-float 1))
    let th random-float 360
    setxy (rr * (sin th)) (rr * (cos th))
  ]
end

;; herds: herd centres are uniform in the disc; each animal joins a random herd and is placed
;; around its centre with a normal scatter of herd-spread km (mean herd size = herd-size)
to place-clustered
  let n-herds max (list 1 (round (population / herd-size)))
  let centres []
  repeat n-herds [
    let rr r3 * (sqrt (random-float 1))
    let th random-float 360
    set centres lput (list (rr * (sin th)) (rr * (cos th))) centres
  ]
  ask turtles [
    let c one-of centres
    let px 0
    let py 0
    let ok false
    while [ not ok ] [
      set px (item 0 c) + (random-normal 0 herd-spread)
      set py (item 1 c) + (random-normal 0 herd-spread)
      set ok ((sqrt ((px * px) + (py * py))) <= r3)
    ]
    setxy px py
  ]
end

;;--------------------------------------------------------------------------
;; CONTACT STRUCTURE AND TRANSMISSION CALIBRATION
;;--------------------------------------------------------------------------
;; contact-mode "nearest-n": each animal has n-contacts potential contacts, its nearest animals
;; (any state).  contact-mode "radius": all animals within contact-radius km (the number of
;; contacts then varies with density).  Daily infection probability from infectious i to
;; susceptible j is beta-daily * K(d_ij), K(d) = exp(-alpha d)  (d in km; Eq. 4).
to build-contacts
  ifelse contact-mode = "radius" [ build-contacts-radius ] [ build-contacts-nearest ]
  set mean-contacts mean [ length nbr-list ] of turtles
end

to build-contacts-nearest
  if n-contacts >= population [
    user-message "n-contacts must be smaller than the population."
    stop
  ]
  let search-r 1.5 * (sqrt (n-contacts / (pi * density)))
  ask turtles [
    let r search-r
    let cand other turtles in-radius r
    while [ count cand < n-contacts ] [
      set r r * 1.5
      set cand other turtles in-radius r
    ]
    let nearest min-n-of n-contacts cand [ distance myself ]
    set nbr-list sort-on [ distance myself ] nearest
    set nbr-k map [ nb -> exp (0 - (alpha * (distance nb))) ] nbr-list
  ]
end

to build-contacts-radius
  ask turtles [
    let cand other turtles in-radius contact-radius
    set nbr-list sort-on [ distance myself ] cand
    set nbr-k map [ nb -> exp (0 - (alpha * (distance nb))) ] nbr-list
  ]
end

;; probability that a given contact j is infected at some time during the infectious life of i:
;; the acute phase, plus (if incl-carrier? is true) the carrier phase of animals that become carriers
to-report pair-infection-probability [ b k incl-carrier? ]
  let q min (list 1 (b * k))
  let p-acute 1 - ((1 - q) ^ n-inf-days)
  let p-carrier 0
  if incl-carrier? [
    set p-carrier 1 - ((1 - (min (list 1 (carrier-infectiousness * q)))) ^ carrier-duration)
  ]
  report 1 - ((1 - p-acute) * (1 - (carrier-probability * p-carrier)))
end

to-report expected-secondary [ b incl-carrier? ]
  report mean map [ ks -> sum map [ k -> pair-infection-probability b k incl-carrier? ] ks ] sample-klists
end

;; Solve for beta-daily (bisection) so that the expected number of secondary infections caused by
;; one infected animal in a fully susceptible herd equals Ro.
;;   calibrate-acute-only? = false : Ro is the LIFETIME value (acute + carrier phase)
;;   calibrate-acute-only? = true  : Ro is the ACUTE-phase value; carriers add transmission on top,
;;                                   so carrier-infectiousness can be varied without changing the
;;                                   acute-phase transmission
to calibrate-beta
  let sample-size min (list 300 (count turtles))
  set sample-klists [ nbr-k ] of (n-of sample-size turtles)
  let incl (not calibrate-acute-only?)
  ifelse (expected-secondary 1 incl) < Ro [
    set beta-daily 1
    show (word "WARNING: Ro = " Ro " cannot be reached with these contacts; beta-daily set to 1.")
  ] [
    let lo 0
    let hi 1
    repeat 50 [
      let mid (lo + hi) / 2
      ifelse (expected-secondary mid incl) < Ro [ set lo mid ] [ set hi mid ]
    ]
    set beta-daily (lo + hi) / 2
  ]
  set calibrated-R0 expected-secondary beta-daily incl
  set lifetime-R0 expected-secondary beta-daily true
  set acute-R0 expected-secondary beta-daily false
end

;;--------------------------------------------------------------------------
;; INDEX CASE(S), VACCINATION, STATE CHANGE
;;--------------------------------------------------------------------------
;; the index animal is the susceptible animal nearest the domain centre (index-location "centre")
;; or a random susceptible animal ("random"); n-foci > 1 adds further random index animals.
;; index-stage is the state of the index animal on day 0 (I_clin, or E for a full infectious life).
to introduce-index
  set focus-xs []
  set focus-ys []
  let cand turtles with [ state = "S" ]
  let src nobody
  ifelse index-location = "random"
    [ set src one-of cand ]
    [ set src min-one-of cand [ distancexy 0 0 ] ]
  if src != nobody [ make-index src ]
  repeat (n-foci - 1) [
    let more turtles with [ state = "S" ]
    if any? more [ make-index (one-of more) ]
  ]
end

to make-index [ a ]
  ask a [
    set-state index-stage
    set initial-source? true
    set size 1.2
    set focus-xs lput xcor focus-xs
    set focus-ys lput ycor focus-ys
  ]
end

to set-state [ new-state ]
  set state new-state
  set days-in-state 0
  if new-state = "S"      [ set color brown ]
  if new-state = "E"      [ set color yellow ]
  if new-state = "I_sub"  [ set color orange ]
  if new-state = "I_clin" [ set color red ]
  if new-state = "C"      [ set color gray ]
  if new-state = "R"      [ set color cyan ]
  if new-state = "V"      [ set color green ]
end

;; primary dose (reactive campaigns): the animal becomes protected (state V) with probability VP,
;; after protection-delay days; otherwise it stays susceptible (vaccine failure).
;; Only susceptible animals are vaccinated.
to vaccinate
  set vaccinated? true
  set doses-given doses-given + 1
  if random-float 1 < VP [ schedule-protection 0 ]
end

to schedule-protection [ k ]
  set dose-round k
  ifelse protection-delay <= 0
    [ become-protected ]
    [ set protect-at sim-day + protection-delay ]
end

to become-protected
  set-state "V"
  set protect-at -1
  ifelse dose-round = 0
    [ set protected-count protected-count + 1 ]
    [ set booster-protected booster-protected + 1 ]
end

;; booster rounds: every booster-interval days after the primary vaccination (n-boosters rounds),
;; a fraction booster-coverage of the animals that received the primary dose and are still
;; susceptible or protected is re-vaccinated.  A protected animal keeps protection and its
;; protection clock is restarted; a susceptible animal is protected with probability VP
;; (after protection-delay days).
to booster-daily
  if (booster-interval > 0) and (first-vacc-day >= 0) and (sim-day > first-vacc-day) [
    let since sim-day - first-vacc-day
    if (since mod booster-interval = 0) and ((since / booster-interval) <= n-boosters) [
      let k round (since / booster-interval)
      ask turtles with [ vaccinated? and (state = "S" or state = "V") ] [
        if random-float 1 < booster-coverage [ boost k ]
      ]
    ]
  ]
end

to boost [ k ]
  set doses-given doses-given + 1
  ifelse state = "V"
    [ set days-in-state 0 ]
    [ if protect-at < 0 [ if random-float 1 < VP [ schedule-protection k ] ] ]
end

;; relative susceptibility of a protected (V) animal d days after protection started:
;; vaccine-leak (0 = fully protected) until waning-window days before the end of protection, then
;; rising linearly to 1 (gradual waning); waning-window = 0 gives abrupt loss of protection
to-report susceptibility-of-protected [ d ]
  ifelse waning-window <= 0
    [ report vaccine-leak ]
    [
      let f (d - (immunity-duration - waning-window)) / waning-window
      report vaccine-leak + ((1 - vaccine-leak) * (max (list 0 (min (list 1 f)))))
    ]
end

to become-exposed
  set-state "E"
  set cumulative-infections cumulative-infections + 1
end

;;--------------------------------------------------------------------------
;; GO  (one tick = one day)
;;   1 strategy-daily (reactive primary vaccination)  2 booster-daily  3 progression
;;   4 transmission  5 outputs
;;--------------------------------------------------------------------------
to go
  if ticks >= simulation-days [ stop ]
  strategy-daily
  booster-daily
  update-states
  transmit
  update-compartment-counts
  update-outputs
  tick
  set sim-day sim-day + 1
  if sim-day = 182 [ set cum-inf-d182 cumulative-infections ]
  if sim-day = 364 [ set cum-inf-d364 cumulative-infections ]
end

;; time-based progression on an integer day counter (fractional durations act as the next
;; whole day); an animal that changes state on a tick does not change again on that tick.
;; Protection lasts exactly immunity-duration days of transmission: the counter of a protected
;; animal is 1 on its first protected day, and protection ends when it exceeds immunity-duration.
to update-states
  ask turtles [
    if protect-at >= 0 and sim-day >= protect-at [
      set protect-at -1
      if state = "S" [ become-protected ]
    ]
    set days-in-state days-in-state + 1
    if state = "V" and days-in-state > immunity-duration [
      set-state "S"
    ]
    if state = "E" and days-in-state >= exposure-period [
      set-state "I_sub"
    ]
    if state = "I_sub" and days-in-state >= subclinical-period [
      set-state "I_clin"
    ]
    if state = "I_clin" and days-in-state >= clinical-period [
      ifelse random-float 1 < carrier-probability [ set-state "C" ] [ set-state "R" ]
    ]
    if state = "C" and days-in-state >= carrier-duration [
      set-state "R"
    ]
  ]
end

;; subclinical and clinical animals infect each susceptible contact with probability
;; beta-daily * K(d); carriers do so with carrier-infectiousness * beta-daily * K(d).
;; Protected (V) animals are infected with probability multiplied by susceptibility-of-protected
;; (0 in the baseline = all-or-nothing protection).  Recovered and already-infected animals
;; cannot be infected.
to transmit
  ask turtles with [ state = "I_sub" or state = "I_clin" or state = "C" ] [
    let from-carrier? (state = "C")
    let rel ifelse-value from-carrier? [ carrier-infectiousness ] [ 1 ]
    if rel > 0 [
      (foreach nbr-list nbr-k [ [ nb k ] ->
        let nb-state [ state ] of nb
        let susc 0
        if nb-state = "S" [ set susc 1 ]
        if nb-state = "V" [ set susc susceptibility-of-protected ([ days-in-state ] of nb) ]
        if susc > 0 [
          if random-float 1 < (rel * beta-daily * k * susc) [
            if nb-state = "V" [ set infections-in-vaccinated infections-in-vaccinated + 1 ]
            if from-carrier? [ set infections-from-carriers infections-from-carriers + 1 ]
            if initial-source? [ set index-offspring index-offspring + 1 ]
            ask nb [ become-exposed ]
          ]
        ]
      ])
    ]
  ]
end

;;--------------------------------------------------------------------------
;; GEOMETRY HELPERS (distance to the index animal(s))
;;--------------------------------------------------------------------------
to-report within-radius-of-focus? [ rad ]
  report member? true (map [ [ fx fy ] -> (distancexy fx fy) <= rad ] focus-xs focus-ys)
end

to-report within-annulus-of-focus? [ rin rout ]
  report member? true (map [ [ fx fy ] -> ((distancexy fx fy) >= rin) and ((distancexy fx fy) <= rout) ] focus-xs focus-ys)
end

;;--------------------------------------------------------------------------
;; OUTPUTS AND CHECKS
;;--------------------------------------------------------------------------
to update-compartment-counts
  set S-count      count turtles with [ state = "S" ]
  set E-count      count turtles with [ state = "E" ]
  set I-sub-count  count turtles with [ state = "I_sub" ]
  set I-clin-count count turtles with [ state = "I_clin" ]
  set C-count      count turtles with [ state = "C" ]
  set R-count      count turtles with [ state = "R" ]
  set V-count      count turtles with [ state = "V" ]
end

to update-outputs
  if I-clin-count > peak-clinical [
    set peak-clinical I-clin-count
    set peak-day ticks
  ]
  if (E-count + I-sub-count + I-clin-count) > 0 [ set epidemic-duration ticks ]
  set attack-rate (cumulative-infections + n-foci) / population
  set major-outbreak ifelse-value (attack-rate > 0.05) [ 1 ] [ 0 ]
  set balance-error population - (S-count + E-count + I-sub-count + I-clin-count + C-count + R-count + V-count)
end

to-report protected-fraction
  report protected-count / population
end

to log-setup
  show (word "Scenario: serotype " current-serotype " | Ro = " Ro " | VP = " VP " | seed = " used-seed)
  show (word "Population = " population " | density = " (precision density 2) " per km2 | mean contacts = " (precision mean-contacts 1) " | alpha = " alpha)
  show (word "beta-daily = " (precision beta-daily 5) " | acute R0 = " (precision acute-R0 3) " | lifetime R0 = " (precision lifetime-R0 3) " (target " Ro ", acute-only calibration = " calibrate-acute-only? ")")
  show (word "Carrier infectiousness = " carrier-infectiousness " | protection = " immunity-duration " d, delay " protection-delay " d, leak " vaccine-leak ", waning window " waning-window " d")
  show (word "Booster interval = " booster-interval " d, rounds = " n-boosters ", coverage = " booster-coverage " | horizon = " sim-days " d")
end

;; press CHECK after SETUP (or at any time) to run the built-in verification tests
to check-model
  update-compartment-counts
  update-outputs
  show (word "Serotype " current-serotype " (chooser: " serotype-choice "): latent " exposure-period " d, subclinical " subclinical-period " d, clinical " clinical-period " d, VP = " VP ", Ro = " Ro " (chooser: " R0-choice ")")
  if current-serotype != serotype-choice or Ro != R0-choice [ show "ERROR: serotype or Ro does not match the choosers" ]
  show (word "Population balance error = " balance-error "  (must be 0)")
  show (word "Calibrated expected secondary infections = " (precision calibrated-R0 3) "  (target Ro = " Ro ")")
  show (word "Doses = " doses-given " | protected by primary dose = " protected-count " | by booster = " booster-protected " | Peff(t=0) = " (precision protected-fraction 4))
  show (word "Mean number of contacts = " (precision mean-contacts 2))
  show (word "Index animals = " (count turtles with [ initial-source? ]) "  (must equal n-foci = " n-foci ")")
end

;; ===== END COMMON ENGINE v3 =====

;;##########################################################################
;; STRATEGY-SPECIFIC BLOCK: prophylactic_VSTAR
;; (procedures required by the engine: strategy-setup, strategy-daily,
;;  strategy-randomize, paint-domain; ring/barrier also in-vaccination-zone?)
;;##########################################################################

to paint-domain
  ask patches [ ifelse (distancexy 0 0) > r3 [ set pcolor gray + 3 ] [ set pcolor white ] ]
end

;; prophylactic vaccination at t = 0 (before the index animal is introduced); protection is
;; assumed to be established by the time the virus arrives, so protection-delay does not apply to
;; the primary dose (it does apply to booster doses given to animals that lost protection)
to strategy-daily
  ;; no reactive vaccination in this strategy
end

to strategy-randomize
  ;; no strategy-specific parameters
end

to prophylactic-vaccination
  set first-vacc-day 0
  ask n-of (round (vaccination-coverage * population)) turtles [
    set vaccinated? true
    set doses-given doses-given + 1
    if random-float 1 < VP [
      set dose-round 0
      become-protected
    ]
  ]
  show (word "Prophylactic vaccination at t = 0: coverage Cv = " (precision vaccination-coverage 4) " | doses = " doses-given " | protected = " protected-count " | Peff = " (precision protected-fraction 4) " (HIT = " (precision HIT 3) ")")
end

;; Cv = V* = HIT / VP (Eq. 2): expected protected fraction Peff = V* x VP = HIT (Eq. 3)
to strategy-setup
  set HIT 1 - (1 / Ro)
  set vaccination-coverage min (list 1 (HIT / VP))
  if (HIT / VP) > 1 [ show "WARNING: HIT / VP > 1; coverage capped at 1, HIT is not attainable." ]
  prophylactic-vaccination
end

@#$#@#$#@
GRAPHICS-WINDOW
215
10
635
430
-1
-1
20.0
1
10
1
1
1
0
0
0
1
-10
10
-10
10
0
0
1
ticks
30.0

BUTTON
5
10
70
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
75
10
140
43
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
145
10
205
43
NIL
check-model
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
645
75
1330
430
FMD dynamics (animals per state)
day
animals
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Susceptible" 1.0 0 -4757638 true "" "plot S-count"
"Exposed" 1.0 0 -4079321 true "" "plot E-count"
"Subclinical" 1.0 0 -955883 true "" "plot I-sub-count"
"Clinical" 1.0 0 -2674135 true "" "plot I-clin-count"
"Carrier" 1.0 0 -7500403 true "" "plot C-count"
"Recovered" 1.0 0 -11221820 true "" "plot R-count"
"Vaccinated" 1.0 0 -10899396 true "" "plot V-count"

MONITOR
645
10
723
55
Day
ticks
0
1
11

MONITOR
728
10
806
55
S
S-count
0
1
11

MONITOR
811
10
889
55
E
E-count
0
1
11

MONITOR
894
10
972
55
I sub
I-sub-count
0
1
11

MONITOR
977
10
1055
55
I clin
I-clin-count
0
1
11

MONITOR
1060
10
1138
55
Carriers
C-count
0
1
11

MONITOR
1143
10
1221
55
R
R-count
0
1
11

MONITOR
1226
10
1304
55
V
V-count
0
1
11

MONITOR
215
440
311
485
Attack rate
attack-rate
3
1
11

MONITOR
315
440
411
485
Cumulative infections
cumulative-infections
0
1
11

MONITOR
415
440
511
485
Peak clinical
peak-clinical
0
1
11

MONITOR
515
440
611
485
Peak day
peak-day
0
1
11

MONITOR
615
440
711
485
Peff (V / N)
protected-fraction
3
1
11

MONITOR
715
440
811
485
Doses given
doses-given
0
1
11

MONITOR
815
440
911
485
Calibrated R0
calibrated-R0
2
1
11

MONITOR
915
440
1011
485
beta (per day)
beta-daily
4
1
11

MONITOR
1015
440
1130
485
Infections from carriers
infections-from-carriers
0
1
11

MONITOR
1135
440
1250
485
Random seed
used-seed
0
1
11

CHOOSER
5
500
205
545
serotype-choice
serotype-choice
"O" "A" "Asia1"
0

CHOOSER
215
500
415
545
R0-choice
R0-choice
2 4
0

CHOOSER
425
500
625
545
density-choice
density-choice
"reference" "moderate" "sparse"
1

SLIDER
635
500
835
533
carrier-infectiousness
carrier-infectiousness
0
1
0.01
0.01
1
NIL
HORIZONTAL

SWITCH
845
500
1045
533
calibrate-acute-only?
calibrate-acute-only?
1
1
-1000

SLIDER
1055
500
1255
533
alpha
alpha
0
5
1.5
0.1
1
NIL
HORIZONTAL

SLIDER
5
550
205
583
n-contacts
n-contacts
5
100
30
1
1
NIL
HORIZONTAL

CHOOSER
215
550
415
595
contact-mode
contact-mode
"nearest-n" "radius"
0

SLIDER
425
550
625
583
contact-radius
contact-radius
0.1
2
0.4
0.05
1
NIL
HORIZONTAL

SLIDER
635
550
835
583
herd-size
herd-size
1
200
1
1
1
NIL
HORIZONTAL

SLIDER
845
550
1045
583
herd-spread
herd-spread
0.02
0.5
0.1
0.01
1
NIL
HORIZONTAL

SLIDER
1055
550
1255
583
immunity-duration
immunity-duration
30
730
182
1
1
NIL
HORIZONTAL

SLIDER
5
600
205
633
protection-delay
protection-delay
0
28
0
1
1
NIL
HORIZONTAL

SLIDER
215
600
415
633
vaccine-leak
vaccine-leak
0
1
0
0.05
1
NIL
HORIZONTAL

SLIDER
425
600
625
633
waning-window
waning-window
0
365
0
5
1
NIL
HORIZONTAL

SLIDER
635
600
835
633
booster-interval
booster-interval
0
730
0
1
1
NIL
HORIZONTAL

SLIDER
845
600
1045
633
n-boosters
n-boosters
0
200
1
1
1
NIL
HORIZONTAL

SLIDER
1055
600
1255
633
booster-coverage
booster-coverage
0
1
1
0.05
1
NIL
HORIZONTAL

SLIDER
5
650
205
683
sim-days
sim-days
100
36500
730
10
1
NIL
HORIZONTAL

CHOOSER
215
650
415
695
index-location
index-location
"centre" "random"
0

CHOOSER
425
650
625
695
index-stage
index-stage
"I_clin" "E"
0

SLIDER
635
650
835
683
n-foci
n-foci
1
3
1
1
1
NIL
HORIZONTAL

SLIDER
845
650
1045
683
run-seed
run-seed
-1
1000000
-1
1
1
NIL
HORIZONTAL

SWITCH
1055
650
1255
683
global-sensitivity?
global-sensitivity?
1
1
-1000

@#$#@#$#@
## WHAT IS IT?

PROPHYLACTIC VACCINATION AT REQUIRED COVERAGE (Cv = V* = HIT / VP). FMD SEIRC-V agent-based model, engine v3 (identical code in all strategy files; only the strategy block differs). Closed cattle population in a circular domain of radius 10 km (reference 18,536, moderate 9,268, sparse 2,001 animals). The same domain, contact network, calibration and disease parameters are used for the prophylactic, ring, barrier and control files, so strategies can be compared directly.

## HOW IT WORKS

Contacts: nearest n-contacts animals (contact-mode nearest-n) or all animals within contact-radius km (radius). Transmission probability beta-daily x exp(-alpha d). beta-daily is calibrated so the expected number of secondary infections per case equals Ro (lifetime = acute + carrier phase; or acute phase only if calibrate-acute-only? is on). Vaccination is all-or-nothing: a dose protects with probability VP (serotype specific) after protection-delay days, for immunity-duration days (waning-window > 0 gives gradual loss of protection; vaccine-leak > 0 gives partial protection). Booster doses are given every booster-interval days (n-boosters rounds) to the animals that received the primary dose. Random seeds are recorded (used-seed); set run-seed >= 0 to reproduce a run.

## HOW TO USE IT

Choose the inputs, press SETUP, then GO; CHECK runs the built-in tests. BehaviorSpace experiments are stored in the file; the run length of each experiment (730 days) must match sim-days.

@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="prophylactic_VSTAR_main" repetitions="100" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
      <value value="&quot;A&quot;"/>
      <value value="&quot;Asia1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;reference&quot;"/>
      <value value="&quot;moderate&quot;"/>
      <value value="&quot;sparse&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="prophylactic_VSTAR_timeseries" repetitions="100" sequentialRunOrder="true" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>doses-given</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
      <value value="&quot;A&quot;"/>
      <value value="&quot;Asia1&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;reference&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_carrier-infectiousness" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0"/>
      <value value="0.01"/>
      <value value="0.1"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_carrier-independent" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0"/>
      <value value="0.01"/>
      <value value="0.1"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_alpha" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="0.5"/>
      <value value="1.5"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_n-contacts" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="10"/>
      <value value="30"/>
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_contact-radius" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;reference&quot;"/>
      <value value="&quot;moderate&quot;"/>
      <value value="&quot;sparse&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
      <value value="&quot;radius&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_herd-clustering" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
      <value value="20"/>
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_index-location" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
      <value value="&quot;random&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_immunity-duration" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="120"/>
      <value value="182"/>
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_vaccine-leak" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
      <value value="0.1"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_booster" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
      <value value="365"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="sens_booster-delay" repetitions="50" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
      <value value="7"/>
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="global_sensitivity" repetitions="500" sequentialRunOrder="true" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="730"/>
    <metric>cumulative-infections</metric>
    <metric>attack-rate</metric>
    <metric>major-outbreak</metric>
    <metric>peak-clinical</metric>
    <metric>peak-day</metric>
    <metric>epidemic-duration</metric>
    <metric>protected-count</metric>
    <metric>booster-protected</metric>
    <metric>doses-given</metric>
    <metric>infections-from-carriers</metric>
    <metric>infections-in-vaccinated</metric>
    <metric>index-offspring</metric>
    <metric>cum-inf-d182</metric>
    <metric>cum-inf-d364</metric>
    <metric>S-count</metric>
    <metric>E-count</metric>
    <metric>I-sub-count</metric>
    <metric>I-clin-count</metric>
    <metric>C-count</metric>
    <metric>R-count</metric>
    <metric>V-count</metric>
    <metric>balance-error</metric>
    <metric>beta-daily</metric>
    <metric>calibrated-R0</metric>
    <metric>acute-R0</metric>
    <metric>lifetime-R0</metric>
    <metric>mean-contacts</metric>
    <metric>used-seed</metric>
    <metric>alpha</metric>
    <metric>n-contacts</metric>
    <metric>immunity-duration</metric>
    <metric>protection-delay</metric>
    <metric>vaccine-leak</metric>
    <metric>carrier-infectiousness</metric>
    <enumeratedValueSet variable="serotype-choice">
      <value value="&quot;O&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="R0-choice">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="density-choice">
      <value value="&quot;moderate&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="carrier-infectiousness">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="calibrate-acute-only?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="alpha">
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-contacts">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-mode">
      <value value="&quot;nearest-n&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contact-radius">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-size">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="herd-spread">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunity-duration">
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="protection-delay">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vaccine-leak">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="waning-window">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-interval">
      <value value="0"/>
      <value value="182"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-boosters">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="booster-coverage">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sim-days">
      <value value="730"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-location">
      <value value="&quot;centre&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="index-stage">
      <value value="&quot;I_clin&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="n-foci">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="run-seed">
      <value value="-1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="global-sensitivity?">
      <value value="true"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
