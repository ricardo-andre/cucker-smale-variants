globals [
  ; discretization time interval
  time-scale
  ; near zero speed
  speed0

  ; runs
  ;; random seed needed for reprudicing known results
  initial-random-seed
  ;; stats
  ; bird furthest from the centre, the actual diameter is at most twice of this
  stat-max-radius
  ; fastest bird from the centre, the actual expansion is at most twice of this
  stat-max-speed

  ; errors
  ; by drifting in the EDO solution
  error-centre
  error-velocity ;in the default case will be zero, because the variation is simmetric
  error-turtles-at-the-wall
]

turtles-own [
  ; NetLogo already keeps current coordinates and heading
  ; we could just keep the speed
  ; but to avoid repeating on and on calculations like
  ;  `dx * speed` or `sin heading * speed`
  ;  we will keep track of velocity and then update the heading
  ; yes speed^2 = vx^2 + vy^2 but we'll waste memory in exchange of speed
  speed ; `speed` seem to cause some bug?
  vx
  vy

  ; the velocity change, deltav
  deltavx
  deltavy

  ;TODO:20230724:ricardo-andre:extend to Störmer-Verlet
  ; to improve the calculation, using something better
  ; than a Euler method for the differential equation
  ; e.g. a Störmer-Verlet method
  ; we will need to keep track of past coordinates and velocity
  ; maybe as lists (or arrays?)
  ;speeds
  ;headings
  ;coordinates
  done ; just a marker
]

links-own [
  weight0 ;the base weigth to be normalized
  weight1 ;psi_i (x_i, x_j)
  weight2 ;psi_j (x_i, x_j), for the non-symmetric case
]


to startup
  ; on first interactive run
  ; loads defaults for...
  ; ... user parameters
  set population 10
  set radius 10
  set topology-type "random"
  set vel-mean 1
  set vel-stdev 0.2
  set dimension "2D"
  set psi-type "standard"
  set beta 1
  set K 1
  ; ... world parameters
  setup
end

to setup
  clear-all
  resize-world -35 35 -35 35
  set-patch-size 7
  ; white backgroud
  ask patches [ set pcolor white ]

  ; globals
  set speed0 10 ^ -3
  set time-scale 0.1

  set initial-random-seed new-seed
  random-seed initial-random-seed

  if (population mod 2 = 1) [set population population + 1]
  create-turtles population [
    ; use all dark-ish colors
    set color 10 * (random 14) + (1 + random-float 4.5)
    set size world-width / 50
  ]

  ;TODO:20230730:ricardo-andre:improve mean & stdev of velocity - currently done in pairs
  ;no need to zero velocities because they are created in pairs
  (ifelse
    topology-type = "random" [
      turtles-random-position
      zero-sum-of-coordinates
      set-turtles-initial-velocity-by-pairs
    ]
    topology-type = "collision" [
      turtles-collision-startup
    ]
    topology-type = "bi-cluster" [
      turtles-bi-startup
  ])


  ; link all turtles, allows to use link-length, and to compute over links
  ;  and possibly link-heading and more later on
  ask turtles [
    create-links-with other turtles [
      set hidden? true
    ]
  ]

  update-stats

  reset-ticks
end

to turtles-random-position
  ask turtles [
    (ifelse
      dimension = "2D" [
        let r radius * sqrt random-float 1
        let theta random-float 360
        setxy r * sin theta r * cos theta
      ]
      dimension = "1D" [
        let x random-float (2 * radius) - radius
        setxy x 0
    ])
  ]
end

;to set-turtles-initial-velocity
  ;TODO:20230724:ricardo-andre:better and chooseable speed distribution
  ;TODO:20231106:ricardo-andre:currentely unused, done in pairs
;  set speed random-gamma ((vel-mean / vel-stdev) ^ 2) (vel-mean / vel-stdev ^ 2)
;  (ifelse
;    dimension = "2D" [
;      set vx speed * sin heading
;      set vy speed * cos heading
;    ]
;    dimension = "1D" [
;      ;TODO:20231106:ricardo-andre:BUG needs to set sign and heading
;      set vx speed
;      set vy 0
;  ])
;end

to set-turtles-initial-velocity-by-pairs
  ask turtles [set done 0]
  ask turtles [
    if (done = 0) [
      set done 1
      set-speed-random-gamma
      let mypair one-of other turtles with [done = 0]
      ask mypair [
        set done 1
        set speed [speed] of myself
        set heading (180 + [heading] of myself)
        set vx (- [vx] of myself)
        set vy (- [vy] of myself)
      ]
    ]
  ]
end

; all turtles face the center!!
to turtles-collision-startup
  ask turtles [set done 0]
  ask turtles [
    if (done = 0) [
      set done 1
      (ifelse
        dimension = "2D" [
          let r radius * sqrt random-float 1
          let theta random-float 180
          setxy r * sin theta r * cos theta
        ]
        dimension = "1D" [
          let x random-float radius
          setxy x 0
      ])
      facexy 0 0
      set-speed-random-gamma
      let mypair one-of other turtles with [done = 0]
      ask mypair [
        set done 1
        set speed [speed] of myself
        set heading (180 + [heading] of myself)
        set xcor (- [xcor] of myself)
        set ycor (- [ycor] of myself)
        set vx (- [vx] of myself)
        set vy (- [vy] of myself)
      ]
    ]
  ]
end

to turtles-bi-startup
  ask turtles [set done 0]
  ask turtles [
    if (done = 0) [
      set done 1
      (ifelse
        dimension = "2D" [
          set xcor (radius / 2) +  random-float (radius / 10)
          set ycor (- radius / 20) + random-float (radius / 10)
        ]
        dimension = "1D" [
          set xcor (radius / 2) +  random-float (radius / 10)
          set ycor 0
      ])
      set heading -45 + random-float 90
      set-speed-random-gamma
      let mypair one-of other turtles with [done = 0]
      ask mypair [
        set done 1
        set speed [speed] of myself
        set heading (180 + [heading] of myself)
        set xcor (- [xcor] of myself)
        set ycor (- [ycor] of myself)
        set vx (- [vx] of myself)
        set vy (- [vy] of myself)
      ]
    ]
  ]
end


to set-speed-random-gamma
  set speed random-gamma ((vel-mean / vel-stdev) ^ 2) (vel-mean / vel-stdev ^ 2)
  (ifelse
    dimension = "2D" [
      set vx speed * sin heading
      set vy speed * cos heading
    ]
    dimension = "1D" [
      set vx speed
      set vy 0
  ])
end

to go
  ;; compute the rate of chage of velocity, dotv
  ; communication weights
  ask links [ calc-link-weights-psi ]
  if (psi-type = "normalized") [
    calc-link-weights-normalize
  ]

  ; new deltav
  ask turtles [
    set deltavx 0
    set deltavy 0
  ]
  ; the links compute each term of deltav and ask turtles to store it
  ; so that we do not have to recompute the (same) deltavx & deltavy from each end
  ask links [ calc-deltav-of-turtles ]
  ; update vx, vy, speed & heading
  ask turtles [
    set vx vx + deltavx * time-scale
    set vy vy + deltavy * time-scale
    set speed sqrt (vx ^ 2 + vy ^ 2)
    set heading atan vx vy
    ; (almost) stopped turtles lose direction
    if (speed < speed0) and (shape = "default") [
      set size size / 2
      set shape "circle"
    ]
  ]

  ; and move! (that is the first EDO)
  ask turtles [ fd speed * time-scale ]

  update-stats
  tick

  ;  ask error-turtles-at-the-wall [
  ;    set vx 0
  ;    set vy 0
  ;  ]

  ; invalid model warnig
  if (any? error-turtles-at-the-wall) [
    ask patches [ set pcolor 139 ]
  ]
  ; stop execution...
  ; ...if all turtles are almost stopped
  ; or
  ; ...if too many turtles hit the walls
  if (
    (not any? turtles with [ speed >= speed0 ])
    or
    (count error-turtles-at-the-wall > population / 10)
  ) [stop]
end



;; mean-velocities
; reports a list with mean of x and mean of y corrdinates of all turtles
to-report mean-velocities
  let meanvx (sum [vx] of turtles) / population
  let meanvy (sum [vy] of turtles) / population
  report list meanvx meanvy
end

;; zero-sum-of-coordinates
; adds a constant to all coordinates so that the sum is zero
; because this expressed at the centre of mass
to zero-sum-of-coordinates
  let mcor mean-coordinates
  let meanx (item 0 mcor)
  let meany (item 1 mcor)
  ask turtles [
    set xcor xcor - meanx
    set ycor ycor - meany
  ]
end

;; mean-coordinates
; reports a list with mean of x and mean of y corrdinates of all turtles
to-report mean-coordinates
  let meanx (sum [xcor] of turtles) / population
  let meany (sum [ycor] of turtles) / population
  report list meanx meany
end

;; calc-deltav-of-turtles
; called by links
; computes deltav and asks each end to update it's own
to calc-deltav-of-turtles
  let dvx [vx] of end2 - [vx] of end1
  let dvy [vy] of end2 - [vy] of end1
  ask end1 [
    set deltavx deltavx + [weight1] of myself * dvx
    set deltavy deltavy + [weight1] of myself * dvy
  ]
  ask end2 [
    set deltavx deltavx - [weight2] of myself * dvx
    set deltavy deltavy - [weight2] of myself * dvy
  ]
end

;; update-stats
; updates the stats to display
to update-stats
  ; global flocking stats
  set stat-max-radius max [sqrt(xcor ^ 2 + ycor ^ 2)] of turtles
  set stat-max-speed max [sqrt(vx ^ 2 + vy ^ 2)] of turtles

  ; error checking stats
  ; if these get far from 1 the the Euler method is drifting away...
  let mcor mean-coordinates
  set error-centre sqrt ((item 0 mcor) ^ 2 + (item 1 mcor) ^ 2)
  let mvel mean-velocities
  set error-velocity sqrt ((item 0 mvel) ^ 2 + (item 1 mvel) ^ 2)

  set error-turtles-at-the-wall turtles with [
    xcor < min-pxcor or xcor > max-pxcor or ycor < min-pycor or ycor > max-pycor
  ]
end


;; communication weights

;; calc-link-weights-psi
; called by links
;  calls one of the specific psi function
;  for "normalized" weights this is only a first step
to calc-link-weights-psi
  (ifelse
    psi-type = "standard" [
      set weight1 k * (psi-standard link-length) / population
      set weight2 weight1
    ]
    psi-type = "normalized" [
      set weight0 psi-standard link-length
    ]
    psi-type = "singular" [
      set weight1 k * (psi-singular link-length) / population
      set weight2 weight1
  ])
end

;; calc-link-weights-normalize
; called by links
;  second step in weights calculation for "normalized" only
to calc-link-weights-normalize
  let S 0
  ask turtles [
    set S 0
    set S S + (sum [weight0] of links with [end1 = myself])
    set S S + (sum [weight0] of links with [end2 = myself])
    ask links with [end1 = myself] [
      set weight1 (k * weight0 / S)
    ]
    ask links with [end2 = myself] [
      set weight2 (k * weight0 / S)
    ]
  ]
end

;; psi-standard
; the standard communication function
; radial, positive, finite
to-report psi-standard [dist]
  report 1 / (1 + dist ^ beta)
end

;; psi-singular
; a (almost) singular psi
; radial, positive, infite at 0
to-report psi-singular [dist]
  let value 0
  carefully [
    set value 1 / (dist ^ beta)
    if (value > 10E10) [set value 10E10]
  ]
  [
    set value 10E10
  ]
  report value
end


;; Helper functions

;; atan2
; a safe atan
; no error on x=0
to-report atan2 [x y]
  let value 0
  ifelse x = 0 [
    if y < 0 [set value 180] ;otherwise 0 stays as the value
  ]
  ;else
  [
   set value atan x y
  ]
  report value
end
@#$#@#$#@
GRAPHICS-WINDOW
305
10
810
516
-1
-1
7.0
1
10
1
1
1
0
0
0
1
-35
35
-35
35
1
1
1
ticks
30.0

BUTTON
15
335
92
368
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
190
335
271
368
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
0

SLIDER
5
10
280
43
population
population
2
100
100.0
2
1
NIL
HORIZONTAL

SLIDER
15
240
187
273
beta
beta
0.1
2.5
1.1
.1
1
NIL
HORIZONTAL

SLIDER
15
277
187
310
K
K
0.1
4
0.1
.1
1
NIL
HORIZONTAL

BUTTON
100
335
181
368
go once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

SLIDER
10
180
145
213
vel-mean
vel-mean
0.1
2
1.0
.1
1
NIL
HORIZONTAL

SLIDER
150
180
285
213
vel-stdev
vel-stdev
0.1
1
0.2
.1
1
NIL
HORIZONTAL

SLIDER
50
65
222
98
radius
radius
5
30
15.0
5
1
NIL
HORIZONTAL

PLOT
830
10
1030
160
max radius
ticks
radius
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot stat-max-radius"

PLOT
1040
10
1240
160
max speed
ticks
speed
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot stat-max-speed"

MONITOR
1040
165
1240
210
max speed
stat-max-speed
17
1
11

MONITOR
830
170
1030
215
max radius
stat-max-radius
17
1
11

PLOT
825
380
1025
530
centre error
ticks
error
0.0
1.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot error-centre"

PLOT
825
545
1025
695
velocity error
ticks
error
0.0
10.0
0.0
0.001
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot error-velocity"

PLOT
1035
380
1235
530
birds at the edge of the world
ticks
birds
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count error-turtles-at-the-wall"

CHOOSER
165
130
257
175
dimension
dimension
"1D" "2D"
1

MONITOR
305
525
427
570
initial-random-seed
initial-random-seed
0
1
11

CHOOSER
195
250
290
295
psi-type
psi-type
"standard" "singular" "normalized"
2

CHOOSER
15
130
153
175
topology-type
topology-type
"random" "collision" "bi-cluster"
0

@#$#@#$#@
## WHAT IS IT?

This model is a discretization of a Cucker-Smale type model, for bird flocking. The model is expressed at the center of mass.

Each bird has a position and a velocity. The movement is expressed by a set of ordinary differential equations (ODEs), this simulation solves the equations, implicitly using a Euler method (intrinsic to the way NetLogo does the time discretization)

## HOW IT WORKS

Each bird starts with a random position and velocity.

Note that the movement we see is relative to the center of mass. We may visualize that the whole set is moving in some direction, e.g. up the screen. So that any bird that seems to be moving down may actually be moving up but slower than the average. Also when birds stop, they are not stopped at mid air, they are flying with a velocity equal to the center of mass.

Each bird tries to align with others given some communication weight, inversely related to the distance to each other. Note that we do not see birds aligning, instead as they move slower (relative to the center of mass) the velocity becomes increasingly dominated by the constant (and invisible in the representation) velocity of the center of mass.

Slow birds, i.e. birds moving very slowly relative to the center of mass (speed < 0.001), will turn their shape to a dot, because by then the (relative) direction of the movement is mostly irrelevant. If all birds become "slow" then the model stops.

For simplicity, the velocity is set in (opposing) pairs. It allows for a simple way of getting zero mean velocity, and to set collision courses. It distorts the velocity distribution, but keeps the main features (mean and variance of each one, global zero mean).


## HOW TO USE IT

Set the POPULATION size the determine the number of birds. Set the mean and standard deviation of the speed, VEL-MEAN and VEL-STDEV. Not the speed is using a beta distribution, so it has a long tail and may generate a few quite fast birds. Choose the RADIUS around the center in which the birds will be generated. Choose BETA and K for the communication function, BETA controls how strong is long range effect, the smaller the BETA the stronger, K is a overall weight. Hit SETUP and the GO, for continuous computation, or GO ONCE to go one tick at a time.

The simulation may take some time to stop, so you might want to increase the simulation speed (using the slider at the top).

The code is not robust, it does almost no error checking. If you change parameters while the model is running it is possible that it breaks. Or not... try anyway. Expect it go wrong if you do.


## THINGS TO NOTICE

Will birds flock or not? Two plots help to quantify that. MAX RADIUS is the distance from the outmost bird to the center, the actual group size is at most twice of that. MAX SPEED is the speed of the fastest bird relative to the center, the actual expansion rate of the group at most twice of that.

How does a bird's position affect its own speed? Will there be one or more flocks?

You should pay attention to error indicators, if they go too high the model is no longer a good approximation. The center of mass and also the average velocity should theoretically be fixed at zero. CENTER and VELOCITY ERROR plots show how large the drift from zero is. The velocity change is symmetric, so that error should remain near zero. 	The position of the center of mass is expected to drift while using a Euler method solution.

The world can not be infinite, so there is also an indication of how many birds have reached THE EDGE OF THE WORLD. As soon as one bird does, the model is no longer valid, and the world will become light pink, to warn about a clear loss of validity. The calculation will nevertheless carry on until there are more than 10% at the edge, because it still gives some visual idea of the interaction.


## THINGS TO TRY

How does BETA affect flocking? Analysis says BETA <= 1 (0.5 for the normalized psi) will cause flocking to occur always. For larger values it will depend on K and on the initial velocities and dispersion of the birds. So change the velocity with VEL-MEAN and VEL-STDEV and change the RADIUS. It should be easier to flock at low speeds and small radius. What happens if the speed variation is low? 

The K factor also influences how fast it takes for flocking to happen. You may check that, and tweak the simulation to get better results, by changing it. If birds are hitting the wall, increase K, if they almost do not move from the center, decrease K.

If not flocking, while the birds disperse as a whole, are there still some (sub) groups which are flocking?


## EXTENDING THE MODEL

Many many thing may be done by changing the function `psi`. Also we could add several breeds of birds, e.g creating leaders, or different species, namely predators (predators would go after prey, prey would escape)

Improve the calculation of the ODEs, e.g. extending to a Störmer-Verlet method.


## NETLOGO FEATURES

We use `links` between the birds / agents although they are not physically linked, except that they may observe each other. Yet the reporters `link-length` and `link-heading` simplify the calculations.

## RELATED MODELS

From the NetLogo library:
* Flocking
* Moths
* Flocking Vee Formation
* Flocking - Alternative Visualizations

## CREDITS AND REFERENCES

Part of a "Movimento de partículas auto-propulsionadas - Modelo de Cucker-Smale e variantes: análise e simulações", a MSc at Universidade Aberta, Portugal. The main reference for the dissertation is:

Young Pil Choi, Seung Yeal Ha e Zhuchun Li. “Emergent dynamics of the cucker–Smale flocking model and its variants”. In: Modeling and Simulation in Science, Engineering and Technology. Active Particles, Volume 1. Ed. por Nicola Bellomo, Pierre Degond e Eitan Tadmor. Birkhäuser, Cham., 2017, pp. 299–331. DOI: 10.1007/978-3-319-49996-3_8

The initial version of this was written while attending the online course on "Introduction to Agent-Based Modeling", at Santa Fe Institute's complexityexplorer.org. See: https://www.complexityexplorer.org/courses/171-introduction-to-agent-based-modeling

## HOW TO CITE

If you mention this model or the NetLogo software in a publication, we ask that you include the citations below.

For the model itself:

* Ricardo André (2023).  NetLogo Cucker-Smale type model.
"Movimento de partículas auto-propulsionadas - Modelo de Cucker-Smale e variantes: análise e simulações", MSc at Universidade Aberta, Portugal

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

## COPYRIGHT AND LICENSE

Copyright 2023 Ricardo André

![CC BY-NC-SA 4.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/4.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.
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

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.3.0
@#$#@#$#@
set population 200
setup
repeat 200 [ go ]
@#$#@#$#@
@#$#@#$#@
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
1
@#$#@#$#@
