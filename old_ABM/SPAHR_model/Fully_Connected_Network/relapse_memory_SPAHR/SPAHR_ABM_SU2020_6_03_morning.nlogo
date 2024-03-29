;;Code for SPAR Agent Based Model
extensions [nw]

globals
[ lambda
  change
  mu_a   ; To be set in linear function mu_a_t
  alpha  ; To be set in piecewise linear function alpha_t
]

turtles-own
[
  class
  relapse-count
  current-class ; temporary variable to avoid overlapping of states in "change?"
  P_neighbor_count ; Temporary variable to count number of neighboring P's
  A_neighbor_count ; Temporary variable to count number of neighboring A's
  H_neighbor_count ; Temporary variable to count number of neighboring H's
  Addiction_History ; Last class of turtle before entering in R class
]

to setup
  ca
  ask patches [set pcolor white]
  crt turtle-count
  set-default-shape turtles "person"
  layout-circle turtles 10
  ask turtles [set color blue]
  ask turtles [set class "S"]
  ask turtles [set size 2]

  let p-count round (turtle-count * initial_P)
  let a-count round (turtle-count * initial_A)
  let r-count round (turtle-count * initial_R)
  let h-count round (turtle-count * initial_H)

  ask n-of p-count turtles with [class = "S"] [Prescribed]
  ask n-of a-count turtles with [class = "S"] [Addicted]
  ask n-of r-count turtles with [class = "S"] [Recovered]
  ask n-of h-count turtles with [class = "S"] [Heroin]

  ask turtles [ create-links-with other turtles ]
  ask links [set color black]

  reset-ticks
end

to go
  if ticks = 100 [stop]

  change?

  ask turtles with [class = "D"] [
    ask my-links [ die ]
    ask turtles with [class = "D"] [ create-links-with other turtles ]
    Susceptible
  ]
  tick
end

to change?
  ;;saves the current status of network
  ask turtles with [class = "S"][set current-class "S" set P_neighbor_count count link-neighbors with [class = "P"] set A_neighbor_count count link-neighbors with [class = "A"] set H_neighbor_count count link-neighbors with [class = "H"]]
  ask turtles with [class = "P"][set current-class "P" set P_neighbor_count count link-neighbors with [class = "P"] set A_neighbor_count count link-neighbors with [class = "A"] set H_neighbor_count count link-neighbors with [class = "H"]]
  ask turtles with [class = "A"][set current-class "A" set P_neighbor_count count link-neighbors with [class = "P"] set A_neighbor_count count link-neighbors with [class = "A"] set H_neighbor_count count link-neighbors with [class = "H"]]
  ask turtles with [class = "R"][set current-class "R" set P_neighbor_count count link-neighbors with [class = "P"] set A_neighbor_count count link-neighbors with [class = "A"] set H_neighbor_count count link-neighbors with [class = "H"]]
  ask turtles with [class = "H"][set current-class "H" set P_neighbor_count count link-neighbors with [class = "P"] set A_neighbor_count count link-neighbors with [class = "A"] set H_neighbor_count count link-neighbors with [class = "H"]]

  ;;for rehab code
  let total_A count turtles with [class = "A"]
  let total_H count turtles with [class = "H"]

  alpha_t ; Updates piecewise linear or constant value of alpha
  mu_a_t  ; Updates linear or constant value of mu_a

  ;;gets lambda
  ; Changes here.........
  ask turtles [

    if (current-class = "S") [ ;;P and A are the link neighbors with class P and A, so beta_A * A is beta_A * link neigh. with class A.
      set lambda (alpha + beta_A * (A_neighbor_count / turtle-count) + beta_P * (P_neighbor_count / turtle-count) + theta_1 * (H_neighbor_count / turtle-count) + mu) ; Note: alpha is alpha(t)

      let prob (1 - exp (- lambda * delta_t))
      let rf random-float 1

      let s_to_p (alpha / lambda)
      let s_to_a ((beta_A * (A_neighbor_count / turtle-count) + beta_P * (P_neighbor_count / turtle-count) ) / lambda)
      let s_to_h (theta_1 * (H_neighbor_count / turtle-count) / lambda)
      let d ( mu / lambda )

      let rf2 random-float 1

      if rf < prob [

        ;;changes S -> P
        if s_to_p > rf2 [
          Prescribed
          ;set s-p (s-p + 1)
        ]

        ;;changes S -> A
        if s_to_a + s_to_p > rf2 and rf2 >= s_to_p [
          Addicted
          ;set s-a (s-a + 1)
        ]

        if (s_to_h + s_to_a + s_to_p) > rf2 and rf2 >= (s_to_a + s_to_p) [
          Heroin
          ;set s-h (s-h + 1)
        ]

        ;;Deads them
        if (s_to_h + s_to_a + s_to_p + d) > rf2 and rf2 >= (s_to_h + s_to_a + s_to_p) [
          Dead
        ]

      ]
    ]

    if (current-class = "P") [

      set lambda (epsilon + gamma + theta_2 * (H_neighbor_count / turtle-count) + mu)

      let prob (1 - exp (- lambda * delta_t))

      let rf random-float 1

      let p_to_s (epsilon / lambda)
      let p_to_a (gamma / lambda)
      let p_to_h (theta_2 * (H_neighbor_count / turtle-count) / lambda)
      let d (mu / lambda)

      let rf2 random-float 1

      if rf < prob [
        ;;changes P -> S
        if p_to_s > rf2 [
          Susceptible
          ;set p-s (p-s + 1)
        ]

        ;;changes P -> A
        if p_to_a + p_to_s > rf2 and rf2 >= p_to_s [ Addicted ]

        ;; changes P -> H
        if (p_to_h + p_to_a + p_to_s) > rf2 and rf2 >= p_to_a + p_to_s [ Heroin ]

        ;; Deads them
        if (p_to_h + p_to_a + p_to_s + d) > rf2 and rf2 >= (p_to_h + p_to_a + p_to_s) [ Dead ]

      ]
    ]

    if (current-class = "A") [

      set lambda (zeta + theta_3 * (H_neighbor_count / turtle-count) + (mu + mu_a) ) ; Note: mu_a is mu_a(t)

      let prob (1 - exp (- lambda * delta_t))

      let rf random-float 1

      let a_to_r (zeta / lambda)
      let a_to_h (theta_3 * (H_neighbor_count / turtle-count) / lambda)
      let d ( (mu + mu_a) / lambda )

      let rf2 random-float 1

      if rf < prob [
        ;;changes A -> R
        if a_to_r > rf2 [ Recovered set Addiction_History "A" ]

        ;;changes A -> H
        if (a_to_r + a_to_h) > rf2 and rf2 >= a_to_r [ Heroin ]

        ;;Deads them
        if (a_to_r + a_to_h + d) > rf2 and rf2 >= (a_to_r + a_to_h) [ Dead ]

      ]
    ]

    ; Recovered case
    if current-class = "R" [ ;; keep global not neighbors...

      set lambda  ( ( (sigma * total_A) / (total_A + total_H + 0.0001) ) + ( (sigma * total_H) / (total_A + total_H + 0.001) ) + mu ) ;; note the 0.0001 is meant to guarantee a non zero denom.

      let prob (1 - exp (- lambda * delta_t))
      let rf random-float 1

      ;let r_to_a ( ( (sigma * total_A) / (total_A + total_H + 0.0001) ) / lambda)
      ;let r_to_h ( ( (sigma * total_H) / (total_A + total_H + 0.0001) ) / lambda)
      let r_to_aorh (sigma / lambda) ;to be used when taking history into account
      let d (mu / lambda)

      let rf2 random-float 1

      if rf < prob [
        ;;changes R -> A
        ;if r_to_a > rf2 [ Addicted ]

        ;;changes R -> H
        ;if r_to_h + r_to_a > rf2 and rf2 >= r_to_a [ Heroin ]

        ;;Deads them
        ;if (d + r_to_h + r_to_a) > rf2 and rf2 >= r_to_h + r_to_a [ Dead ]



;;;;TO Be implemented after behavior space. History is next steps. to replace ln 171 -178 ish.
        ;;changes R -> A or H
        if r_to_aorh > rf2 [ relapse ]

        ;;Deads them
        if (d + r_to_aorh ) > rf2 and rf2 >= r_to_aorh [ Dead ]


        ;does current-class change my pens? or no since the current-class and GUI match up at the end of a tick?
        ;does H procedure below belong under R chek or now? what should be our logic now. - Q for meeting with owen.
      ]
    ]

    if current-class = "H" [

      set lambda ( v_nu + mu + mu_h )

      let prob (1 - exp (- lambda * delta_t))
      let rf random-float 1

      let h_to_r (v_nu / lambda)
      let d ( (mu + mu_h) / lambda )

      let rf2 random-float 1

      if rf < prob [
       ;;changes H -> R
        if rf2 < h_to_r [ Recovered set Addiction_History "H"]
       ;;Deads them
        if h_to_r + d > rf2 and rf2 > h_to_r [ Dead ]
      ]

    ]
  ]

end

; Implements the piecewise linear alpha(t) function
to alpha_t

  let t (ticks * delta_t)

  ;If alpha is constant
  if (constant_alpha = 1) [ set alpha alpha_c ]

  ;If alpha is non-constant
  if (constant_alpha != 1) [
    if (t <= 3.25) [
      set alpha ( m_tilde * t + b_tilde )
    ]

    ; Equivalent to an "else" statment
    if (t > 3.25) [
      set alpha ( m_tilde * (3.25) + b_tilde + c_tilde * (t - 3.25) )
    ]
  ]

end

; Implements the linear mu_a(t) function
to mu_a_t

  let t (ticks * delta_t)

  ;If mu_a is constant
  if (constant_mu_a = 1) [ set mu_a mu_a_c ]

  ;If mu_a is non-constant
  if (constant_mu_a != 1) [
    set mu_a ( d_tilde * t + e_tilde )
  ]
end

to Susceptible
  set class "S"
  set color blue
end

to Prescribed
  set class "P"
  set color yellow
end

to Addicted
  set class "A"
  set color red
end

to Recovered
  set class "R"
  set color green
end

to relapse
  if Addiction_History = "H" [Heroin]
  if Addiction_History = "A" [Addicted]
end

to Dead
  set class "D"
  set color black
  set Addiction_History "None"
end

to Heroin
  set class "H"
  set color violet
end
@#$#@#$#@
GRAPHICS-WINDOW
210
10
647
448
-1
-1
13.0
1
10
1
1
1
0
0
0
1
-16
16
-16
16
1
1
1
ticks
30.0

INPUTBOX
25
32
180
92
turtle-count
500.0
1
0
Number

BUTTON
67
101
136
134
Set Up
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
25
150
88
183
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

BUTTON
114
150
189
183
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
1

INPUTBOX
447
501
557
561
initial_P
0.1
1
0
Number

INPUTBOX
651
500
760
560
initial_A
0.016
1
0
Number

PLOT
661
12
1005
283
Group Proportions 
Ticks
Class Proportion
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"Suseptible" 1.0 0 -10649926 true "" "plot ( count turtles with [ class = \"S\" ] / turtle-count )"
"Prescribed" 1.0 0 -1184463 true "" "plot ( count turtles with [ class = \"P\" ] / turtle-count )"
"Addicted" 1.0 0 -2139308 true "" "plot ( count turtles with [ class = \"A\" ] / turtle-count )"
"Heroin" 1.0 0 -6917194 true "" "plot ( count turtles with [ class = \"H\" ] / turtle-count )"
"Recovered" 1.0 0 -13840069 true "" "plot ( count turtles with [ class = \"R\" ] / turtle-count )"

INPUTBOX
1250
511
1360
571
initial_R
0.014
1
0
Number

SLIDER
415
605
587
638
epsilon
epsilon
0
100
2.53
1
1
NIL
HORIZONTAL

SLIDER
415
573
587
606
gamma
gamma
0
100
0.4
1
1
NIL
HORIZONTAL

SLIDER
616
573
788
606
zeta
zeta
0
100
0.198
1
1
NIL
HORIZONTAL

SLIDER
1219
580
1391
613
sigma
sigma
0
100
0.102
1
1
NIL
HORIZONTAL

SLIDER
19
229
191
262
mu
mu
0
100
0.0071
1
1
NIL
HORIZONTAL

SLIDER
20
603
192
636
beta_P
beta_P
0
100
0.00654
1
1
NIL
HORIZONTAL

SLIDER
20
570
192
603
beta_A
beta_A
0
100
8.78E-4
1
1
NIL
HORIZONTAL

SLIDER
19
262
191
295
delta_t
delta_t
0.01
1
0.1
0.01
1
NIL
HORIZONTAL

MONITOR
673
299
730
344
S
count turtles with [class = \"S\"]
17
1
11

MONITOR
739
299
796
344
P
count turtles with [class = \"P\"]
17
1
11

MONITOR
804
299
861
344
A
count turtles with [class = \"A\"]
17
1
11

MONITOR
939
299
996
344
R
count turtles with [class = \"R\"]
17
1
11

INPUTBOX
51
498
160
558
initial_S
0.85
1
0
Number

SLIDER
20
636
192
669
theta_1
theta_1
0
1
0.222
0.01
1
NIL
HORIZONTAL

SLIDER
415
638
587
671
theta_2
theta_2
0
100
0.236
1
1
NIL
HORIZONTAL

SLIDER
616
606
788
639
theta_3
theta_3
0
100
5.0
1
1
NIL
HORIZONTAL

SLIDER
1021
579
1193
612
v_nu
v_nu
0
100
5.31E-4
1
1
NIL
HORIZONTAL

SLIDER
1021
612
1193
645
mu_h
mu_h
0
100
0.0466
1
1
NIL
HORIZONTAL

TEXTBOX
1060
472
1164
490
Heroin Parameters
11
0.0
1

INPUTBOX
1052
510
1162
570
initial_H
0.02
1
0
Number

MONITOR
871
299
928
344
H
count turtles with [class = \"H\"]
17
1
11

SLIDER
216
571
388
604
m_tilde
m_tilde
0
100
0.0
1
1
NIL
HORIZONTAL

SLIDER
216
604
388
637
b_tilde
b_tilde
0
100
0.0
1
1
NIL
HORIZONTAL

SLIDER
216
637
388
670
c_tilde
c_tilde
0
100
0.0
1
1
NIL
HORIZONTAL

SLIDER
816
573
988
606
d_tilde
d_tilde
0
100
0.0
1
1
NIL
HORIZONTAL

SLIDER
816
606
988
639
e_tilde
e_tilde
0
100
1.0
1
1
NIL
HORIZONTAL

TEXTBOX
221
497
390
525
alpha(t) = prescription rate per person per year\n
11
0.0
1

TEXTBOX
819
494
983
521
mu_a(t) = overdose death rate of addicted users
11
0.0
1

SLIDER
217
716
389
749
alpha_c
alpha_c
0
100
0.4
1
1
NIL
HORIZONTAL

SLIDER
817
689
989
722
mu_a_c
mu_a_c
0
100
0.00883
1
1
NIL
HORIZONTAL

TEXTBOX
141
470
291
488
Susceptible Parameters
11
0.0
1

SLIDER
243
676
359
709
constant_alpha
constant_alpha
0
1
1.0
1
1
NIL
HORIZONTAL

SLIDER
843
647
960
680
constant_mu_a
constant_mu_a
0
1
1.0
1
1
NIL
HORIZONTAL

TEXTBOX
444
474
567
492
Prescribed Parameters
11
0.0
1

TEXTBOX
747
468
897
486
Addicted Parameters
11
0.0
1

TEXTBOX
1245
474
1395
492
Recovered Parameters
11
0.0
1

TEXTBOX
1224
631
1374
649
sigma = natural relapse
11
0.0
1

TEXTBOX
53
205
203
223
Overall Parameters
11
0.0
1

TEXTBOX
418
679
590
721
gamma = prescription-induced addiction
11
0.0
1

TEXTBOX
418
718
593
746
epsilon = end prescription w/o addiction
11
0.0
1

TEXTBOX
622
647
790
675
zeta = rate of A entry into rehab
11
0.0
1

TEXTBOX
1028
652
1178
680
v_nu = rate of H entry into rehab
11
0.0
1

TEXTBOX
28
680
178
698
beta_A = illicit addiction
11
0.0
1

TEXTBOX
28
706
178
724
beta_P = illicit addiction
11
0.0
1

TEXTBOX
1027
687
1193
729
mu_h = overdose death rate for heroin/fentanyl users
11
0.0
1

TEXTBOX
29
731
179
759
theta_1 = illicit heroin addiction
11
0.0
1

TEXTBOX
622
686
781
714
theta_3 = rate of A addiction to heroin
11
0.0
1

TEXTBOX
209
533
413
575
m_tilde, b_tilde, c_tilde := parameters in piecewise linear alpha(t)
11
0.0
1

TEXTBOX
819
529
991
571
d_tilde, e_tilde := parameters in linear mu_a(t)
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
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
NetLogo 6.1.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="heroin_run" repetitions="300" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="100"/>
    <metric>turtle-count</metric>
    <metric>mu</metric>
    <metric>delta_t</metric>
    <metric>initial_S</metric>
    <metric>beta_A</metric>
    <metric>beta_P</metric>
    <metric>theta_1</metric>
    <metric>constant_alpha</metric>
    <metric>alpha_c</metric>
    <metric>m_tilde</metric>
    <metric>b_tilde</metric>
    <metric>c_tilde</metric>
    <metric>initial_P</metric>
    <metric>gamma</metric>
    <metric>epsilon</metric>
    <metric>theta_2</metric>
    <metric>initial_A</metric>
    <metric>zeta</metric>
    <metric>theta_3</metric>
    <metric>constant_mu_a</metric>
    <metric>mu_a_c</metric>
    <metric>d_tilde</metric>
    <metric>e_tilde</metric>
    <metric>initial_H</metric>
    <metric>v_nu</metric>
    <metric>mu_h</metric>
    <metric>initial_R</metric>
    <metric>sigma</metric>
    <metric>count turtles with [class = "S"] / turtle-count</metric>
    <metric>count turtles with [class = "P"] / turtle-count</metric>
    <metric>count turtles with [class = "A"] / turtle-count</metric>
    <metric>count turtles with [class = "R"] / turtle-count</metric>
    <metric>count turtles with [class = "H"] / turtle-count</metric>
  </experiment>
  <experiment name="heroin_run_test" repetitions="300" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="100"/>
    <metric>count turtles with [class = "S"] / turtle-count</metric>
    <metric>count turtles with [class = "P"] / turtle-count</metric>
    <metric>count turtles with [class = "A"] / turtle-count</metric>
    <metric>count turtles with [class = "H"] / turtle-count</metric>
    <metric>count turtles with [class = "R"] / turtle-count</metric>
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
