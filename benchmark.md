Without dealing with absorption layers


  100x100x2500
    python compiled version: 5.35 s  lol tried 3 times
    julia compiled version:  0.41 s

  200x200x2500
    python compiled version: 1.85 s
    julia compiled version:  1.52 s

  200x200x1000
    python compiled version: 0.72 s
    julia compiled version:  0.62 s

  400x400x2500
    python compiled version: 5.7 s
    julia compiled version:  6.5 s

  600x600x1000
    python compiled version: 5.4 s
    julia compiled version:  5.4 s

Observed:
  Extra allocated memory is a function of time


with "if only_seis", always calculating coefs
321 321 5000 46 s
with "if only_seis", calculating coefs only once
321 321 5000 31 s
321 321 2500 15.7 s
without "ifs", calculating coefs only once
321 321 5000 32 s (wow)


Propagation times (refference) (!compilation)
(321,321,3000) pure:                5.5 s
(321,321,3000) 0 tpr, !att.:        6.0 s   -> 8% slwr than pure
(400,400,3000) pure:                8.5 s
(321,321,3000) 60 tpr, !att.:       9.5 s   -> 11% slwr than pure
(321,321,3000) 60 tpr, att.:        12.3 s  -> 29% slwr due to att

new w. 60 tpr
(321,321,3000) !saving:             12.4 s -> turned to 17 s one day before :o
(321,321,3000) saving seis to SSD:  12.5 s
(321,321,3000) saving seis to HDD:  12,5 s
(321,321,3000) saving snaps to SSD: 18,8 s (21 s before)
(321,321,3000) saving snaps to HDD: 19,8 s

So I'm using my HDD :S
