# disclaimer

The "@time" macro accounts for parsing, compiling and running time. To only consider run times you must use "@btime", from the package BenchmarkTools.


# comparisons with previous python script

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


# calculate attenuation coefficientes just once?

with "if only_seis", always calculating coefs
321 321 5000 46 s
with "if only_seis", calculating coefs only once
321 321 5000 31 s
321 321 2500 15.7 s
without "ifs", calculating coefs only once
321 321 5000 32 s (wow)


# benchmarks for old, purer functions

Propagation times (refference) (!compilation)
(321,321,3000) pure:                5.5 s
(321,321,3000) 0 tpr, !att.:        6.0 s   -> 8% slwr than pure
(400,400,3000) pure:                8.5 s
(321,321,3000) 60 tpr, !att.:       9.5 s   -> 11% slwr than pure
(321,321,3000) 60 tpr, att.:        12.3 s  -> 29% slwr due to att


# benchmarks for functions before I knew metaprogramming

new w. 60 tpr
(321,321,3000) !saving:             12.4 s -> turned to 17 s one day before :o
(321,321,3000) saving seis to SSD:  12.5 s
(321,321,3000) saving seis to HDD:  12,5 s
(321,321,3000) saving snaps to SSD: 18,8 s (21 s before)
(321,321,3000) saving snaps to HDD: 19,8 s

So I'm using my HDD, as there is no difference in time wasted.


# benchmarks for functions with metaprogramming

propagate with metaprogramming
(321,321,3000) !saving:             11.7 s (16,5 s)
(321,321,3000) saving snaps to HDD: 19,9 s


# 20200822 update: rewriting propagate without sectors

(321,321,3000) !saving ORDER=2  6,7 s
(321,321,3000) !saving ORDER=6  7,3 s


# 20200823 update: with one only spatial loop

(321,321,3000) !saving ORDER=2  6,2 s
(321,321,3000) !saving ORDER=4  
(321,321,3000) !saving ORDER=6  
(321,321,3000) !saving ORDER=8  
