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
