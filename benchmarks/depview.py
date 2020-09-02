import numpy as np
import matplotlib.pyplot as plt


x = [2, 4, 6, 8]
y = [6.378, 6.822, 7.632, 7.963]

a, b = np.polyfit(x, y, 1)
eq = f'y_ls(x) = {b:.2g} + {a:.2g}x'

f = np.poly1d((a, b))

plt.plot(x, y, label="y")
plt.plot(x, f(x), label=eq, alpha=.8)
plt.legend()
plt.show()


print(eq)
print(f'{f(16)=}')
