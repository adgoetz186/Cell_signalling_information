import numpy as np
import seaborn as sns
import matplotlib.ticker as ticker
from matplotlib import pyplot as plt

mtrx = np.random.rand(55,55)
ax = sns.heatmap(mtrx)

#ticker_to_use.tick_values(10,11)
#plt.xticks([0,10,34])


ax.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(3))
plt.xticks([0,,40],labels=["a","b","c"])
plt.title("Information")
plt.xlabel("$log_{10}k_{on}s^{-1}M^{-1}$")
plt.show()