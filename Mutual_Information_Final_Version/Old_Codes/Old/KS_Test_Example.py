import numpy as np
from matplotlib import pyplot as plt
import scipy.stats as st


x_val = np.linspace(0,100,1000)

norm1pdf = st.norm.pdf(x_val,scale = 10,loc = 40)
norm2pdf = st.norm.pdf(x_val,scale = 10,loc = 45)

norm1cdf = st.norm.cdf(x_val,scale = 10,loc = 40)

norm110draw = np.random.normal(50,15,100)

percentile_list = []
for i in x_val:
	percentile_list.append(st.percentileofscore(np.sort(norm110draw),i)/100)
	
norm110draw2 = np.random.normal(40,10,100)

percentile_list2 = []
for i in x_val:
	percentile_list2.append(st.percentileofscore(np.sort(norm110draw2),i)/100)

plt.plot(x_val,norm1pdf)
plt.title("Probability density function")
plt.ylabel("PDF")
plt.xlabel("x")
plt.show()
plt.plot(x_val,norm1cdf)
plt.title("cumulative distribution function")
plt.ylabel("CDF")
plt.xlabel("x")
plt.show()




kstestr = st.kstest(norm110draw,norm110draw2,alternative='two-sided')
max_diff = 0
xmax_ind = 0
bgger = "CDF"
for i in range(len(norm1cdf)):
	diff = percentile_list2[i] - percentile_list[i]
	if abs(diff) > max_diff:
		max_diff = abs(diff)
		xmax_ind = i
		if diff > 0:
			bgger = "CDF"
		else:
			bgger = "Sample"
plt.plot(x_val,percentile_list2,label = "Sample 1: mean = 40, std = 10")
plt.plot(x_val,percentile_list,label = "Sample 2: mean = 50, std = 15")
if bgger == "CDF":
	plt.vlines(x_val[xmax_ind],percentile_list[xmax_ind],percentile_list2[xmax_ind],linewidth=4,color="green",label="Max CDF Distance")
if bgger == "Sample":
	plt.vlines(x_val[xmax_ind],percentile_list2[xmax_ind],percentile_list[xmax_ind],linewidth=4,color="green",label="Max CDF Distance")
plt.title(f"100 samples from very different distributions\nMax distance between cdfs = {round(kstestr[0],3)}, p score = {round(kstestr[1],3)}")
plt.ylabel("CDF")
plt.xlabel("x")
plt.legend()
plt.show()