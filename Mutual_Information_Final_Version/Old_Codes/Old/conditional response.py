import numpy as np
from scipy.stats import nbinom
from scipy.stats import gamma
from matplotlib import pyplot as plt


def mutual_information_from_matrix(relative_input_vector, conditional_probability_matrix):
    #Takes in the relative input vector (1D Numpy array) and the systems conditional probablility matrix (2D Numpy array)
    #Returns the mutual information between the input signal and the response
    signal_ratios = [relative_input_vector]
    response_values = np.matmul(signal_ratios, conditional_probability_matrix)
    signal_matrix = np.transpose(np.repeat(signal_ratios, repeats=np.shape(conditional_probability_matrix)[1], axis=0))
    response_matrix = np.repeat(response_values, repeats=np.shape(conditional_probability_matrix)[0], axis=0)
    ratio = np.divide(conditional_probability_matrix, response_matrix, out=np.zeros_like(conditional_probability_matrix), where=response_matrix != 0)
    log_of_ratio = np.log2(ratio, out=np.zeros_like(ratio), where=ratio != 0)
    mutual_info_matrix = np.multiply(np.multiply(signal_matrix, conditional_probability_matrix), log_of_ratio)
    return np.sum(np.nan_to_num(mutual_info_matrix, nan=0))

domaincont = np.arange(0,100,0.1)
domain = np.arange(0,100)
mean = 0.1
variance = 45
scale0 = variance / mean
shape0 = mean / scale0
print(f"Low signal: Shape {shape0}, Scale {scale0}")
distribution = gamma.pdf(domaincont, shape0, scale=scale0)

mean = 0.5
variance = 55
scale1 = variance / mean
shape1 = mean / scale1
print(f"High signal: Shape {shape1}, Scale {scale1}")
distribution1 = gamma.pdf(domaincont, shape1, scale=scale1)






mean = 0.1
variance = 45
p = mean / variance
n = mean * p / (1 - p)
distribution2 = nbinom.pmf(domain, n, p)


mean = 0.5
variance = 55
p = mean / variance
n = mean * p / (1 - p)
distribution3 = nbinom.pmf(domain, n, p)
mutual_info = 0.721179
mutual_info_flooded = 0.274641
plt.plot(domaincont,distribution,alpha=0.5,label = "Low Signal")
plt.plot(domaincont,distribution1,alpha=0.5, label = "High Signal")
plt.title(f"Gamma Distributions\nInfinite Density Case\nMutual Information: {round(mutual_info_flooded,5)}")
plt.xlabel("Response Value")
plt.ylabel("Probability Density Of Response")
plt.legend()
plt.show()
signal = np.array([.5,.5])
crm = np.vstack((distribution2,distribution3))
r = np.matmul(signal,crm)
print(r)
mutual_info = mutual_information_from_matrix(signal,crm)
print(mutual_info)
plt.bar(domain,distribution2,width=1,alpha=0.5,label = "Low Signal")
plt.bar(domain,distribution3,width=1,alpha=0.5, label = "High Signal")
plt.title(f"Negative Binomial Distributions\nMutual Information: {round(mutual_info,5)}")
plt.xlabel("Response Value")
plt.ylabel("Probability Of Response")
plt.legend()
plt.show()



bin_size = 10
domaincont = np.arange(0,100,bin_size)
smoothdomain = np.arange(0,100,0.1)
domain = np.arange(0,100)
mean = 0.1
variance = 45
scale0 = variance / mean
shape0 = mean / scale0
print(f"Low signal: Shape {shape0}, Scale {scale0}")
distributionpre0 = gamma.pdf(domaincont, shape0, scale=scale0)
distribution0 = distributionpre0/np.sum(distributionpre0)

mean = 0.5
variance = 55
scale1 = variance / mean
shape1 = mean / scale1
print(f"High signal: Shape {shape1}, Scale {scale1}")
distributionpre1 = gamma.pdf(domaincont, shape1, scale=scale1)
distribution1 = distributionpre1/np.sum(distributionpre1)


alt_distribution_1 = gamma.cdf(domaincont, shape0, scale=scale0)
alt_distribution_1 = alt_distribution_1[1:]-alt_distribution_1[0:-1]
alt_distribution_1 = alt_distribution_1/np.sum(alt_distribution_1)

alt_distribution_2 = gamma.cdf(domaincont, shape1, scale=scale1)
alt_distribution_2 = alt_distribution_2[1:]-alt_distribution_2[0:-1]
alt_distribution_2 = alt_distribution_2/np.sum(alt_distribution_2)





signal = np.array([.5,.5])
crm = np.vstack((distribution0,distribution1))
mutual_info = mutual_information_from_matrix(signal,crm)

alt_crm = np.vstack((alt_distribution_1,alt_distribution_2))
alt_mutual_info = mutual_information_from_matrix(signal,alt_crm)

print(mutual_info)
print(alt_mutual_info)



domaincontshift = (domaincont+bin_size/2)[1:]
print(domaincontshift)
print(alt_distribution_1)
print(alt_distribution_2)
plt.bar(domaincontshift,alt_distribution_1,width=bin_size/bin_size,alpha=0.5,label = "Low Signal",color = "C0")
plt.bar(domaincontshift,alt_distribution_2,width=bin_size/bin_size,alpha=0.5, label = "High Signal",color = "C1")
plt.title(f"Integrative Method\nBin Size: {bin_size}\nMutual Information: {round(alt_mutual_info,5)}")
plt.xlabel("Response Value")
plt.ylabel("Probability Of Response")
plt.show()



input()



plt.plot(smoothdomain,gamma.pdf(smoothdomain, shape0, scale=scale0))
plt.plot(smoothdomain,gamma.pdf(smoothdomain, shape1, scale=scale1))
plt.title(f"Step 1: Use Moments to Obtain Gamma Distributions")
plt.xlabel("Response Value")
plt.ylabel("Probability Density Of Response")
plt.show()


plt.plot(smoothdomain,gamma.pdf(smoothdomain, shape0, scale=scale0))
plt.plot(smoothdomain,gamma.pdf(smoothdomain, shape1, scale=scale1))
plt.scatter(domaincont,distributionpre0)
plt.scatter(domaincont,distributionpre1)
plt.title(f"Step 2: Evaluate Gamma Distribution\nAt Several Points Equal Distances Apart")
plt.xlabel("Response Value")
plt.ylabel("Probability Density Of Response")
plt.show()


plt.plot(smoothdomain,gamma.pdf(smoothdomain, shape0, scale=scale0))
plt.plot(smoothdomain,gamma.pdf(smoothdomain, shape1, scale=scale1))
plt.scatter(domaincont,distributionpre0)
plt.scatter(domaincont,distributionpre1)
plt.bar(domaincont,distributionpre0,width=bin_size/bin_size,alpha=0.5,label = "Low Signal",color = "C0")
plt.bar(domaincont,distributionpre1,width=bin_size/bin_size,alpha=0.5, label = "High Signal",color = "C1")
plt.title(f"Step 3: Use these points to define a discrete distribution")
plt.xlabel("Response Value")
plt.ylabel("Probability Density Of Response")
plt.show()

plt.bar(domaincont,distribution0,width=bin_size/bin_size,alpha=0.5,label = "Low Signal",color = "C0")
plt.bar(domaincont,distribution1,width=bin_size/bin_size,alpha=0.5, label = "High Signal",color = "C1")
plt.title(f"Step 4: Normalize The Distribution")
plt.xlabel("Response Value")
plt.ylabel("Probability Of Response")
plt.show()

plt.bar(domaincont,distribution0,width=bin_size/bin_size,alpha=0.5,label = "Low Signal",color = "blue")
plt.bar(domaincont,distribution1,width=bin_size/bin_size,alpha=0.5, label = "High Signal",color = "orange")
plt.show()


plt.bar(domaincont,distribution0,width=bin_size/bin_size/10,alpha=0.5,label = "Low Signal")
plt.bar(domaincont,distribution1,width=bin_size/bin_size/10,alpha=0.5, label = "High Signal")
plt.title(f"PDF Evaluation Approach\nBin Size: {bin_size}\nMutual Information: {round(mutual_info,5)}")
plt.xlabel("Response Value")
plt.ylabel("Probability Of Response")
plt.legend()
plt.show()
