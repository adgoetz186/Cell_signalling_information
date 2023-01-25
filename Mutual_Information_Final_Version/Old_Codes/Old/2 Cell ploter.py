import matplotlib.pyplot as plt
listtoplot1 = [[10, 2.6089723661964332], [11, 2.5952157897436354], [12, 2.5267519258361713], [13, 2.4574710602373933], [14, 2.399632227158036], [15, 2.3555087752620922], [16, 2.3245237968060164], [17, 2.298108717299646], [18, 2.2807998782244163], [19, 2.287301386609937], [20, 2.2805190999314773]]
listplot1 = []
x = []
for i in listtoplot1:
    listplot1.append(i[1])
    x.append(i[0])
plt.plot(x,listplot1,color = "red",label = "Old Method For Population MI Calculation")


listtoplot2 = [[10, 2.6041006922238408], [11, 2.6758743446477093], [12, 2.735408561277724], [13, 2.7860361559976585], [14, 2.8267224558980772], [15, 2.8768121507816558], [16, 2.9153202007365446], [17, 2.970149184618223], [18, 3.006838276464214], [19, 3.041763598408023], [20, 3.0668375450611967]]
listplot2 = []
for i in listtoplot2:
    listplot2.append(i[1])
plt.plot(x,listplot2,color = "orange",label = "Cell B's MI")

listplot4 = []
for i in listplot1:
    listplot4.append(2.6041006922238408)
plt.plot(x,listplot4,color = "blue",label = "Cell A's MI")

plt.legend()
plt.xlabel("$\\theta$ for Cell B")
plt.ylabel("Mutual Information")

plt.show()


listtoplot3 = [[10, 2.607385938327967], [11, 2.640991165943947], [12, 2.6706695558819566], [13, 2.691816147217573], [14, 2.718743004580422], [15, 2.744753376824464], [16, 2.7707675147197524], [17, 2.7886258838342086], [18, 2.8057455263069984], [19, 2.8239890234213627], [20, 2.842183499830375]]
listplot3 = []
for i in listtoplot3:
    listplot3.append(i[1])
plt.plot(x,listplot3,color = "green",label = "New Method For Population MI Calculation")

listplot4 = []
for i in listplot1:
    listplot4.append(2.6041006922238408)
plt.plot(x,listplot4,color = "blue",label = "Cell A's MI")

listtoplot2 = [[10, 2.6041006922238408], [11, 2.6758743446477093], [12, 2.735408561277724], [13, 2.7860361559976585], [14, 2.8267224558980772], [15, 2.8768121507816558], [16, 2.9153202007365446], [17, 2.970149184618223], [18, 3.006838276464214], [19, 3.041763598408023], [20, 3.0668375450611967]]
listplot2 = []
for i in listtoplot2:
    listplot2.append(i[1])
plt.plot(x,listplot2,color = "orange",label = "Cell B's MI")

plt.legend()
plt.xlabel("$\\theta$ for Cell B")
plt.ylabel("Mutual Information")

plt.show()