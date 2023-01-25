import testwhat

GAN_order = [i+1 for i in range(72)]
testwhat.shuffle(GAN_order)
Rank = [[i+1,j+1] for i in range(8) for j in range(9)]
print(Rank)
count = [0,0,0]
for i in range(len(GAN_order)):
    if GAN_order[i] <25:
        Rank[i].append("GAN")
        count[0] =1 + count[0]
    elif GAN_order[i] < 49:
        Rank[i].append("VAE")
        count[1] = 1 + count[1]
    else:
        Rank[i].append("Real")
        count[2] = 1 + count[2]
print(GAN_order)
print(Rank)
print(" df ")
for i in range(len(Rank)):
    if i % 6 == 5:
        print(Rank[i])
    else:
        print(Rank[i], end="")
print(count)