import numpy as np

x = np.linspace(-1000, 1000, 40)
y = np.linspace(-1000, 1000, 40)

X, Y = np.meshgrid(x, y)

freqs = [10]

emdf = open('dongguashan.emd', 'w')

print("# frequencies", file=emdf)
print("%d" % (len(freqs)), file=emdf)
for f in freqs:
	print(" %.5E" % (f), file=emdf)

print("# transmitters", file=emdf)
print("%d" % (2), file=emdf)
print(0.0, 0.0, 800.0, 0.0, 90.0, 1.0, 200.0, file=emdf)
print(0.0, 0.0, 1400.0, 0.0, 90.0, 1.0, 200.0, file=emdf)

print("# receivers", file=emdf)
print("%d" % (X.shape[0]*X.shape[1]), file=emdf)
for i in range(X.shape[0]):
    for j in range(X.shape[1]):
        print("% .4E % .4E % .4E" % (X[i, j], Y[i, j], 0.1), file=emdf)
print('', file=emdf)
