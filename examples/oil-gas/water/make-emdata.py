import numpy as np

sites = np.linspace(-6000, 6000, 121)

freqs = [10, 100]

emdf = open('oil-gas.emd', 'w')

print("# frequencies", file=emdf)
print("%d" % (len(freqs)), file=emdf)
for f in freqs:
	print(" %.5E" % (f), file=emdf)

print("# transmitters", file=emdf)
print("%d" % (1), file=emdf)
print(0.0, 0.0, 3000.0, 0.0, 90.0, 1.0, 1000.0, file=emdf)

print("# receivers", file=emdf)
print("%d" % (len(sites)), file=emdf)
for i in range(0, len(sites)):
	print("% .4E % .4E % .4E" % (0.0, sites[i], 0.1), file=emdf)
print('', file=emdf)
