import numpy as np

sites = np.linspace(-500, 500, 26)

freqs = [100, 1000]

emdf = open('borehole.emd', 'w')

print("# frequencies", file=emdf)
print("%d" % (len(freqs)), file=emdf)
for f in freqs:
	print(" %.5E" % (f), file=emdf)

print("# transmitters", file=emdf)
print("%d" % (1), file=emdf)
print(0.0, 0.0, 100.0, 0.0, 90.0, 1.0, 100.0, file=emdf)

print("# receivers", file=emdf)
print("%d" % (len(sites)), file=emdf)
for i in range(0, len(sites)):
	print("% .4E % .4E % .4E" % (sites[i], 0.0, 0.1), file=emdf)
print('', file=emdf)
