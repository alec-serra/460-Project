import numpy as np

t1 = np.loadtxt('tree_test.dat')
t2 = np.load('meow.npy')
print(np.shape(t1), np.shape(t2))

print(np.sqrt(np.sum(t1-t2)**2)/len(t1))
