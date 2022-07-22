# importing the required modules
import matplotlib.pyplot as plt
import numpy as np

x = []
y = []
count = 0

for line in open('toPlot.txt', 'r'):
    lines = [i for i in line.split()]
    x.append(count)
    y.append(int(lines[0]))
    count = count+1
	
plt.title("Best per Threshold")
plt.xlabel('Generation')
plt.ylabel('Fitness')
plt.yticks(y)
plt.plot(x, y, marker = 'o', c = 'g')

plt.show()
