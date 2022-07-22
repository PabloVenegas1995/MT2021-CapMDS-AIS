#mediaCalculation
from prettytable import PrettyTable
import os, re, math

solutions_quality_list = []
solutions_times_list = []
table_data = []
#mean_sol_list = []
#toTable = tuple()
#mean_time_list = []
#st_dev_list = []
row_names = [   "100, 300, t = 0.75m","100, 300, t = 0.80m","100, 300, t = 0.85m", "100, 300, Derivativ",
                "100, 600, t = 0.75m","100, 600, t = 0.80m","100, 600, t = 0.85m", "100, 600, Derivativ",
                "100, 800, t = 0.75m","100, 800, t = 0.80m","100, 800, t = 0.85m", "100, 800, Derivativ",
                "200, 300, t = 0.75m","200, 300, t = 0.80m","200, 300, t = 0.85m", "200, 300, Derivativ",
                "200, 600, t = 0.75m","200, 600, t = 0.80m","200, 600, t = 0.85m", "200, 600, Derivativ",
                "200, 800, t = 0.75m","200, 800, t = 0.80m","200, 800, t = 0.85m", "200, 800, Derivativ"]

# Specify the Column Names while initializing the Table
myTable = PrettyTable(["Configuration + Threshold", "Mean Quality of Solutions", "Standar Deviation", "Mean time"])
solutions_dir = ['T1','T2','T3','T4']

for y in range(len(solutions_dir)):

    content_dir = os.listdir(solutions_dir[y])

    for i in range(len(content_dir)):

        path = solutions_dir[y] + '/' + sorted(content_dir)[i]
        #print(path)
        with open(path, 'r') as f:
            last_line = f.readlines()[-1]
        solutions_quality_list.append((re.split(r'\t+', last_line.rstrip('\t'))[0]))
        solutions_times_list.append((re.split(r'\t+', last_line.rstrip('\n'))[1]))
        
    solutions_quality_list = list(map(int, solutions_quality_list))
    solutions_times_list = list(map(float, solutions_times_list))

    #print(solutions_quality_list)
    #print('done')
    #print(solutions_times_list)

    sum_sol = 0.0
    sum_time = 0.0
    mean_sol = 0.0
    mean_time = 0.0
    var = 0.0
    st_dev = 0.0
    start = 0

    for i in range(6):
        for j in range(10):
            sum_sol += solutions_quality_list[start+j]
            sum_time += solutions_times_list[start+j]

        mean_sol = sum_sol/10.0
        mean_time = sum_time/10.0

        for j in range(10):
            var += (solutions_quality_list[start+j] - mean_sol)**2
        
        st_dev = math.sqrt(var/10.0)
        
        table_data.append((mean_sol, round(st_dev,2), round(mean_time,2)))
        sum_sol  = 0.0
        sum_time = 0.0
        mean_sol = 0.0
        mean_time = 0.0
        var = 0.0
        st_dev = 0.0
        start+= 10
    solutions_quality_list.clear()
    solutions_times_list.clear()


# Add rows
count = 0
for k in range(6):
    for kk in range(4):
        myTable.add_row([row_names[count], table_data[k + 6*kk][0], table_data[k + 6*kk][1],table_data[k + 6*kk][2]])
        count+=1
 
print(myTable)

#para cada problema (configuración+threshold): media de solución, desviación estandar solución, media tiempo 




