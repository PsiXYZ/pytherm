import csv

d = {}
with open('2.csv') as csvfile:
    reader = csv.reader(csvfile, delimiter=';')
    for line in reader:
        # print(line)
        d[f'{line[1]+line[2]}'] = float(line[3])

for i in d:
    print(f'    \'{i}\': {d[i]},')
