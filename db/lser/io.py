import csv
import os

path = os.path.abspath(__file__)[:-5:]


def read_subs():
    subs = []
    with open(path + '/subs.csv',  'r') as csvfile:
        data = csv.reader(csvfile)
        for line in data:
            subs.append(line)
    return subs[1:]

def write_subs(data):
    with open(path + '/out.csv', 'w') as outfile:
        wr = csv.writer(outfile)
        wr.writerow(['Formula', 'E', 'S', 'A', 'B', 'V', 'L'])
        wr.writerows(data)
